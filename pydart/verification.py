import numpy as np
from sys import exit
from scipy import signal as sg
from scipy.signal import convolve2d

def contingency(fcst,obs,threshold):
   """
   Calculate contingency table values.
   obs and fcst have dimensions of [ny,nx]
   threshold is a double

   Return:
      a = correct hit
      b = False Alarm
      c = Miss
      d = correct negative
   """
   a = float(np.sum(np.where((fcst > threshold) & (obs > threshold),1,0)))    # Correct Hit
   b = float(np.sum(np.where((fcst > threshold) & (obs < threshold),1,0)))    # False Alarm
   c = float(np.sum(np.where((fcst < threshold) & (obs > threshold),1,0)))    # Miss
   d = float(np.sum(np.where((fcst < threshold) & (obs < threshold),1,0)))    # Correct Negative
   return a,b,c,d

def define_kernel(radius, dropoff=0.0):
   """
   Create a convolution kernel based upon provided radius

   Required Inputs
     dropoff  --  value between 0 and 1 to reduce influence of points near edge of neighborhood
                  (values closer to 1 will result in greater dropoff with range).
   Returns:
      kernel - The outputted kernel
   """
   if radius < 0:
      exit('ERROR: Invalid neighborhood radius -- radius cannot be negative!')

   kernel = np.zeros((int((2*np.floor(radius))+1), int(2*(np.floor(radius))+1    )))
   for i in range(0,kernel.shape[0]):
      for j in range(0,kernel.shape[1]):
         icoord = abs(i - (np.floor(radius)))
         jcoord = abs(j - (np.floor(radius)))
         if ((icoord*icoord) + (jcoord*jcoord) <= (radius*radius)):
            kernel[i, j] = 1.0 - ( ((icoord*icoord) + (jcoord*jcoord)) / (radius*radius) ) * dropoff
   return(kernel)

def gen_buffer(field,kernel_radius,threshold):
   """

   Requires: field, kernel_radius,threshold

   Returns: conv - The convolution array that

   """
   kernel = define_kernel(kernel_radius)
   conv = sg.convolve2d(field > threshold, kernel, mode='same', boundary='symm')
   conv[conv>0] = threshold+1.
   return conv

def performance_scores(fcst,obs,threshold,kernel=None):

   if kernel is not None:
     fcst = gen_buffer(fcst,kernel,threshold)
     obs = gen_buffer(obs,kernel,threshold)
   
   a,b,c,d = contingency(fcst,obs,threshold)
   FOH  = float(a)/float(a+b)
   POD  = float(a)/float(a+c)
   CSI  = float(a)/float(a+b+c)
   BIAS = float(a+b)/float(a+c)
   return FOH,POD,CSI,BIAS
