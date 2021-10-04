import numpy as np
from sys import exit
from scipy import signal as sg
from scipy.signal import convolve2d
from scipy.ndimage import gaussian_filter
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

def gen_buffer(field,kernel_radius,threshold,ones=False):
   """
   Create a buffer around with a specified kerenel radius around all gridpoints
   that exceed a threshold.

   Requires: 
      field:          The field buffer is applied to [ny,nx]
      kernel_radius:  The influence radius (float)
      threshold:      The verified field

   Optional:
      ones:  A bool values.  If true the return array will just be 1/0's 
                             1s suggest the threshold is hit (important for 
                             probabilistic forecasts)

                             If false, the array will be 0's and the threshold + 1 value

   Returns: conv - The convolution array that

   """
   kernel = define_kernel(kernel_radius)
   conv = sg.convolve2d(field > threshold, kernel, mode='same', boundary='symm')
   if ones:
      conv[conv>0] = 1.
   else:
      conv[conv>0] = threshold+1.
   return conv

def performance_scores(fcst,obs,threshold,kernel=None):
   """
   Get the scores used to plot a performance diagram 
 
   Required inputs:
         fcst:  Predicted field [ny,nx]
          obs:  Observed field  [ny,nx]
    threshold:  The value threshold to verify (float)
 
   Optional inpurs:
      kernel:  The kernel radius of influence for forecasts
               and observations. (Ignored if not specified)
   """
   if kernel is not None:
     fcst = gen_buffer(fcst,kernel,threshold)
     obs = gen_buffer(obs,kernel,threshold)
   
   a,b,c,d = contingency(fcst,obs,threshold)
   FOH  = float(a)/float(a+b)
   POD  = float(a)/float(a+c)
   CSI  = float(a)/float(a+b+c)
   BIAS = float(a+b)/float(a+c)
   return FOH,POD,CSI,BIAS

def nmep(fcst,kernel,threshold,smooth=True):
   """
   Calculate the neighborhood maximum ensemble probability of an event
  
   Required inputs:
         fcst:  Ensemble forecast predicted field [nen,ny,nx]
       kernel:  The neighborhood search radiu
    threshold:  The value threshold to verify (float)

   Optional Inputs:
     smooth:  Apply smoothing via a guassian filter
              kernel radius is same as kernel
   """
   nen = fcst.shape[0]
   binary_ens = np.zeros(fcst.shape)

   for mindex in range(0,nen): #--- Get a binary array of 1 (exceed) or 0 (less than) threshold value
      binary_ens[mindex] = gen_buffer(fcst[mindex],kernel,threshold,ones=True)
   NMEP = np.sum(binary_ens,axis=0)/float(nen)

   if smooth: #--- Apply smoothing if desired
      NMEP =  gaussian_filter(NMEP,sigma = kernel)

   return NMEP


def createReliability(NEP,obs,threshold,perc = False,min_val=5.0):
   """
   Goal:  Make a reliability curve
           based upon selected swath of information
   
     Need: NEP, Obs, Threshold
     Note: NEP/Obs Shape= [ny,nx]
          Threshold  = single val
   
     Produced May 14, 2015
     Author: Jon Labriola/Nate Snook
   
    Modifications: None
   """

   if perc == True:
     tmp = np.where(obs<min_val,nan,obs)
     thresh_tmp = np.nanpercentile(tmp,threshold)
     threshold = thresh_tmp
     print('Threshold= '+str(threshold))

   [ny,nx] = NEP.shape
   xmin = 0
   ymin = 0
   xmax = nx
   ymax = ny

   bin_width = 0.05
   reliability_bins = np.arange(0.00, 1.00 + bin_width, bin_width)
   num_bins = len(reliability_bins) - 1

   #Reliability Diagram Variables
   obs_frequency = np.zeros((num_bins))
   bin_centers = np.zeros((num_bins))
   bin_count = np.zeros((num_bins))


   total_points = 0
   occur_points = 0
   for ix in np.arange(xmin, xmax, 1):
      for jy in np.arange(ymin, ymax, 1):
         total_points = total_points + 1
         if (obs[jy][ix] > threshold):
            occur_points = occur_points + 1

   obs_climo = float(occur_points) / float(total_points)

   for cur_bin in np.arange(0, num_bins, 1):
      samples = 0
      obs_yes = 0
      obs_no = 0

      for i in np.arange(xmin, xmax, 1):
         for j in np.arange(ymin, ymax, 1):
            if( (NEP[j,i] >= reliability_bins[cur_bin]) and (NEP[j,i] <= reliability_bins[cur_bin + 1])):
               samples = samples + 1   #We found a forecast within the current bin
               if (obs[j,i] >= threshold):
                  obs_yes = obs_yes + 1 #The criteria was met in the obs. (positive observation)
               elif (obs[j,i] < threshold):
                  obs_no = obs_no + 1   #The criteria was not met in the obs. (negative observation)

      #Calculate observed frequency and center of forecast probability bin
      bin_centers[cur_bin] = (reliability_bins[cur_bin] + reliability_bins[cur_bin + 1]) / 2.0

      if(samples > 0):
         obs_frequency[cur_bin] = float(obs_yes) / float(obs_yes + obs_no)
      else:
         obs_frequency[cur_bin] = -999 #No samples -- no data!  Avoid divide-by-zero.
      bin_count[cur_bin] = samples

   sample_climo = np.zeros((1001))
   no_skill = np.zeros((1001))
   bin_climo = np.zeros((num_bins))

   for i in np.arange(0, num_bins, 1):
      bin_climo[i] = obs_climo

   for index, i in enumerate(np.arange(0, 1.001, 0.001)):
      sample_climo[index] = obs_climo
      no_skill[index] = 0.5 * (obs_climo + i)

   return sample_climo,no_skill,obs_frequency,bin_centers,bin_climo,reliability_bins,bin_count


def pmmean(var2D_ens):
   """
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Goal:  Make the probability Matched Mean of
   !        an ensemble  quantity
   !
   ! Need: Ensemble Vals
   ! Note: Ens Shape= [n_ens,ny,nx,]
   !
   ! Produced June 5, 2015
   ! Author: Jon Labriola
   !
   ! Modifications: Flipped N_ens to be the first
   !                coordinate - Easier to put in
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   """

   n_ens = var2D_ens.shape[0]
   ny = var2D_ens.shape[1]
   nx = var2D_ens.shape[2]

   ensmean = np.mean(var2D_ens, axis=0) #Calculate ensemble mean

   all_values = []
   for jy in np.arange(0, ny):
      for ix in np.arange(0, nx):
         for member in np.arange(0, n_ens):
            all_values.append(var2D_ens[member, jy, ix])

   #Construct a tethered array to track where our ranked values are.
   tethered_array = np.zeros((ny, nx))
   for jy in np.arange(0, ny):
      for ix in np.arange(0, nx):
         tethered_array[jy, ix] = jy + (ny * ix)

   #Construct 1-D arrays to rank values
   ensmean_1d = []
   tethered_1d = []
   for jy in np.arange(0, ny):
      for ix in np.arange(0, nx):
         ensmean_1d.append(ensmean[jy, ix])
         tethered_1d.append(tethered_array[jy, ix])

   #Prepare the all_values array by sorting it from high to low:
   all_values.sort()
   all_values.reverse()

   #Tethering array elements to their position as tuples ensures the sort is reversible.
   position_tuples = zip(tethered_1d, ensmean_1d)
   #The lambda performs a "mini-function" -- google "python lambda notation" for info.
   sorted_tuples = sorted(position_tuples, key = lambda x: x[1])
   sorted_tuples.reverse()  #Put largest values first
   sorted_tuples = zip(*sorted_tuples)
   tethered_1d = sorted_tuples[0]
   ensmean_1d = sorted_tuples[1]

   pm_mean = np.zeros((ny, nx))

   for point in np.arange(0, nx*ny, 1):
      #The tethered array stores information on the position of the ranked points in ensmean.
      #We must first convert this data back into i- and j-locations, then assign the corresponding data.
      loc_jrank = np.mod(tethered_1d[point], ny)
      loc_irank = np.floor(int(tethered_1d[point] / ny))

      pm_mean[loc_jrank, loc_irank] = all_values[0 + (n_ens * point)]

   return pm_mean

def ROC(NEP,obs,threshold):
   """
   Calculate the Relative Operating Characteristic
   Requires: NEP,obs,threshold


   Returns: POD, POFD

   """
   probthresh = np.arange(0.00, 1.00, 0.01) #define your probability thresholds
   pod = np.zeros(len(probthresh) + 1)
   pofd = np.zeros(len(probthresh) + 1)  #False Alarm Rate

   for counter, curprob in enumerate(probthresh):  # Counter is 1,2,3,4... curprob = the probability threshold
      #fcst_yes = ma.masked_where(NEP < curprob, NEP) # Forecasted
      #fcst_no  = ma.masked_where(NEP >= curprob, NEP)

      fcst_yes = np.ma.masked_where(NEP < curprob, NEP) # Forecasted
      fcst_no  = np.ma.masked_where(NEP >= curprob, NEP)

      hit_array  = np.ma.masked_where(obs < threshold, fcst_yes)
      miss_array = np.ma.masked_where(obs < threshold, fcst_no)
      fa_array   = np.ma.masked_where(obs >= threshold, fcst_yes)
      cn_array   = np.ma.masked_where(obs >= threshold, fcst_no)

      #Calculate probability of detection (pod) and false alarm rate (far)
      #Total events = hits + misses
      #Total non-events = false alarms + correct nos

      hits = np.ma.count(hit_array)
      false_alarms = np.ma.count(fa_array)
      misses = np.ma.count(miss_array)
      correct_no = np.ma.count(cn_array)

      pod[counter] = float(hits) / float(hits + misses)    #Probability of detection: (hits) / (total events)
      pofd[counter] = float(false_alarms) / float(correct_no + false_alarms)

   return pod,pofd


def auc(NEP,obs,threshold):
   """
   Calculate the area underneath the ROC Curve

   """
   probthresh = np.arange(0.00, 1.00, 0.01) #define your probability thresholds
   pod = np.zeros(len(probthresh) + 1)
   pofd = np.zeros(len(probthresh) + 1)  #False Alarm Rate

   for counter, curprob in enumerate(probthresh):  # Counter is 1,2,3,4... curprob = the probability threshold
      #fcst_yes = ma.masked_where(NEP < curprob, NEP) # Forecasted
      #fcst_no  = ma.masked_where(NEP >= curprob, NEP)

      fcst_yes = np.ma.masked_where(NEP < curprob, NEP) # Forecasted
      fcst_no  = np.ma.masked_where(NEP >= curprob, NEP)

      hit_array  = np.ma.masked_where(obs < threshold, fcst_yes)
      miss_array = np.ma.masked_where(obs < threshold, fcst_no)
      fa_array   = np.ma.masked_where(obs >= threshold, fcst_yes)
      cn_array   = np.ma.masked_where(obs >= threshold, fcst_no)

      #Calculate probability of detection (pod) and false alarm rate (far)
      #Total events = hits + misses
      #Total non-events = false alarms + correct nos

      hits = np.ma.count(hit_array)
      false_alarms = np.ma.count(fa_array)
      misses = np.ma.count(miss_array)
      correct_no = np.ma.count(cn_array)

      pod[counter] = float(hits) / float(hits + misses)    #Probability of detection: (hits) / (total events)
      pofd[counter] = float(false_alarms) / float(correct_no + false_alarms)  #Probability of false detection: (false alarms / total non-events)


   #Now calculate AUC:
   AUC = -np.trapz(pod, pofd) #Area under curve using trapezoidal approximation.
   return AUC,pod,pofd

def BSS(NEP,obs,threshold):
   """
   Calculate the Briers Skill Score

   NEP = [ny x nx] Array
   obs = [ny x nx] Array
   threshold = Float

   return BSS (float)
   """


   [ny,nx] = obs.shape
   obs = np.where(obs>threshold,1,0)
   BS = np.mean((NEP - obs)**2)
   climo_freq = float(np.sum(obs))/float(ny*nx)
   BSref = np.mean((climo_freq - obs) **2)
   return 1- BS/float(BSref)
