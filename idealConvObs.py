import numpy as np
import argparse
import interp,fwdop
from netCDF4 import Dataset
import plotting
#============================================================================#
#  This program is designed to generate conventional observations from CM1 
#  output to prepare for data assimilation.  Has the capability to simulate...
#  1.) Radiosonde
#  2.) Surface
#  3.) Profiler Observations
#  
#  Once complete the observations and necesary meta data will be saved 
#  within a netCDF file for later analysis


#--- These Values Only once an Experiment and can therfore be changed manually to avoid mistakes
sim_snd = True   # Simulate soundings
sim_sfc = True   # Simulate surface obs
sim_pro = True   # Simulate profilers

if sim_snd:
   #--- Block 1: Sounding
   snd_locx = [30000,250000,75000]   #--- The Location of the Radar (either Longitude [WRF] or Distance [CM1]
   snd_locy = [30000,250000,75000]   #--- The Location of the Radar (either Latitude  [WRF] or Distance [CM1]
   snd_hgt  = [415,0.0,0.0]          #--- The height of the radar station of sea level (JDL Working on This)
   snd_obs  = ['u','v','p','T','qv']    #--- Evaluated Sounding Observations
   snd_error = [2.0,2.0,4.0,1.0,0.005]  #--- Observation Errors for Each Type (Can Grow More Complex With Time)

if sim_sfc:
   #--- Block 2: Surface
   sfc_locx = [30000,250000,75000]   #--- The Location of the Radar (either Longitude [WRF] or Distance [CM1]
   sfc_locy = [30000,250000,75000]   #--- The Location of the Radar (either Latitude  [WRF] or Distance [CM1]
   sfc_hgt  = [415,0.0,0.0]          #--- The height of the radar station of sea level (JDL Working on This)
   sfc_obs  = ['u','v','p','T','qv']    #--- Evaluated Sounding Observations
   sfc_error = [2.0,2.0,4.0,1.0,0.005]  #--- Observation Errors for Each Type (Can Grow More Complex With Time)

if sim_pro:
   #--- Block 3: Profiler
   snd_locx = [30000,250000,75000]   #--- The Location of the Radar (either Longitude [WRF] or Distance [CM1]
   snd_locy = [30000,250000,75000]   #--- The Location of the Radar (either Latitude  [WRF] or Distance [CM1]
   snd_hgt  = [415,0.0,0.0]          #--- The height of the radar station of sea level (JDL Working on This)
   snd_obs  = ['u','v','p','T','qv']    #--- Evaluated Sounding Observations
   snd_error = [2.0,2.0,4.0,1.0,0.005]  #--- Observation Errors for Each Type (Can Grow More Complex With Time)

#--- Argparse provides file name
parser = argparse.ArgumentParser()
parser.add_argument("obs_path", type=str,            help = 'The location of the model file used to extract obs')
arguments = parser.parse_args()

if sim snd:
  nobs = len(snd_locx)
  for ob in nobs:
   
