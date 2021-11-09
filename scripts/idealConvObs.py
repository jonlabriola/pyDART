import numpy as np
import argparse
#import interp,fwdop
#from netCDF4 import Dataset
#import observations,readvar
#import plotting
import pydart
import pickle
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
sim_snd = False   # Simulate soundings
sim_sfc = False   # Simulate surface obs
sim_pro = True    # Simulate profilers
output_path = '../output/' #--- Save conventional observations here 
missing_value = -9999.0
if sim_snd:
   #--- Block 1: Sounding
   xsnd = [105000]#,250000,75000]   #--- The Location of the Radar (either Longitude [WRF] or Distance [CM1]
   ysnd = [105000]#,250000,75000]   #--- The Location of the Radar (either Latitude  [WRF] or Distance [CM1]
   zsnd  = [0.0] #,0.0,0.0]          #--- The height of the radar station of sea level (JDL Working on This)
   snd_obs  = ['RADIOSONDE_U_WIND_COMPONENT','RADIOSONDE_V_WIND_COMPONENT',
               'RADIOSONDE_TEMPERATURE','RADIOSONDE_SPECIFIC_HUMIDITY']    #--- Evaluated Sounding Observations
   snd_err = [2.0,2.0,4.0,0.005]  #--- Observation Errors for Each Type (Can Grow More Complex With Time)

   #--- Vertical Ob Spacing
   zmax_snd = 15000.  #--- Max observation height is 15000 m
   zincr_snd = 500.   #--- Observation collected every 500 m

if sim_sfc:
   #--- Block 2: Surface
   xsfc = [50000, 100000, 150000, 200000,
           50000, 100000, 150000, 200000,
           50000, 100000, 150000, 200000,
           50000, 100000, 150000, 200000]    #--- The Location of the Radar (either Longitude [WRF] or Distance [CM1]

   ysfc = [ 25000, 25000, 25000, 25000,
            75000, 75000, 75000, 75000,
           125000,125000,125000,125000,
           175000,175000,175000,175000]      #--- The Location of the Radar (either Latitude  [WRF] or Distance [CM1]

   zsfc  = [0.0,0.0,0.0,0.0,
            0.0,0.0,0.0,0.0,
            0.0,0.0,0.0,0.0,
            0.0,0.0,0.0,0.0]                 #--- The height of the radar station of sea level (JDL Working on This)

   sfc_obs  = ['TEMPERATURE_2M','U_WIND_10M','V_WIND_10M',
               'SURFACE_PRESSURE','SPECIFIC_HUMIDITY_2M']    #--- Evaluated Sounding Observations

   sfc_err = [2.0,2.0,4.0,1.0,0.005]  #--- Observation Errors for Each Type (Can Grow More Complex With Time)


if sim_pro:

   xpro = [50000, 100000, 150000, 200000,
           50000, 100000, 150000, 200000,
           50000, 100000, 150000, 200000,
           50000, 100000, 150000, 200000]    #--- The Location of the Radar (either Longitude [WRF] or Distance [CM1]

   ypro = [ 25000, 25000, 25000, 25000,
            75000, 75000, 75000, 75000,
           125000,125000,125000,125000,
           175000,175000,175000,175000]      #--- The Location of the Radar (either Latitude  [WRF] or Distance [CM1]

   zpro  = [0.0,0.0,0.0,0.0,
            0.0,0.0,0.0,0.0,
            0.0,0.0,0.0,0.0,
            0.0,0.0,0.0,0.0]                 #--- The height of the radar station of sea level (JDL Working on This)
   
   pro_obs  = ['RADIOSONDE_U_WIND_COMPONENT','RADIOSONDE_V_WIND_COMPONENT',
               'RADIOSONDE_TEMPERATURE','RADIOSONDE_SPECIFIC_HUMIDITY']

   pro_err = [2.0,2.0,4.0,0.0005]  #--- Observation Errors for Each Type (Can Grow More Complex With Time)

   zmax_pro  = 3000.    #--- Maximum profiler height
   zincr_pro = 250.     #--- Observation increment

#--- Argparse provides file name
parser = argparse.ArgumentParser()
parser.add_argument("obs_path", type=str,            help = 'The location of the model file used to extract obs')
arguments = parser.parse_args()

#--- Step 0: Read in Observation Codes for DA
observation_codes = pydart.readvar.obcode() 
#--- Step 1: Generate Conventional Observation Class
convob = pydart.observations.obs_plat('conv')
#--- Step 2: Read in model information
convob.readmod('cm1',arguments.obs_path)

min_hgt = np.amin(convob.model['zh'])  #--- min observation height

#--- Simulate sounding Observations
if sim_snd:
   nobs = len(xsnd)
   for ob in range(0,nobs):
      #---STEP 3:  Define Sounding Information
      platform_name = 'snd_%03d'%(ob+1)
      convob.estab_platform(platform_name,xsnd[ob],ysnd[ob],zsnd[ob])
      #--- Loop over the sounding profiles for each observation
      for oindex,obvar in enumerate(snd_obs):
         #--- Define observation locations
         hgts = np.arange(zsnd[ob],zmax_snd+zincr_snd,zincr_snd)
         hgts = hgts[np.where(hgts > min_hgt)] #--- Throw out all heights less than min model height 
         y_locations = np.ones(hgts.shape)*ysnd[ob]
         x_locations = np.ones(hgts.shape)*xsnd[ob]
         #--- Get Observations
         convob.conv_operator(obvar,xloc=x_locations,yloc=y_locations,zloc=hgts)        

         #--- Define error / obcode
         oberr =  np.ones(hgts.shape)*snd_err[oindex]
         obtype = observation_codes[obvar]
         convob.addob(obvar,oberr)

#--- Simulate Profiler Observations
if sim_pro:
   nobs = len(xpro)
   for ob in range(0,nobs):
      #---STEP 3:  Define Sounding Information
      platform_name = 'pro_%03d'%(ob+1)
      convob.estab_platform(platform_name,xpro[ob],ypro[ob],zpro[ob])
      #--- Loop over the sounding profiles for each observation
      for oindex,obvar in enumerate(pro_obs):

         #--- Define observation locations
         hgts = np.arange(convob.zloc,zmax_pro+zincr_pro,zincr_pro)
         hgts = hgts[np.where(hgts > min_hgt)] #--- Throw out all heights less than min model height
         y_locations = np.ones(hgts.shape)*ypro[ob]
         x_locations = np.ones(hgts.shape)*xpro[ob]

         #--- Get Observations for heights
         convob.conv_operator(obvar,xloc=x_locations,yloc=y_locations,zloc=hgts)

         #--- Define error / DART observation code
         oberr =  np.ones(hgts.shape)*pro_err[oindex]
         obtype = observation_codes[obvar]
         convob.addob(obvar,oberr,seed=(ob*100)+oindex)

#--- Simulate Surface Observations
if sim_sfc:
   nobs = len(xsfc)
   for ob in range(0,nobs):
      #---STEP 3:  Define Sounding Information
      platform_name = 'sfc_%03d'%(ob+1)
      convob.estab_platform(platform_name,xsfc[ob],ysfc[ob],zsfc[ob])
      #--- Loop over the sounding profiles for each observation
      for oindex,obvar in enumerate(sfc_obs):

         #--- Define Height above ground 
         if obvar in ['TEMPERATURE_2M','SPECIFIC_HUMIDITY_2M']:
            if (2.+convob.zloc) > min_hgt:
               hgt = [2.+convob.zloc]
            else:  
               hgt = [min_hgt]

         elif obvar in ['SURFACE_PRESSURE']:
            if convob.zloc > min_hgt:
               hgt = [convob.zloc]
            else:
               hgt = [min_hgt]
 
         elif obvar in ['U_WIND_10M','V_WIND_10M']:
            if (10.+convob.zloc) > min_hgt:
               hgt = [10.+convob.zloc]
            else:
               hgt = [min_hgt]
         y_locations = np.ones(hgts.shape)*ysfc[ob]
         x_locations = np.ones(hgts.shape)*xsfc[ob]

         #--- Get Observations
         convob.conv_operator(obvar,xloc=x_locations,yloc=y_locations,zloc=hgts)

         #--- Define error / obcode
         oberr =  np.ones(hgts.shape)*sfc_err[oindex]
         #obtype = observation_codes[obvar]
         convob.addob(obvar,oberr)  

#--- Short Snippet to Plot out Variables
#print('The date is ...',convob.obs[ob]['date'])
obs_plat = convob.obs.keys()
for platform in obs_plat:
   print('Platform = ',platform)
#   data = convob.obs[platform].keys()
#   for data_ind in data:
#     print('stored data includes...',data_ind) 

output_path = '%s/conv_obs_%04d%02d%02d%02d%02d%02d.pickle'%(output_path,convob.model['year'],convob.model['month'],convob.model['day'],convob.model['hour'],convob.model['minute'],convob.model['second'])   
pickle.dump(convob,open(output_path, "wb" ) )
