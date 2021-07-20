import numpy as np
import argparse
#import observations, readvar
#import interp,fwdop,observations
#from netCDF4 import Dataset
#import plotting
import pydart
import pickle
#--- These Values Only once an Experiment and can therfore be changed manually to avoid mistakes

#---- Static Parameters Throughout DA ----#

case_type = 'idealized' #-- Either Idealized (i.e., create your own observations) or Real (i.e., read observation from file)
if case_type == 'idealized':     #--- Must Explicitly List All Parameters
   xloc = [-100000,250000,75000]   #--- The Location of the Radar (either Longitude [WRF] or Distance [CM1]
   yloc = [-75000, 250000,0]   #--- The Location of the Radar (either Latitude  [WRF] or Distance [CM1]
   hgt  = [415,0.0,0.0]          #--- The height of the radar station of sea level (JDL Working on This)
   tilts = np.radians(np.array([0.5, 0.9, 1.3, 1.8, 2.4, 3.1, 
           4.0, 5.1, 6.4, 8.0, 10.0, 12.5, 15.6, 19.5])) #--- NEXRAD Tilts
   case_type = 'idealized'       #-- Real Date or Idealized
   clearair_dbz = 10.            #--- Threshold for clear air observations
   #missing = -999.               #--- Value of Missing Observation
   vr_error = 4.0                #--- Vr Error (m/s)
   dbz_error = 6.0               #--- dBZ Error (dBZ)
   save_fine_obs = True          #--- Flag to save fine observations (in addition to coarse obs)
   nrdr = len(xloc)
   output_path = './output/' #--- Output Radar Observation Object into a Pickle File
else:  #--- A Real-Data Case
   print('Work on this at a later time')


#--- Superobbing Information
superobs = True        #--- If True Apply Cressman Filter
nskip = 4              #--- The number of grid points to skip
roi = 6000.            #--- Cressman Function influence radius

#--- Thinning Observations (Need To Update if Desired)
#clear_thin = 6    #--- The number of grid points to skip where Z < clearir_dbz
#storm_thin = 3    #--- The number of grid points to skip where Z > clearair_dbz

#--- Argparse provides file name
parser = argparse.ArgumentParser()
parser.add_argument("obs_path", type=str,            help = 'The location of the radar observation file')
arguments = parser.parse_args()

#=================================#
#       Start Working Code        #
#=================================#
#--- Get Observation Codes
observation_codes = pydart.readvar.obcode()

#--- Create Observation Class for Radars
radobs = pydart.observations.obs_plat('radar',clear_air = clearair_dbz)

#--- Open Model State Variables
radobs.readmod('cm1',arguments.obs_path)

for rdr in range(0,nrdr): #--- Loop over State Variables

   #---Step 1: Define the Observation Platform
   platform_name = 'radar_%03d'%(rdr+1)
   radobs.estab_platform(platform_name,xloc[rdr],yloc[rdr],hgt[rdr],tilts)

   #--- Step 2: Calculate Radar Observation Locations in volume
   radobs.obloc()  

   #--- Step 3: Call Forward Operator
   radobs.radar_operator()

   if superobs : #---Perform Super Obbing if Desired 
      if save_fine_obs: #--- Save High-Res Obs if Superobbing
         #--- Save the Fine Observations (dbz)
         dbztype = np.zeros(radobs.obdbz.shape).fill(13)         #--- REFLECTIVITY CODE for DART
         dbzerr = np.zeros(radobs.obdbz.shape).fill(dbz_error)   #--- REFLECTIVITY ERROR ASSUMPTIONS for DA
         radobs.addob("fine_z", dbztype, dbzerr, obdbz=True) 

         #--- Save the Fine Observations (vr)
         vrtype = np.zeros(radobs.obvr.shape).fill(11)       #--- VELOCITY CODE for DART
         vrerr = np.zeros(radobs.obvr.shape).fill(vr_error)  #--- VELOCITY ERROR ASSUMPTIONS FOR DA
         radobs.addob("fine_vr", vrtype, vrerr, obvr=True)

      #--- Step 4: Filter Observations to Coarse Grid
      radobs.superobs(nskip,roi=roi)


   #--- JDL Eventually Split Up - REFLECTIVITY/CLEAR AIR REFLECITIVITY
   varname = 'RADAR_REFLECTIVITY'
   code = observation_codes[varname]
   dbztype = np.ones(radobs.obdbz.shape)*code         #--- REFLECTIVITY CODE for DART
   dbzerr  = np.ones(radobs.obdbz.shape)*dbz_error   #--- REFLECTIVITY ERROR ASSUMPTIONS for DA
   radobs.addob(varname, dbztype, dbzerr, obdbz=True)

   varname = 'DOPPLER_RADIAL_VELOCITY'
   code = observation_codes[varname]    
   vrtype = np.zeros(radobs.obvr.shape)*code       #--- VELOCITY CODE for DART
   vrerr = np.zeros(radobs.obvr.shape)*vr_error  #--- VELOCITY ERROR ASSUMPTIONS FOR DA
   radobs.addob("DOPPLER_RADIAL_VELOCITY", vrtype, vrerr, obvr=True) 

output_path = '%s/radar_obs_%04d%02d%02d%02d%02d.pickle'%(output_path,radobs.date['year'],radobs.date['month'],radobs.date['day'],radobs.date['hour'],radobs.date['minute'])
pickle.dump(radobs,open(output_path, "wb" ) )
