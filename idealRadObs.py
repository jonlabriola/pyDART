import numpy as np
import argparse
import interp,fwdop,inout
from netCDF4 import Dataset
import plotting
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
   output_path = 'radar_obs.pickle' #--- Output Radar Observation Object into a Pickle File
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

#--- Create Observation Class for Radars
radobs = inout.obs_plat('radar',clear_air = clearair_dbz)

#--- Open Model State Variables
radobs.readmod('cm1',arguments.obs_path)

for rdr in range(0,nrdr): #--- Loop over State Variables

   #---Step 1: Define the Observation Platform
   platform_name = 'radar%02d'%(rdr+1)
   radobs.estab_platform(platform_name,xloc[rdr],yloc[rdr],hgt[rdr],tilts)

   #--- Step 2: Calculate Radar Observation Locations in volume
   radobs.obloc()  

   #--- Step 3: Call Forward Operator
   radobs.forward_operator('zvr')
   if save_fine_obs: #--- Save fine Observations
      #--- Save the Fine Observations (dbz)
      dbztype = np.zeros(radobs.obdbz.shape).fill(13)         #--- REFLECTIVITY CODE for DART
      dbzerr = np.zeros(radobs.obdbz.shape).fill(dbz_error)   #--- REFLECTIVITY ERROR ASSUMPTIONS for DA
      radobs.addob("fine_z", dbztype, dbzerr, obdbz=True) 

      #--- Save the Fine Observations (vr)
      vrtype = np.zeros(radobs.obvr.shape).fill(11)       #--- VELOCITY CODE for DART
      vrerr = np.zeros(radobs.obvr.shape).fill(vr_error)  #--- VELOCITY ERROR ASSUMPTIONS FOR DA
      radobs.addob("fine_vr", vrtype, vrerr, obvr=True)

   #--- Step 4: Filter Observations to Coarse Grid
   if superobs:
#      #--- First Save the High-Resolution Observations
       radobs.superobs(nskip,roi=roi)
   
       dbztype = np.zeros(radobs.obdbz.shape).fill(13)         #--- REFLECTIVITY CODE for DART
       dbzerr = np.zeros(radobs.obdbz.shape).fill(dbz_error)   #--- REFLECTIVITY ERROR ASSUMPTIONS for DA
       radobs.addob("superob_z", dbztype, dbzerr, obdbz=True)

       vrtype = np.zeros(radobs.obvr.shape).fill(11)       #--- VELOCITY CODE for DART
       vrerr = np.zeros(radobs.obvr.shape).fill(vr_error)  #--- VELOCITY ERROR ASSUMPTIONS FOR DA
       radobs.addob("superob_vr", vrtype, vrerr, obvr=True) 


pickle.dump(radobs,open(output_path, "wb" ) )
#test_obj = pickle.load(open("radarobs.p", "rb" ))


#---Step 6 Create NetCDF file
#fn = 'radar_obs.nc'
#wrtfile = Dataset(fn, 'w', format='NETCDF4')
#if 'nx2' in mem:
#  wrtfile.setncattr('nx',mem['nx2'])
#  wrtfile.setncattr('ny',mem['ny2'])
#else:
#  wrtfile.setncattr('nx',mem['nx'])
#  wrtfile.setncattr('ny',mem['ny'])
#wrtfile.setncattr('ntilt',ntilt)
#wrtfile.setncattr('nrdr',nrdr)
#if 'nx2' in mem:
#   wrtfile.createDimension('yh', mem['ny2'])
#   wrtfile.createDimension('xh', mem['nx2'])
#else:
#   wrtfile.createDimension('yh', mem['ny'])
#   wrtfile.createDimension('xh', mem['nx'])
#wrtfile.createDimension('tilts', ntilt)
#wrtfile.createDimension('radars', nrdr)
#
#for key in rdrobs.keys():
#   wrtfile.createVariable(key, 'f4', ('radars','tilts','yh','xh'))
#   wrtfile.variables[key][:,:,:,:] = rdrobs[key][:,:,:,:]
#   if key in ['dbz','dbzerr']:
#      wrtfile.variables[key].units = 'dBZ'
#   elif key in ['vr','vrerr']:
#      wrtfile.variables[key].units = 'm s^{-1}'
#   elif key in ['rdrx','rdry','rdrz','x','y','z']:
#      wrtfile.variables[key].units = 'm'
#   elif key in ['az','elv']:
#      wrtfile.variables[key].units = 'radians'
##
#wrtfile.close()
