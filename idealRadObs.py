import numpy as np
import argparse
import os
import location,fwdop
from netCDF4 import Dataset
import matplotlib
import matplotlib.pyplot as plt
import plotting

#--- These Values Only once an Experiment and can therfore be changed manually to avoid mistakes

#---- Static Parameters Throughout DA ----#

case_type = 'idealized' #-- Either Idealized (i.e., create your own observations) or Real (i.e., read observation from file)
if case_type == 'idealized': #--- Must Explicitly List All Parameters
   locx = [30000,250000,75000]  #--- The Location of the Radar (either Longitude [WRF] or Distance [CM1]
   locy = [30000,250000,75000]  #--- The Location of the Radar (either Latitude  [WRF] or Distance [CM1]
   hgt  = [415,0.0,0.0]         #--- The height of the radar station of sea level (JDL Working on This)
   tilts = np.radians(np.array([0.5, 0.9, 1.3, 1.8, 2.4, 3.1, 4.0, 5.1, 6.4, 8.0, 10.0, 12.5, 15.6, 19.5])) #--- NEXRAD Tilts
   #tilts = np.radians(np.array([0.5]))
   case_type = 'idealized' #-- Real Date or Idealized
   var=['vr','refl'] # The Observations to extract
   clearair_dbz = 10.
   missing = -999.
   vr_error = 4.0
   dbz_error = 6.0
   #clear_thin = 6 #--- The number of grid points to skip 
   #storm_thin = 3
   smoothing = True #--- Perform cressman filter on output'
else:  #--- A Real-Data Case
   print('Work on this at a later time')


#--- Superobbing Information
superobs = True #---> If True Apply Cressman Filter
new_grid = 5000. #---> New Horizontal Grid Spacing
influence_rad = 6000. #--- Cressman Function influence radius

#--- Argparse provides file name
parser = argparse.ArgumentParser()
parser.add_argument("obs_path", type=str,            help = 'The location of the radar observation file')
arguments = parser.parse_args()


#--- Step 1: Grab Model Data (Currently CM1)
mem = {}
dumpfile = Dataset(arguments.obs_path,"r",fortmat="NETCDF4")
varname = ['ua','va','wa','dbz','rho','xh','yh','zh']
for var in varname:
   if var in ['xh','yh','zh']:
      var_tmp = np.squeeze(dumpfile.variables[var][:])
   elif var in ['ua']:
      var_tmp = location.unstagger_grid(np.squeeze(dumpfile.variables[var][0,:,:,:]),2)
   elif var in ['va']:
      var_tmp = location.unstagger_grid(np.squeeze(dumpfile.variables[var][0,:,:,:]),1)
   elif var in ['wa']:
      var_tmp = location.unstagger_grid(np.squeeze(dumpfile.variables[var][0,:,:,:]),0)
   else:
       var_tmp = np.squeeze(dumpfile.variables[var][0,:,:,:])
   mem[var] = var_tmp
mem['dy'] = mem['yh'][1]- mem['yh'][0]
mem['ny'] = mem['yh'].shape[0]
mem['dx'] = mem['xh'][1]- mem['xh'][0]
mem['nx'] = mem['xh'].shape[0]
mem['nz'] = mem['zh'].shape[0]


#--- Observation Location Info
nrdr = len(locx)
ntilt = tilts.shape[0]
obloc = {}  #--- radar observation locations

#--- Determin Grid of Coarse Obs
if superobs:
   #--- Coarsened data is saved with a 2 symbol
   nskip_x = int(np.floor(new_grid/mem['dx']))
   nskip_y = int(np.floor(new_grid/mem['dy']))
   mem['xh2'] = mem['xh'][::nskip_x]
   mem['yh2'] = mem['yh'][::nskip_y]  
   mem['nx2'] = mem['xh2'].shape[0]
   mem['ny2'] = mem['yh2'].shape[0]
   mem['dx2'] = mem['xh2'][1]-mem['xh2'][0]
   mem['dy2'] = mem['yh2'][1]-mem['yh2'][0]


#--- JDL If your decide to further thin observations make sure 
#--- to update this block


#--- Define your final entry variables to save
rdrobs = {}
for varname in ['dbz','vr','x','y','z','elv','az','rdrx','rdry','rdrz','dbzerr','vrerr']:
   if 'xh2' in mem:
      rdrobs[varname] = np.zeros((nrdr,ntilt,mem['ny2'],mem['nx2']))
   else: 
      rdrobs[varname]  = np.zeros((nrdr,ntilt,mem['ny'],mem['nx'])) 


for rdr in range(0,nrdr):
   for rtilt,radtilt in enumerate(tilts):
      tmpob = {}
      #--- Step 2: Calculate Radar Observation Locations for Tilts
      tmpob['x'],tmpob['y'],tmpob['z'],tmpob['elv'],tmpob['az'] = location.rad_obs_loc(mem,locx[rdr],locy[rdr],hgt[rdr],radtilt)

      #--- Step 3: Call Forward Operator
      obtype = np.zeros((np.size(tmpob['x'])))
      obtype[:] = 11 
      tmpob['vr'],tmpob['dbz'] = fwdop.calcHx(mem, obtype, tmpob['x'], tmpob['y'], tmpob['z'], tmpob['elv'], tmpob['az'], 
                                       clear_dbz= clearair_dbz,cartesian=True)


      #dbz_threed = np.reshape(tmpob['dbz'],(mem['ny'],mem['nx']),order='F')
      #plotting.rough_plot(dbz_threed,'dbz')

      #--- Plot observations if desired
      #vr_threed = np.reshape(obs_vr,(mem['ny'],mem['nx']),order='F') 
      #plotting.rough_plot(vr_threed,'vr')

      #--- Step 4: Filter Observations to Coarse Grid
      if superobs:
         #--- First Save the High-Resolution Observations
         fine_x = tmpob['x']
         fine_y = tmpob['y']
         fine_dbz = tmpob['dbz']
         fine_vr = tmpob['vr']

         #--- Coarsen Observation Grid by nskip (Need to Spatially Regrid To Do This) 
         for key in tmpob.keys():
            tmpob[key] = np.reshape(tmpob[key],(mem['ny'],mem['nx']),order='F')[::nskip_y,::nskip_x].flatten(order='F')
         
         #dbz_threed = np.reshape(tmpob['dbz'],(mem['ny2'],mem['nx2']),order='F')
         #plotting.rough_plot(dbz_threed,"dbz")

         #--- Cressman Filter #--- JDL Need to check the cressman function is not complete
         print('Applying Cressman Filter to dBZ')
         tmpob['dbz'] = location.cressman(fine_x,fine_y,fine_dbz,tmpob['x'],tmpob['y'],influence_rad)

         print('Applying Cressman Filter to Vr')
         tmpob['vr'] = location.cressman(fine_x,fine_y,fine_vr,tmpob['x'],tmpob['y'],influence_rad)
     

      #--- Conduct Test Plot
      #dbz_threed = np.reshape(tmpob['dbz'],(mem['ny2'],mem['nx2']),order='F')
      #plotting.rough_plot(dbz_threed,"dbz") 

      #--- Step 5 Thin Radar Data (If Desired)

      
      #--- Step 6 Return to a 4-Dimensional Array (nrdr,tilt,ny,nx) to prepare for storage
      if 'nx2' in mem:
         ny_save = mem['ny2']
         nx_save = mem['nx2']
      else:
         ny_save = mem['ny']
         nx_save = mem['nx']

      for key in rdrobs.keys():
         if key == 'rdrx':    rdrobs[key][rdr,rtilt] = locx[rdr] #--- Radar X-Location
         elif key == 'rdry':  rdrobs[key][rdr,rtilt] = locy[rdr] #--- Radar Y-Location
         elif key == 'rdrz':  rdrobs[key][rdr,rtilt] = hgt[rdr]  #--- Radar Z-Location
         elif key == 'dbzerr': rdrobs[key][rdr,rtilt] = dbz_error #--- Observation Error
         elif key == 'vrerr':  rdrobs[key][rdr,rtilt] = vr_error
         else:
            rdrobs[key][rdr,rtilt] = np.reshape(tmpob[key],(ny_save,nx_save),order='F')


#---Step 6 Create NetCDF file
fn = 'radar_obs.nc'
wrtfile = Dataset(fn, 'w', format='NETCDF4')
if 'nx2' in mem:
  wrtfile.setncattr('nx',mem['nx2'])
  wrtfile.setncattr('ny',mem['ny2'])
else:
  wrtfile.setncattr('nx',mem['nx'])
  wrtfile.setncattr('ny',mem['ny'])
wrtfile.setncattr('ntilt',ntilt)
wrtfile.setncattr('nrdr',nrdr)
if 'nx2' in mem:
   wrtfile.createDimension('yh', mem['ny2'])
   wrtfile.createDimension('xh', mem['nx2'])
else:
   wrtfile.createDimension('yh', mem['ny'])
   wrtfile.createDimension('xh', mem['nx'])
wrtfile.createDimension('tilts', ntilt)
wrtfile.createDimension('radars', nrdr)

#wrtfile.variables['tilts'], 'f4', ('tilts'))
for key in rdrobs.keys():
   wrtfile.createVariable(key, 'f4', ('radars','tilts','yh','xh'))
   wrtfile.variables[key][:,:,:,:] = rdrobs[key][:,:,:,:]
   
   if key in ['dbz','dbzerr']:
      wrtfile.variables[key].units = 'dBZ'
   elif key in ['vr','vrerr']:
      wrtfile.variables[key].units = 'm s^{-1}'
   elif key in ['rdrx','rdry','rdrz','x','y','z']:
      wrtfile.variables[key].units = 'm'
   elif key in ['az','elv']:
      wrtfile.variables[key].units = 'radians'

wrtfile.close()
