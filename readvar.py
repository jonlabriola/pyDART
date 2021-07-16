import datetime as dt
from netCDF4 import Dataset
import interp
import numpy as np

def read_cm1(path,time_str=None):
   """
   Read in model output that can be used to extract simulated observations
   """
   #--- Create a dictionary Array of Model Output
   #varnames = {'output name in model':'desired save name'}
   varnames = {'ua':'u','va':'v','wa':'w','dbz':'dbz','rho':'rho',
               'prs':'p','theta':'pt','qv':'qv','xh':'xh','yh':'yh','zh':'zh'}
   model = {}
   #--- Read In
   #xmin  =150 #xmax = 200 #ymin = 150 #ymax = 200
   xmin = 0
   xmax = 300
   ymin = 0
   ymax = 300
   dumpfile = Dataset(path,"r",fortmat='NETCDF4')
   for var in varnames.keys():
      if var in ['xh','yh','zh']: #--- 1-D variables
         var_tmp = np.squeeze(dumpfile.variables[var][:])
         if var in ['xh']: var_tmp = var_tmp[xmin:xmax]
         elif var in ['yh']: var_tmp = var_tmp[ymin:ymax]
      elif var in ['ua','va','wa']: #--- Unstagger the grid for certain vars
         if var in ['ua']: unstag_ax = 2 #--- x-axis
         if var in ['va']: unstag_ax = 1 #--- y-axis
         if var in ['wa']: unstag_ax = 0 #--- z-axis
         var_tmp = interp.unstagger_grid(np.squeeze(dumpfile.variables[var][0,:,:,:]),unstag_ax)
         var_tmp = var_tmp[:,ymin:ymax,xmin:xmax]
      else: #--- No Staggered Grids
         var_tmp = np.squeeze(dumpfile.variables[var][0,:,:,:])
         var_tmp = var_tmp[:,ymin:ymax,xmin:xmax]
      model[varnames[var]] = var_tmp

   #--- Other Grid Parameters
   model['dx'] = model['xh'][1] - model['xh'][0]
   model['nx'] = model['xh'].shape[0]
   model['dy'] = model['yh'][1]-model['xh'][0]
   model['ny'] = model['yh'].shape[0]
   model['nz'] = model['zh'].shape[0]

   #--- Making 2D Planes of x/y coordinates
   model['xh2d'] = np.tile(model['xh'], (model['ny'],1))
   model['yh2d'] = np.transpose(np.tile(model['yh'], (model['nx'],1)))


   if time_str is not None:
      #--- Add Time (CM1 Only Includes Initial Time So Must Add Model Time)
      for time in time_str:
         model[time] = dumpfile.getncattr(time)
      epoch_time = int(dumpfile.variables['time'][0])
      init_time = dt.datetime(model['year'],model['month'],model['day'],
                                model['hour'],model['minute'],model['second'])
      mod_time = init_time + dt.timedelta(seconds=epoch_time)
      mod_time.timetuple()

      for tindex,time in enumerate(time_str):
         model[time] = mod_time.timetuple()[tindex]


   dumpfile.close()
   return model

def obcode():
   """
   Grab the observation codes that can be used in DART Data Assimilation

   Returns
     obcode - A dictionary that contains the observation code numbers
   """
   lines = open('./data/DART_obs.csv','r')
   obcode = {}
   for lindex, line in enumerate(lines):
      if lindex > 0:
         ln = line.split(', ')
         obcode[ln[0]] = int(ln[-1])
   return obcode
