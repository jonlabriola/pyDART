import datetime as dt
from netCDF4 import Dataset
import pydart
import numpy as np
import pkg_resources
#--- This module is designed to read in external information including model output,
#--- and text files

def read_cm1(path,time_str=None):
   """
   Read in model output that can be used to extract simulated observations
   You will try to read in one of two file types, Restart and Forecast

   Output:
      model - A dictionary that contains the saved model output
   """

   try:
      #--- Restart File
      varnames = {'ua':'u','va':'v','wa':'w','dbz':'dbz','rho':'rho',
                  'prs':'p','theta':'pt','qv':'qv','xh':'xh','yh':'yh','zh':'zh'}
      model = gen_model_output(path,varnames,time_str,rstfile=True)
      print('Working with a Restart File')
   except:
      #--- Forecast File 
      varnames = {'u':'u','v':'v','w':'w','dbz':'dbz','rho':'rho',
                  'prs':'p','th':'pt','qv':'qv','xh':'xh','yh':'yh','zh':'zh'}
      model = gen_model_output(path,varnames,time_str,rstfile=False)
      print('Working with a Forecast File')

   return model


def gen_model_output(path,varnames,time_str,rstfile):
   """
   This function is used to read in model output 
   with the provided variable names.

   The variables are stored in the model dictionary which is passed down
   to make observations.

   Required input:
      path     - A string point to the model output file to be opened
      varname  - A dictionary containing all the variables you want to read
                 and what they will be refered to in the model dictionary
      time_str - A string that is used to denote the evaluated time.  If time_str is None
                 than the variable will be ignored 
      rstfile  - A Boolean used to determine if you are working with a restart file
                 if rstfile == false, than you are working with a forecast file

   Output:
      model - A dictionary that contains the saved model output.
   """
   model = {}
   dumpfile = Dataset(path,"r",fortmat='NETCDF4')
   if rstfile:
     scale_factor =  1
   else:
     scale_factor = 1000. #--- Forecast file are written in units of km, want to convert to m

   for var in varnames.keys():
      if var in ['xh','yh','zh']: #--- 1-D variables
         var_tmp = np.squeeze(dumpfile.variables[var][:])#*scale_factor
         #--- JDL Check to Make Meters
         if np.amax(var_tmp) < 5000:
            var_tmp = var_tmp * scale_factor
         if var in ['xh']: var_tmp = var_tmp[:]
         elif var in ['yh']: var_tmp = var_tmp[:]
         #--- Stagger grids in x/y directions, while deleting edges  (JDL NEW)
         #--- This is done so forecast and obs grids align in horizontal
         if var in ['xh','yh']:
            var_tmp = (var_tmp[0:-1] + var_tmp[1:])/2.
         if var in ['zh']:
            if var_tmp.ndim > 1: var_tmp = var_tmp[:,0,0] # Remove Extra Dimensions added by DART
      else:
         #--- We want all grids to be staggered in the horizontal
         #--- This makes for easy forecast verification when switching between grids
         #--- For w - unstagger the grid in the vertical
         var_tmp = np.squeeze(dumpfile.variables[var][0,:,:,:])
         if var in ['ua','u']:
            #--- Remove border gridpoints on edges of staggered x-grid and stagger y-axis
            var_tmp = pydart.interp.shift_grid(np.squeeze(dumpfile.variables[var][0,:,:,1:-1]),1) 
         elif var in ['va','v']:
            #--- Remove border gridpoints on edges of staggered y-grid and stagger x-axis
            var_tmp = pydart.interp.shift_grid(np.squeeze(dumpfile.variables[var][0,:,1:-1,:]),2) 
         else: #--- Grids initially unstaggered in Horizontal
            var_tmp = pydart.interp.shift_grid(var_tmp[:,:,:],2) #--- stagger x-axis
            var_tmp = pydart.interp.shift_grid(var_tmp[:,:,:],1) #--- stagger y-axis
            if var in ['wa','w']:
               var_tmp = pydart.interp.shift_grid(var_tmp[:,:,:],0) #--- Unstagger Z-axis
      model[varnames[var]] = var_tmp

   #--- Other Grid Parameters
   model['dx'] = (model['xh'][1] - model['xh'][0])
   model['nx'] = model['xh'].shape[0]
   model['dy'] = (model['yh'][1]-model['yh'][0]) 
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
   return(model)


def obcode():
   """
   Grab the observation codes that can be used in DART Data Assimilation

   Returns
     obcode - A dictionary that contains the observation code numbers
   """
   #lines = open('./data/DART_obs.csv','r')
   filepath = pkg_resources.resource_filename('pydart','data/DART_obs.csv')
   lines = open(filepath)
   obcode = {}
   for lindex, line in enumerate(lines):
      if lindex > 0:
         ln = line.split(', ')
         obcode[ln[0]] = int(ln[-1])
   return obcode
