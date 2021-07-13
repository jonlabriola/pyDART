import interp
from netCDF4 import Dataset
import numpy as np
#--- A list of functions used to read in or output variables

def read_cm1(path,varname):
   """
   A simple function designed to extract existing CM1 model state variables
   and put all variables on the save grid

   inputs:
      path:  The filepath to the existing data file
      varnames:  A list of all the existing variable names to extract
   """

   output = {}
   dumpfile = Dataset(path,"r",fortmat="NETCDF4")
   for var in varname:
      if var in ['xh','yh','zh']: #--- 1-D variables 
         var_tmp = np.squeeze(dumpfile.variables[var][:])
      elif var in ['ua']: #--- Unstagger on the x-axis
         var_tmp = interp.unstagger_grid(np.squeeze(dumpfile.variables[var][0,:,:,:]),2)
      elif var in ['va']: #--- Unstaggger on the y-axis
         var_tmp = interp.unstagger_grid(np.squeeze(dumpfile.variables[var][0,:,:,:]),1)
      elif var in ['wa']: #--- Unstagger on the z - axis
         var_tmp = interp.unstagger_grid(np.squeeze(dumpfile.variables[var][0,:,:,:]),0)
      else: #--- No need to worry about unatggering
          var_tmp = np.squeeze(dumpfile.variables[var][0,:,:,:])
      output[var] = var_tmp
   output['dy'] = output['yh'][1]- output['yh'][0]
   output['ny'] = output['yh'].shape[0]
   output['dx'] = output['xh'][1]- output['xh'][0]
   output['nx'] = output['xh'].shape[0]
   output['nz'] = output['zh'].shape[0]

   return output
