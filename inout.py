import interp,fwdop
from netCDF4 import Dataset
import numpy as np

class obs_plat(object):
   """
   This class is designed to handle osberavtional data in an organized manner

   The obs platform class generates two seperate dictionaries:
      1.) model - This dictionaary contains all model information including 
                  grid dimensions, model state variables, etc.
 
      2.) obs -   A nested dictionary that says all relvent observation information
                  via the addob function
 
                  Each observation platform has its unique nested dictionary, and then
                  each observation type has its own dictionary nested within the platform

                  This gives way to a dictionary structure that resembles
                  obs['PLATFORM_NAME']['OBS_NAME']['obx','oby','obz','elv','obs','obtype'...]

   Prior to saving informtion in the "obs" dictionary, the observation information is stored
   in the class itself for temporary modifications.
   """
   #--- Initialize Model Output
   def __init__(self,obtype,**kwargs):
      self.obtype = obtype
      self.obs = {}  #--- The longterm storage dictionaries for observations
      self.obx = None
      self.oby = None
      self.obz = None
      self.missing = -9999.
      #--- Adding Radar Specific Constants
      if obtype =='radar':
        #--- Add a variable to delineate clear_air
        if 'clear_air' in kwargs:
          self.clearair = kwargs['clear_air'] 
        else:  
          self.clearair = 0.
        self.az = None
        self.elv = None
       

   #--- Read Model Output (For Forward Operator)
   def readmod(self,model,path):
      """
      Read in model output that can be used to extract simulated observations
      """
      #--- Create a dictionary Array of Model Output
      if model.upper() == 'CM1':
         #varnames = {'output name in model':'desired save name'}
         varnames = {'ua':'u','va':'v','wa':'w','dbz':'dbz','rho':'rho',
                     'prs':'p','theta':'pt','qv':'qv','xh':'xh','yh':'yh','zh':'zh'} 
         self.model = {}
         #--- Read In 
         #xmin  =150 
         #xmax = 200
         #ymin = 150 
         #ymax = 200
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
            self.model[varnames[var]] = var_tmp
      
         dumpfile.close() 

         #--- Other Grid Parameters
         self.model['dx'] = self.model['xh'][1] - self.model['xh'][0] 
         self.model['nx'] = self.model['xh'].shape[0]
         self.model['dy'] = self.model['yh'][1]-self.model['xh'][0]
         self.model['ny'] = self.model['yh'].shape[0]
         self.model['nz'] = self.model['zh'].shape[0]

         #--- Making 2D Planes of x/y coordinates
         self.model['xh2d'] = np.tile(self.model['xh'], (self.model['ny'],1))
         self.model['yh2d'] = np.transpose(np.tile(self.model['yh'], (self.model['nx'],1)))

   def estab_platform(self,plat_name,x,y,z,tilts=None):
       """
       Define x,y,z observation locations for platform
       Defines platform name (for long-term storage in dictionary)
 
       Required inputs:
          plat_name:   The name that the platform will be saved under str
          x        :   The x-location of the platform  int
          y        :   The y-location of the platform  int
          z        :   The z-location of the platform  int

       Optional inputs:
          tilts    :   Array of defined radar tiles [ntilt]                       
       """
       self.platform = {}
       self.xloc = x
       self.yloc = y
       self.zloc = z
       self.plat_name = plat_name
   
       #--- Save Platform Information in Dictionary      
       self.obs[self.plat_name] = {}
       self.obs[self.plat_name]['plat_x']  = self.xloc
       self.obs[self.plat_name]['plat_y']  = self.yloc
       self.obs[self.plat_name]['plat_z']  = self.zloc  

       #--- Adding Radar Capabilities
       if self.obtype == 'radar':
          self.tilts = tilts
          self.ntilt = len(tilts)
          self.obs[self.plat_name]['tilts'] = self.tilts

   def obloc(self):
      """
      Calculating radar observation locations that will be saved to the object
      Must run estab_platform to define radar location first
      """
      if self.obtype == 'radar':
         self.obx,self.oby,self.obz,self.elv,self.az = interp.rad_obs_loc(self.model,self.xloc,self.yloc,self.zloc,self.tilts)     
      else:
         print('Other observations types are not accomodated yet')


   def conv_operator(self,varname,**kwargs):
      """
      Interpolate observations of the same type (profileof the same type (profiof the same type (profiler,sounding,surface))   
      Required arguments:
         varname:  A string of the variable name   

      Optional arguments:
         xloc = array of x-location(s)  [nobs]
         yloc = array of y-location(s)  [nobs]
         zloc = array of z-location(s)  [nobs]
   
      If any parameter is not defined than refer to station location
      """
      #--- Select Observation location
      if 'xloc' in kwargs:
          self.obx = kwargs['xloc']
      else:
          self.obx = [self.xloc]

      if 'yloc' in kwargs:
          self.oby = kwargs['yloc']
      else:
          self.oby = [self.xloc]

      if 'zloc' in kwargs:
          self.obz = kwargs['zloc']
      else:
          self.obz = [self.xloc]

      self.ob = np.zeros((len(self.obz)))

      #--- Loop Through observation locations
      for zindex,zloc in enumerate(self.obz):
         xloc = self.obx[zindex]
         yloc = self.oby[zindex]
   
         if varname.upper() in ['RADIOSONDE_U_WIND_COMPONENT','U_WIND_10M']:
            self.ob[zindex] = interp.point_interp(self.model,'u',xloc,yloc,zloc)

         elif varname.upper() in ['RADIOSONDE_V_WIND_COMPONENT','V_WIND_10M']:
            self.ob[zindex] = interp.point_interp(self.model,'v',xloc,yloc,zloc)    

         elif varname.upper() in ['RADIOSONDE_TEMPERATURE','TEMPERATURE_2M']:
            p = interp.point_interp(self.model,'p',xloc,yloc,zloc)
            pt = interp.point_interp(self.model,'pt',xloc,yloc,zloc)
            self.ob[zindex] = fwdop.theta_to_temp(pt,p)
            print('JDL Come Back to Make sure you dont have to keep Temperature in Kelvin')

         elif varname.upper() in ['RADIOSONDE_SURFACE_PRESSURE','SURFACE_PRESSURE']:
            self.ob[zindex] = interp.point_interp(self.model,'p',xloc,yloc,zloc)
            print('JDL Come Back to Make Sure the Assimilated units of pressure is correct')

         elif varname.upper() in ['RADIOSONDE_SPECIFIC_HUMIDITY','SPECIFIC_HUMIDITY_2M']: #--- JDL Does CM1 use qv or specific humidity?
            self.ob[zindex] = interp.point_interp(self.model,'qv',xloc,yloc,zloc)
            print('JDL Come Back to Make sure you dont have to convert to specific humidity') 

         else:
            print('Observation Unknown: %s'%varname)


   def radar_operator(self,dbz_only=False):
      """
      Observation operator given model state variables
    
      Optional Arguments:
         dbz_only: A string of the observation to generatee
      """
      if dbz_only:
         self.obdbz = fwdop.calcHx(self.model,self.obx,self.oby,self.obz,self.elv,self.az,
                                            clear_dbz = self.clearair,cartesian=True,dbzOnly=True)
      else:
         self.obvr,self.obdbz = fwdop.calcHx(self.model,self.obx,self.oby,self.obz,self.elv,self.az,
                                            clear_dbz = self.clearair,cartesian=True) 




   def superobs(self,nskip,roi):
      """
      A function used to interpolate observed data to a new grid

      The function is currently only compatible with radar data
      """
      #--- Save the fine observations
      fine_dbz = self.obdbz
      fine_vr = self.obvr
      fine_x = self.obx
      fine_y = self.oby

      #--- Coarsen existing grid
      self.obx = self.obx[:,::nskip,::nskip]
      self.oby = self.oby[:,::nskip,::nskip]
      self.obz = self.obz[:,::nskip,::nskip]
      self.az  = self.az[:,::nskip,::nskip]
      self.elv = self.elv[:,::nskip,::nskip]

      #--- Perform cressman interpolation for each radar tilt
      self.obdbz = interp.cressman(fine_x,fine_y,fine_dbz,self.obx,self.oby,roi)
      self.obvr  = interp.cressman(fine_x,fine_y,fine_vr,self.obx,self.oby,roi)


   def addob(self,obname,obtype,error,obdbz=False,obvr=False):
      """   
      The saved variables in a nested dictionary array include: 
         obname   :  The name the the observations will be stored under 
                     (same for one type of platform)
         obs      :  The observed values
         obtype   :  The DART code for the observation
         xloc     :  The x-location of the observation
         yloc     :  The y-location of the observation
         zloc     :  The z-location of the observation
         error    :  The observation error
         missing_flag   : A flag to tell DART what observations to ignore

         Radar Exclusive Variables:
         elevation     :  The radar tilt of the radar    (if needed)
         azimuth       :  The azimuth angle of the radar (if needed)
         clear_air     :  The threshold for clear air obs

      All Variables should have the same dimension except missing_flag

      The structure of the nested dictionary is as follows:
         dictionary['OBS_PLATFORM']['OBSNAME']['obsx','obsy','obsz'...]

      """
      #--- Saving the observation platform location
      self.obs[self.plat_name][obname]  = {}

      #--- Saving Captured Observations (special excpetions for dbz/vr)
      if obdbz :
         self.obs[self.plat_name][obname]['obs']     = self.obdbz
      elif obvr:
         self.obs[self.plat_name][obname]['obs']     = self.obvr
      else:
         self.obs[self.plat_name][obname]['obs']     = self.ob

      self.obs[self.plat_name][obname]['xloc']    = self.obx
      self.obs[self.plat_name][obname]['yloc']    = self.oby
      self.obs[self.plat_name][obname]['zloc']    = self.obz
      self.obs[self.plat_name][obname]['obstype'] = obtype
      self.obs[self.plat_name][obname]['error']   = error
      self.obs[self.plat_name][obname]['missing_flag'] = self.missing

      #--- Radar Important Information
      if obdbz or obvr:
         try:
            self.obs[self.plat_name][obname]['azimuth'] = self.az
         except:
            print('Unable to save azimuth information')
         try:
            self.obs[self.plat_name][obname]['elevation'] = self.elv
         except:
            print('Unable to save elevation information')
         try:
            self.obs[self.plat_name][obname]['clearair'] = self.clearair 
         except:
            print('Unable to save clearair information') 

def obcode():
   """
   Grab the observation codes that can be used in DART Data Assimilation

   Returns
     obcode - A dictionary that contains the observation code numbers
   """
   lines = open('/work/jonathan.labriola/python_scripts/pyDART/data/DART_obs.csv','r')
   obcode = {}
   for lindex, line in enumerate(lines):
      if lindex > 0:
         ln = line.split(', ')
         obcode[ln[0]] = int(ln[-1])
   return obcode

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
