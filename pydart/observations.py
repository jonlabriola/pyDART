#import interp,fwdop,readvar
import pydart
import datetime as dt
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
      self.time_str = ['year','month','day','hour','minute','second'] #--- order time is read
      self.obx = None
      self.oby = None
      self.obz = None
      self.max_hgt = None
      self.missing = 99999999.

      #--- Adding Radar Specific Constants
      if obtype =='radar':
        #--- Add a variable to delineate clear_air
        if 'clear_air' in kwargs:
          self.clearair = kwargs['clear_air'] 
        else:  
          self.clearair = 0.
        self.az = None
        self.elv = None
        self.nyquist = None 
   #--- Read Model Output (For Forward Operator)
   def readmod(self,model,path):
      """
      Read in model output that can be used to extract simulated observations
      """
      #--- Create a dictionary Array of Model Output
      if model.upper() == 'CM1':
         self.model = pydart.readvar.read_cm1(path,self.time_str)
      else:
         print('Add Capabilities to read other models')

   def estab_platform(self,plat_name,x,y,z,tilts=None,nyquist=None,date=None,max_hgt=None):
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
          nyquist  :   The nyquist velocity of the radar
          date     :   An array with defined time information
                       [year,month,day,hour,minute,second]
          max_hgt  :   The maximum height of an observation
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
          self.obs[self.plat_name]['nyquist']  = nyquist 
       #--- Setting observation time
       self.obs[self.plat_name]['date'] = {}
       for tindex,time in enumerate(self.time_str):
         if date is not None:          #--- Manually Input Date
            self.obs[self.plat_name]['date'][time] = date[tindex]
         elif hasattr(self,'model'):    #--- Grab time from model output
            self.obs[self.plat_name]['date'][time] = self.model[time]
         else:                         #--- Establish with None
            self.obs[self.plat_name]['date'][time] = None

       if max_hgt is None:
          self.max_hgt = 1E50  #--- Some very high observation (arbitrary)
       else:
          self.max_hgt = max_hgt

   def obloc(self):
      """
      Calculating radar observation locations that will be saved to the object
      Must run estab_platform to define radar location first
      """
      if self.obtype == 'radar':
         self.obx,self.oby,self.obz,self.elv,self.az = pydart.interp.rad_obs_loc(self.model,self.xloc,self.yloc,self.zloc,self.tilts,rad_top = self.max_hgt)    
      else:
         print('Other observations types are not accomodated yet')


   def conv_operator(self,varname,cloud_base_limit=False,refl_limit=False,**kwargs):
      """
      Interpolate observations of the same type (profileof the same type (profiof the same type (profiler,sounding,surface))   
      Required arguments:
         varname:  A string of the variable name   

      Optional arguments:
         cloud_base_limit: Boolean to consider obs above cloud base.  Currently Only limits Temperature and Qv fields (AERI)
         refl_limit: Boolean to consider obs where reflectivity is minimal (automatically called in cloud_base_limit called (e.g., DWL)
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
          self.oby = [self.yloc]

      if 'zloc' in kwargs:
          self.obz = kwargs['zloc']
      else:
          self.obz = [self.xloc]

      self.ob = np.zeros((len(self.obz)))

      
      #--- Loop Through observation locations
      for zindex,zloc in enumerate(self.obz):
         xloc = self.obx[zindex]
         yloc = self.oby[zindex]

         #----------------------------------#
         # This is the thresholding section #
         #----------------------------------#

         #--- Keep Obs Beneath Cloud Base (If Desired)
         if cloud_base_limit:
            cloud_conc = np.zeros(self.model['zh'].shape)
            for hindex, hgt in enumerate(self.model['zh']):
               if hindex > 0: cloud_conc[hindex] = pydart.interp.point_interp(self.model,'qc',xloc,yloc,hgt)
            if np.nanmax(cloud_conc) > 1E-5:
               indices = np.where(cloud_conc>1E-5)
               cld_base = np.nanmin(self.model['zh'][indices])
               if zindex == 0:
                  print(cloud_conc[0:44])
                  print(cld_base)
            else:
               cld_base = 1E100
         else:
            cld_base = 1E100    
 
         #--- Next Make Sure Reflectivity is Low (If Desired)
         if refl_limit or cloud_base_limit: #--- Perfomed in refl_limit or cloud_base_limit called
            refl_vertical = np.zeros(self.model['zh'].shape)
            refl_threshold = 5.
            for hindex, hgt in enumerate(self.model['zh']):
               if hindex > 0: refl_vertical[hindex] = pydart.interp.point_interp(self.model,'dbz',xloc,yloc,hgt)
            if np.nanmax(refl_vertical) > refl_threshold:
               indices = np.where(refl_vertical>refl_threshold)
               refl_base = np.nanmin(self.model['zh'][indices])
               if zindex == 0:
                  print("Don't Assimilate Observations Above ... ",refl_base)
            else:
               refl_base = 1E100   
         else:
            refl_base = 1E100

         #----------------------------------#
         #   End the thresholding section   #
         #----------------------------------#

         if zloc > cld_base or zloc > refl_base:
            self.ob[zindex] = np.nan
         else:
            if varname.upper() in ['RADIOSONDE_U_WIND_COMPONENT','U_WIND_10M','DROPSONDE_U_WIND_COMPONENT']:
              self.ob[zindex] = pydart.interp.point_interp(self.model,'u',xloc,yloc,zloc)

            elif varname.upper() in ['RADIOSONDE_V_WIND_COMPONENT','V_WIND_10M','DROPSONDE_V_WIND_COMPONENT']:
               self.ob[zindex] = pydart.interp.point_interp(self.model,'v',xloc,yloc,zloc)    

            elif varname.upper() in ['RADIOSONDE_TEMPERATURE','TEMPERATURE_2M','DROPSONDE_TEMPERATURE']:
               #self.model['T'] = pydart.fwdop.theta_to_temp(self.model['pt'],self.model['p'])
               #self.ob[zindex] = pydart.interp.point_interp(self.model,'T',xloc,yloc,zloc)
               p = pydart.interp.point_interp(self.model,'p',xloc,yloc,zloc)
               pt = pydart.interp.point_interp(self.model,'pt',xloc,yloc,zloc)
               self.ob[zindex] = pydart.fwdop.theta_to_temp(pt,p)

            elif varname.upper() in ['RADIOSONDE_SURFACE_PRESSURE','SURFACE_PRESSURE','DROPSONDE_SURFACE_PRESSURE']:
               self.ob[zindex] = pydart.interp.point_interp(self.model,'p',xloc,yloc,zloc)

            elif varname.upper() in ['RADIOSONDE_SPECIFIC_HUMIDITY','SPECIFIC_HUMIDITY_2M','DROPSONDE_SPECIFIC_HUMIDITY']: #--- JDL Does CM1 use qv or specific humidity?
               #self.model['hum'] = pydart.fwdop.qv_to_spechum(self.model['qv'])
               #self.ob[zindex] = pydart.interp.point_interp(self.model,'hum',xloc,yloc,zloc)
               qv = pydart.interp.point_interp(self.model,'qv',xloc,yloc,zloc)
               self.ob[zindex] =pydart.fwdop.qv_to_spechum(qv) 

            elif varname.upper() in ['RADIOSONDE_DEWPOINT','DROPSONDE_DEWPOINT']: #--- JDL Does CM1 use qv or specific humidity?
               qv = pydart.interp.point_interp(self.model,'qv',xloc,yloc,zloc)
               p = pydart.interp.point_interp(self.model,'p',xloc,yloc,zloc)
               self.ob[zindex] =pydart.fwdop.cal_td(qv,p)
          
            elif varname.upper() in ['RADIOSONDE_RELATIVE_HUMIDITY','DROPSONDE_RELATIVE_HUMIDITY']: #--- JDL Does CM1 use qv or specific humidity?
               qv = pydart.interp.point_interp(self.model,'qv',xloc,yloc,zloc)
               p  = pydart.interp.point_interp(self.model,'p',xloc,yloc,zloc)
               pt = pydart.interp.point_interp(self.model,'pt',xloc,yloc,zloc)
               self.ob[zindex] =pydart.fwdop.cal_rh(qv,p,pt)

            else:
               print('Observation Unknown: %s'%varname)


   def radar_operator(self,dbz_only=False):
      """
      Observation operator given model state variables
    
      Optional Arguments:
         dbz_only: A string of the observation to generatee
      """
      if dbz_only:
         self.obdbz = pydart.fwdop.calcHx_fast(self.model,self.obx,self.oby,self.obz,self.elv,self.az,
                                          self.clearair,dbzOnly=True)
      else:
         self.obvr,self.obdbz = pydart.fwdop.calcHx_fast(self.model,self.obx,self.oby,self.obz,self.elv,self.az,
                                                    self.clearair) 




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
      #--- The nskip-1 helps to start at 2 km
      self.obx = self.obx[:,nskip-1::nskip,nskip-1::nskip]
      self.oby = self.oby[:,nskip-1::nskip,nskip-1::nskip]
      self.obz = self.obz[:,nskip-1::nskip,nskip-1::nskip]
      self.az  =  self.az[:,nskip-1::nskip,nskip-1::nskip]
      self.elv = self.elv[:,nskip-1::nskip,nskip-1::nskip]

      #--- Perform cressman interpolation for each radar tilt
      self.obdbz = pydart.interp.cressman(fine_x,fine_y,fine_dbz,self.obx,self.oby,roi)
      self.obvr  = pydart.interp.cressman(fine_x,fine_y,fine_vr,self.obx,self.oby,roi)


   def addob(self,obname,error,obtype=None,obdbz=False,obvr=False,seed=None):
      """ 
      Required Inputs:
         obname - The name of the observation
         error  - The observation error variance
 
      Optional Inputs include:
         obtype - The name of the observation platform (e.g., RADAR_REFLECTIVITY)
         obdbz -  Working with Reflectivity Observations
         obvr -   Working with Radial Velocity Observations
         seed -   A seed to generate random perturbations 
 
      The saved variables in a nested dictionary array include: 
         obname   :  The name the the observations will be stored under 
                     (same for one type of platform)
         truth    :  The observed values
         obs      :  The observed values (with noise)
         xloc     :  The x-location of the observation
         yloc     :  The y-location of the observation
         zloc     :  The z-location of the observation
         error    :  The observation error variance
         obtype   :  The DART code for the observation (int)
         missing_flag   : A flag to tell DART what observations to ignore (float)

         Radar Exclusive Variables:
         elevation     :  The radar tilt of the radar    (if needed)
         azimuth       :  The azimuth angle of the radar (if needed)
         clear_air     :  The threshold for clear air obs


      All Variables should have the same dimension except missing_flag and obtype

      The structure of the nested dictionary is as follows:
         dictionary['OBS_PLATFORM']['OBSNAME']['obsx','obsy','obsz'...]

      """
      #--- Saving the observation platform location
      self.obs[self.plat_name][obname]  = {}
      if obtype is None:
         obcode = pydart.readvar.obcode()
         try:
            obtype = obcode[obname]
         except:
            print('Warning %s Does Not Exist... Setting to -1'%obname)
            obtype = -1

      if seed is None:
         pass
      else:
         np.random.seed(seed)

      #--- Saving Captured Observations (special excpetions for dbz/vr)
      if obdbz :
         self.obs[self.plat_name][obname]['truth']   = self.obdbz
         #--- Adding random errors to observations
         obnoise = np.random.normal(loc=0.0,scale=np.sqrt(np.mean(error)),size=self.obdbz.shape)
         noisy_obs = self.obdbz + obnoise
         noisy_obs[noisy_obs<0] = 0
         self.obs[self.plat_name][obname]['obs']     = noisy_obs 
      elif obvr:
         self.obs[self.plat_name][obname]['truth']   = self.obvr
         #--- Add random errors to Obs
         obnoise = np.random.normal(loc=0.0,scale=np.sqrt(np.mean(error)),size=self.obvr.shape)
         self.obs[self.plat_name][obname]['obs']     = self.obvr + obnoise
      else:
         self.obs[self.plat_name][obname]['truth']     = self.ob
         #--- Add random rrrors to Obs
         obnoise = np.random.normal(loc=0.0,scale=np.sqrt(np.mean(error)),size=self.ob.shape)
         self.obs[self.plat_name][obname]['obs']     = self.ob + obnoise

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
   lines = open('./data/DART_obs.csv','r')
   obcode = {}
   for lindex, line in enumerate(lines):
      if lindex > 0:
         ln = line.split(', ')
         obcode[ln[0]] = int(ln[-1])
   return obcode
