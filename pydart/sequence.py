#--- This module is used to read and output obs_sequence files
import os
import numpy as np
class create_sequence():
   def __init__(self,filepath):
      try:
         os.system('rm %s'%filepath)
      except:
         print('obs_sequence file does not exist')
      self.txtfile = open(filepath,"w")
   

   def add_header(self,obname,obcode,nob):
      """
      This function is used to generate the text header 
      in the observation sequence file 

      Example Header Output - 

      obsequence
         2
         1 RADIOSONDE_U_WIND_COMPONENT
         42 RADAR_REFLECTIVITY   
         num_copies:            2  num_qc:            1
         num_obs:            3  max_num_obs:            3
      observations
      truth
      TEST

      Required inputs:
         txtfile - An opened ascii text file
         obname  - The listed available observations
         obcode  - An array of the DART observation codes
         nob     - The total number of observations (assimilated and flagged)

      """
      nobtype = len(obname)
      self.nob = nob
      self.obnames = obname
      self.txtfile.write("obs_sequence\n")
      self.txtfile.write("obs_kind_definitions\n")
      self.txtfile.write("           %s\n"%len(obname))
   
      print('The number of obs =',nobtype)

      #--- Writing out the assimilated ob types   
      for index in range(0,nobtype):
         self.txtfile.write("          %s %s\n"%(obcode[index],obname[index]))

      #--- Bottom half 
      self.txtfile.write("  num_copies:            %d  num_qc:            %d\n"%(2,1))  #--- Don't know exactly the purpose of this line
      self.txtfile.write("  num_obs:            %d  max_num_obs:            %d\n"%(nob,nob))
      self.txtfile.write("observations\n")
      self.txtfile.write("truth\n")
      self.txtfile.write("Platform Number\n") #--- Reference Number for Ob Platform
      self.txtfile.write("Quality Control\n")    #--- Quality Control
      self.txtfile.write("  first:            %d  last:       %d\n"%(1,nob))

   def add_entries(self,observation):
      """ 
      This function is used to add an observation entry to the obs_sequence file. 
   
      Required Entries 
         txtfile     -  The ASCII txtfile used to print out info
         observation -  A nested dictionary created by the observation class used
                        to provide observation information
      """ 
      #--- Radial Wind Observaitions Require Special Output
      #--- JDL You May need to make a seeded random number generator
      nob = 0
      for platform  in observation.keys():
         for obname in observation[platform].keys():
            if obname in self.obnames:
                observations = observation[platform][obname]['obs'].flatten() 
                xloc = observation[platform][obname]['xloc'].flatten()
                yloc = observation[platform][obname]['yloc'].flatten()
                zloc = observation[platform][obname]['zloc'].flatten()
                obtype = observation[platform][obname]['obstype'].flatten()
                error = observation[platform][obname]['error'].flatten()
                missing = observation[platform][obname]['missing_flag'] #---  Value for missing flag
                for oindex,ob in enumerate(observations):
                   if np.isnan(ob): 
                      ob = 0
                      obx= 0
                      oby = 0
                      obz = 0
                      flag = missing
                   else:
                      obx = xloc[oindex]
                      oby = yloc[oindex]
                      obz = zloc[oindex]  
                      flag = 0.0

                   self.txtfile.write(" OBS            %d\n"%nob)
                   self.txtfile.write(" %d\n"%ob) #--- JDK One of these includes errors find out which one, and how to treat
                   self.txtfile.write(" %d\n"%ob) #--- JDL One of these includes errors find out which one, and how to treat this
                   self.txtfile.write(" %d\n"%int(platform[-3:])) #--- Platform Number e.g., RADAR_002 = 002
                   self.txtfile.write(" %d\n"%flag) #--- Quality Control Flag     
                   #--- ADDING THREE LOCATION MARKER
                   if nob == 0:
                      self.txtfile.write("      %d          %d          %d\n"%(-1,2,-1)) 
                   elif nob == (self.nob - 1): 
                      self.txtfile.write("      %d          %d          %d\n"%(nob,-1,-1))
                   else:
                      self.txtfile.write("      %d          %d          %d\n"%(nob,nob+1,-1))
                   self.txtfile.write("obdef\n")
                   self.txtfile.write("loc3d\n")
                   self.txtfile.write("     %d        %d         %d      3\n"%(obx,oby,obz))
                   self.txtfile.write("kind\n")
                   self.txtfile.write("          %d\n"%int(obtype[oindex])) 
                   nob += 1    
            #if observation[platform][obname]['obstype'] == 36:
      #else:
      self.txtfile.close()





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
