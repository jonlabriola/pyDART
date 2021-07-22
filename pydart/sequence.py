#--- This module is used to read and output obs_sequence files
import os
import numpy as np
import datetime as dt
import pydart


class create_sequence():
   def __init__(self,filepath):
      try:
         os.system('rm %s'%filepath)
      except:
         print('obs_sequence file does not exist')
      self.txtfile = open(filepath,"w")
 
      #--- A list of all available DART Observations
      self.obs_codes = pydart.readvar.obcode() #--- Listed Observations
      self.dart_obs  = self.obs_codes.keys()   #--- Observation Names

      #--- Obs information
      self.nob = None
      self.obname = None
      self.obcode = None 

 
   def obinfo(self,observations):
      """
      Get observation information for the sequence file including the number of
      assimilated observations
   
      Required input:
         Observations - A dictionary of the different assimilated observations

      """

      nobs = 0
      self.obname = []
      self.obcode = []

      for platform in observations.keys():          #--- Loop over platforms
         for key in observations[platform]:         #--- Loop over saved values
            if key in self.dart_obs:
               nobs += np.size(observations[platform][key]['obs'])
               if key in self.obname:
                  pass
               else:
                  self.obname.append(key)
                  self.obcode.append(self.obs_codes[key]) 
      self.nob = nobs

   def add_header(self,observations):#,obname,obcode):
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
         self.obname  - The listed available observations
         self.obcode  - An array of the DART observation codes
         self.nob     - The total number of observations (assimilated and flagged)

      """
      #--- Gather Observation Information If not Available
      if self.obname is None:
         self.obinfo(observations)


      nobtype = len(self.obname)
      self.txtfile.write("obs_sequence\n")
      self.txtfile.write("obs_kind_definitions\n")
      self.txtfile.write("           %s\n"%nobtype)
   
      print('The number of obs =',nobtype)

      #--- Writing out the assimilated ob types   
      for index in range(0,nobtype):
         self.txtfile.write("          %s %s\n"%(self.obcode[index],self.obname[index]))

      #--- Bottom half 
      self.txtfile.write("  num_copies:            %d  num_qc:            %d\n"%(2,1))  #--- Don't know exactly the purpose of this line
      self.txtfile.write("  num_obs:            %d  max_num_obs:            %d\n"%(self.nob,self.nob))
      self.txtfile.write("observations\n")
      self.txtfile.write("truth\n")
      self.txtfile.write("Platform Number\n") #--- Reference Number for Ob Platform
      self.txtfile.write("Quality Control\n")    #--- Quality Control
      self.txtfile.write("  first:            %d  last:       %d\n"%(1,self.nob))

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
      obnum = 0
      init_time = dt.datetime(1601,1,1,0,0,0)
      vr_code = self.obs_codes['DOPPLER_RADIAL_VELOCITY']
      for platform  in observation.keys():
         #--- Grab Observation Date, Calculate Epoch Time
         date = observation[platform]['date']
         curr_time = dt.datetime(date['year'],date['month'],date['day'],date['hour'],date['minute'],date['second'])
         delta = curr_time - init_time
         days = int(delta.days)
         seconds = int(delta.seconds)

         for obtype in observation[platform].keys():
            if obtype in self.obname:
                #--- Gathering Information For Each Observation
                missing = observation[platform][obtype]['missing_flag']
                obs = observation[platform][obtype]['obs'] 
                obs[obs==np.nan] = missing


                tmp  = observation[platform][obtype] 
                missing = tmp['missing_flag']
                obcode = tmp['obstype'] 
                #--- Flatten 3D Arrays 
                obs    = tmp['obs'].flatten()
                flag   = np.where(np.isnan(obs),missing,0.)
                xloc   = tmp['xloc'].flatten()
                yloc   = tmp['yloc'].flatten()
                zloc   = tmp['zloc'].flatten()
                error  = tmp['error'].flatten()
                flag   = np.where(np.isnan(obs),missing,0.)
             
                #--- Remove NaN's
                indices = np.where(np.isnan(obs))
                flag   = np.where(np.isnan(obs),missing,0.)
                obs[indices]  = missing
                xloc[indices] = missing
                yloc[indices] = missing
                zloc[indices] = missing               

                if obcode == self.obs_codes['DOPPLER_RADIAL_VELOCITY']:
                   nyquist =  observation[platform]['nyquist']
                   azimuth = tmp['azimuth'].flatten()
                   elevation = tmp['elevation'].flatten()
                   azimuth[indices] = missing
                   elevation[indices] = missing

                for oindex,ob in enumerate(obs):

                   #--- Writing Out for Each Observations
                   self.txtfile.write(" OBS            %d\n"%obnum)
                   self.txtfile.write(" %20.14f\n"%ob) #--- JDL One of these includes errors find out which one, and how to treat
                   self.txtfile.write(" %20.14f\n"%ob) #--- JDL One of these includes errors find out which one, and how to treat this
                   self.txtfile.write(" %d\n"%int(platform[-3:])) #--- Platform Number e.g., RADAR_002 = 002
                   self.txtfile.write(" %d\n"%flag[oindex])               #--- Quality Control Flag     
                   #--- ADDING THREE LOCATION MARKER
                   if obnum == 0:                 #--- First Observation
                      self.txtfile.write("      %d          %d          %d\n"%(-1,2,-1)) 
                   elif obnum == (self.nob - 1):  #--- Final Observation
                      self.txtfile.write("      %d          %d          %d\n"%(obnum,-1,-1))
                   else:                          #--- Every Other Observation
                      self.txtfile.write("      %d          %d          %d\n"%(obnum,obnum+1,-1))
                   self.txtfile.write("obdef\n")
                   self.txtfile.write("loc3d\n")
                   self.txtfile.write("    %20.14f          %20.14f        %20.14f    3\n"%(xloc[oindex],yloc[oindex],zloc[oindex])) #--- Obs Location
                   self.txtfile.write("kind\n")
                   self.txtfile.write("          %d\n"%int(obcode))                  #--- Obs Code

                   #--- Add Special input for Radial Velocity Observations
                   if obcode == self.obs_codes['DOPPLER_RADIAL_VELOCITY']: 
                      #--- Grabbed From DART Code
                      #beam_direction(1) = sin(az) * cos(el)
                      #beam_direction(2) = cos(az) * cos(el)
                      #beam_direction(3) = sin(el)
                      beam_dir1 = np.sin(azimuth[oindex]) * np.cos(elevation[oindex])
                      beam_dir2 = np.cos(azimuth[oindex]) * np.cos(elevation[oindex])
                      beam_dir3 = np.sin(elevation[oindex])
                      self.txtfile.write("    %20.14f          %20.14f        %20.14f\n" % (beam_dir1, beam_dir2, beam_dir3) ) 
                      self.txtfile.write("    %20.14f     \n" %nyquist)


                   self.txtfile.write("    %d          %d     \n" % (seconds, days))         #--- Obs Date


                   obnum += 1    







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
