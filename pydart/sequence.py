#--- This module is used to read and output obs_sequence files
import os
import numpy as np
import datetime as dt
import pydart

class create_sequence():
   """
   This class is used to generate an observation sequence file given observation classes that are generated 
   via pyDart

   Require Inputs: 
      filepath - The path of the obs_sequence file
      observaitions - the observation class

   Output:
      An ascii text file that emulates obs_sequence files created by DART 
   """

   def __init__(self,filepath,observations):
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

      self.add_header(observations)
      self.add_entries(observations)
 
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
      self.txtfile.write("  num_copies:            %d  num_qc:            %d\n"%(2,2))  #--- May need to change num_qc if number of flags changese
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
      vr_ob_count = 1 #--- Counting the number of vr observations
      for platform  in observation.keys():
         #--- Grab Observation Date, Calculate Epoch Time
         date = observation[platform]['date']
         curr_time = dt.datetime(date['year'],date['month'],date['day'],date['hour'],date['minute'],date['second'])
         delta = curr_time - init_time
         days = int(delta.days)
         seconds = int(delta.seconds)

         #--- Location
         platform_x = observation[platform]['plat_x'] 
         platform_y = observation[platform]['plat_y'] 
         platform_z = observation[platform]['plat_z'] 
            

         for obtype in observation[platform].keys():
            if obtype in self.obname:
                #--- Gathering Information For Each Observation
                missing = observation[platform][obtype]['missing_flag']
                obs = observation[platform][obtype]['obs'] 
                obs[obs==np.nan] = missing


                tmp  = observation[platform][obtype] 
                missing = tmp['missing_flag']
                obcode = tmp['obstype'] 

                #--- JDL Hard Code to account for changes in ob code
                #--- between DART versions
                if obcode == 41: obcode = 45 #--- Switch Vr Code
                if obcode == 42: obcode = 46 #--- Switch REFL Code

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
                   self.txtfile.write(" OBS            %d\n"%(obnum+1))
                   self.txtfile.write(" %20.14E\n"%ob) #--- JDL One of these includes errors find out which one, and how to treat
                   self.txtfile.write(" %20.14E\n"%ob) #--- JDL One of these includes errors find out which one, and how to treat this
                   self.txtfile.write(" %20.14E\n"%int(platform[-3:])) #--- Platform Number e.g., RADAR_002 = 002
                   #--- JDL This is done to avoid assimilating any observations on the edge boundaries (Results in a DART Error)
                   #--- This is currently hard coded for the QLCS Case 
                   if xloc[oindex] < 1000. or xloc[oindex] > 199000. or yloc[oindex] < 1000. or yloc[oindex] > 199000. or zloc[oindex] < 150. or zloc[oindex] >15000:
                      self.txtfile.write(" %20.14E\n"%np.absolute(99999.))
                   #--- JDL This is what you would do other wise
                   else:
                      self.txtfile.write(" %20.14E\n"%np.absolute(flag[oindex]))              #--- Quality Control Flag     
                   #--- ADDING THREE LOCATION MARKER
                   if self.nob == 1:
                      self.txtfile.write("      %d          %d          %d\n"%(-1,-1,-1))
                   elif obnum == 0:                 #--- First Observation
                      self.txtfile.write("      %d          %d          %d\n"%(-1,2,-1)) 
                   elif obnum == (self.nob - 1):  #--- Final Observation
                      self.txtfile.write("      %d          %d          %d\n"%(obnum,-1,-1))
                   else:                          #--- Every Other Observation
                      self.txtfile.write("      %d          %d          %d\n"%(obnum,obnum+2,-1))
                   self.txtfile.write("obdef\n")
                   self.txtfile.write("loc3Dxyz\n")
                   self.txtfile.write("    %20.14f          %20.14f        %20.14f    \n"%(xloc[oindex],yloc[oindex],zloc[oindex])) #--- Obs Location
                   self.txtfile.write("kind\n")
                   self.txtfile.write("          %d  \n"%int(obcode))                  #--- Obs Code

                   #--- Add Special input for Radial Velocity Observations
                   if obcode == self.obs_codes['DOPPLER_RADIAL_VELOCITY']: 
                      #--- Grabbed From DART Code
                      #beam_direction(1) = sin(az) * cos(el)
                      #beam_direction(2) = cos(az) * cos(el)
                      #beam_direction(3) = sin(el)

                      self.txtfile.write("platform\n")
                      self.txtfile.write("loc3Dxyz\n")
                      self.txtfile.write("    %20.14f          %20.14f        %20.14f\n" % (platform_x,platform_y,platform_z) )
                      self.txtfile.write("dir3d\n")
                      beam_dir1 = np.sin(azimuth[oindex]) * np.cos(elevation[oindex])
                      beam_dir2 = np.cos(azimuth[oindex]) * np.cos(elevation[oindex])
                      beam_dir3 = np.sin(elevation[oindex])
                      self.txtfile.write("    %20.14f          %20.14f        %20.14f\n" % (beam_dir1, beam_dir2, beam_dir3) ) 
                      self.txtfile.write("    %20.14f     \n" %nyquist)
                      self.txtfile.write("    %d     \n" %int(vr_ob_count)) # JDL - For some reason new DART version requires vr to save vr ob count
                      vr_ob_count += 1


                   self.txtfile.write("    %d          %d     \n" % (seconds, days))         #--- Obs Date
                   self.txtfile.write("    %20.14f    \n" % (error[oindex]))         #--- Obs Error

                   obnum += 1    

      self.txtfile.close()


class read_sequence():
   """
   This class is used to read in observation sequence files 
   """
   def __init__(self,filepath):
      self.read_header(filepath)
      self.read_radObs(filepath)

   def read_header(self,filepath):
      """
      The goal of this function is to read an ascii output file and save the relevant information into an observation class

      Required inputs: Filepath
      """
      lines = open(filepath)
      self.obcode = []
      self.obname = []

      self.copy = []
      self.qc = []

      self.obindex =[]

      #--- Adding Holders
      num_ob_types = 999.
      copy_start = 999.
      copy_start = 999.
      qc_start   = 999.
      qc_stop    = 999.

     
      for lindex, line in enumerate(lines):
         if lindex == 2: 
            #--- Number of Observation Types
            num_ob_types = int(line)

         elif lindex > 2 and lindex <= 2+num_ob_types:
            #--- Read the observation  types and there coresponding codes
            ln = line.split()
            self.obcode.append(int(ln[0]))
            self.obname.append(ln[1])

         elif lindex == 3+num_ob_types:
            #--- Number of Ob copies/QC Flags
            ln = line.split()
            self.ncopies = int(ln[1])
            self.nqc = int(ln[3])

            #--- Updating the boundary limits
            copy_start = 4+num_ob_types
            copy_stop = copy_start+self.ncopies
            qc_start = copy_stop
            qc_stop = qc_start+self.nqc

         elif lindex == 4+num_ob_types:
            #--- Read in the number of obs 
            ln = line.split()
            self.nob = int(ln[1])  

         elif lindex > copy_start and lindex <= copy_stop:
            #--- Read in all the observation copies
            self.copy.append(line)            

         elif lindex > qc_start and lindex <= qc_stop:
            #--- Read in the QC flags
            self.qc.append(line)
         
         elif lindex > qc_stop:
            #--- Break Statement to save time
            if 'OBS' in line:
               self.obindex.append(lindex)
            #self.obstart = lindex+1
            #break



   def read_radObs(self,filepath):
      """
      How to read observation types
      
      Should I read by observation type?
      """
      lines = open(filepath)
      test = lines.readlines()
      plat_index   = len(self.copy) + 1    #--- Observation Platform
      flag_index   = len(self.copy) + 2    #--- JDL Flag   
      dartqc_index = len(self.copy) + 3    #--- DART QC Flag
      loc_index    = len(self.qc) + len(self.copy) + 4 #--- Location of the Observation locations
      kindex       = len(self.qc) + len(self.copy) + 6 #--- Observation Kind


      nobs = 1
      nplat = 1
      nobs_radar = []
      dbz_ob_index = []
      vr_ob_index = []

      #---Step 1 Get The Number of Observations For Each Radar
      for index in self.obindex:
         platform = float(test[plat_index+index].split()[0])
         if int(test[kindex+index]) == 45:
            dbz_ob_index.append(index)
            if platform > nplat:
               nobs_radar.append(nobs) 
               nobs = 1
               nplat = platform
            else:
               nobs += 1

         elif int(test[kindex+index]) == 46:
            vr_ob_index.append(index)

      #--- Add Final Platform
      nobs_radar.append(nobs)

      #--- Step 2: Create the Radar Dictionary
      obs = {}
      print('JDL numplat = ',nobs_radar)
      for nrdr,nob in enumerate(nobs_radar):
         plat_name = 'radar_%03d'%(nrdr+1)
         obs[plat_name]={}
         #print('Number of obs = ',nobs_radar)
         #--------------------------#
         #--- Radar Reflectivity ---#
         #--------------------------#
         obs[plat_name]['RADAR_REFLECTIVITY'] = {}
         obs[plat_name]['RADAR_REFLECTIVITY']['obx']     = np.zeros((nob))
         obs[plat_name]['RADAR_REFLECTIVITY']['oby']     = np.zeros((nob))
         obs[plat_name]['RADAR_REFLECTIVITY']['obz']     = np.zeros((nob))
         obs[plat_name]['RADAR_REFLECTIVITY']['flag']    = np.zeros((nob))  
         obs[plat_name]['RADAR_REFLECTIVITY']['DART_QC'] = np.zeros((nob))
         for cindex,copy in enumerate(self.copy):
            #--- Loop Over the Different Copies
            obs[plat_name]['RADAR_REFLECTIVITY'][copy] = np.zeros((nob)) 
            obcount = 0 
            for oindex in dbz_ob_index:
               #--- Make sure platform number is the same as radar number
               if float(test[oindex+plat_index]) == nrdr+1: 
                  ob_loc = test[oindex+loc_index].split()
                  obs[plat_name]['RADAR_REFLECTIVITY'][copy][obcount]      = float(test[oindex+cindex+1])
                  if cindex == 0:
                     obs[plat_name]['RADAR_REFLECTIVITY']['obx'][obcount]     = float(ob_loc[0])
                     obs[plat_name]['RADAR_REFLECTIVITY']['oby'][obcount]     = float(ob_loc[1])
                     obs[plat_name]['RADAR_REFLECTIVITY']['obz'][obcount]     = float(ob_loc[2]) 
                     obs[plat_name]['RADAR_REFLECTIVITY']['flag'][obcount]    = float(test[oindex+flag_index])
                     obs[plat_name]['RADAR_REFLECTIVITY']['DART_QC'][obcount] = float(test[oindex+dartqc_index])
                  obcount += 1

         #--------------------------#
         #---  Radial Velocity   ---#
         #--------------------------#
         obs[plat_name]['DOPPLER_RADIAL_VELOCITY'] = {}
         obs[plat_name]['DOPPLER_RADIAL_VELOCITY']['obx']     = np.zeros((nob))
         obs[plat_name]['DOPPLER_RADIAL_VELOCITY']['oby']     = np.zeros((nob))
         obs[plat_name]['DOPPLER_RADIAL_VELOCITY']['obz']     = np.zeros((nob))
         obs[plat_name]['DOPPLER_RADIAL_VELOCITY']['flag']    = np.zeros((nob))  
         obs[plat_name]['DOPPLER_RADIAL_VELOCITY']['DART_QC'] = np.zeros((nob))
         for cindex,copy in enumerate(self.copy):
            #--- Loop Over the Different Copies
            obs[plat_name]['DOPPLER_RADIAL_VELOCITY'][copy] = np.zeros((nob)) 
            obcount = 0 
            for oindex in vr_ob_index:
               #--- Make sure platform number is the same as radar number
               if float(test[oindex+plat_index]) == nrdr+1:       
                  ob_loc = test[oindex+loc_index].split()
                  obs[plat_name]['DOPPLER_RADIAL_VELOCITY'][copy][obcount]      = float(test[oindex+cindex+1])
                  if cindex == 0:
                     obs[plat_name]['DOPPLER_RADIAL_VELOCITY']['obx'][obcount]     = float(ob_loc[0])
                     obs[plat_name]['DOPPLER_RADIAL_VELOCITY']['oby'][obcount]     = float(ob_loc[1])
                     obs[plat_name]['DOPPLER_RADIAL_VELOCITY']['obz'][obcount]     = float(ob_loc[2]) 
                     obs[plat_name]['DOPPLER_RADIAL_VELOCITY']['flag'][obcount]    = float(test[oindex+flag_index])
                     obs[plat_name]['DOPPLER_RADIAL_VELOCITY']['DART_QC'][obcount] = float(test[oindex+dartqc_index])
                  obcount += 1


      #---- JDL Final Addition - Need to reshape the different dictionary entries to occupy the space of ntilt,ny,nx
      #---- This could be done in a seperate function

      return obs
               
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
