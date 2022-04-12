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
                #if obcode == 41: obcode = 45 #--- Switch Vr Code
                #if obcode == 42: obcode = 46 #--- Switch REFL Code

                #--- Flatten 3D Arrays 
                truth  = tmp['truth'].flatten(order='F')
                print('JDL truth max = ',np.nanmax(truth))
                obs    = tmp['obs'].flatten(order='F')
                print('JDL obs max = ',np.nanmax(obs))
                print('Obcode = ',obcode)
                print('XLOC = ',tmp['xloc'].flatten(order='F'))
                print('YLOC = ',tmp['yloc'].flatten(order='F'))
                flag   = np.where(np.isnan(obs),missing,0.)
                xloc   = tmp['xloc'].flatten(order='F')
                yloc   = tmp['yloc'].flatten(order='F')
                zloc   = tmp['zloc'].flatten(order='F')
                error  = tmp['error'].flatten(order='F')
                #--- JDL SWITCH ERROR
                #error[error == 6.0] = 5.0
                #error[error == 4.0] = 3.0
                flag   = np.where(np.isnan(obs),missing,0.)
             
                #--- Remove NaN's
                indices = np.where(np.isnan(obs))
                flag   = np.where(np.isnan(obs),missing,0.)
                obs[indices]    = missing
                truth[indices]  = missing
                xloc[indices]   = missing
                yloc[indices]   = missing
                zloc[indices]   = missing               

                if obcode == self.obs_codes['DOPPLER_RADIAL_VELOCITY']:
                   nyquist =  observation[platform]['nyquist']
                   azimuth = tmp['azimuth'].flatten(order='F')
                   elevation = tmp['elevation'].flatten(order='F')
                   azimuth[indices] = missing
                   elevation[indices] = missing

                for oindex,ob in enumerate(obs):

                   #--- Writing Out for Each Observations
                   self.txtfile.write(" OBS            %d\n"%(obnum+1))
                   self.txtfile.write(" %20.14E\n"%ob) #--- JDL One of these includes errors find out which one, and how to treat
                   self.txtfile.write(" %20.14E\n"%truth[oindex]) #--- JDL One of these includes errors find out which one, and how to treat this
                   self.txtfile.write(" %20.14E\n"%int(platform[-3:])) #--- Platform Number e.g., RADAR_002 = 002
                   #--- JDL This is done to avoid assimilating any observations on the edge boundaries (Results in a DART Error)
                   #--- This is currently hard coded for the QLCS Case 
                   #if xloc[oindex] < 1000. or xloc[oindex] > 199000. or yloc[oindex] < 1000. or yloc[oindex] > 199000. or zloc[oindex] < 0. or zloc[oindex] >14000:
                   #if xloc[oindex] < 1500. or xloc[oindex] > 196500. or yloc[oindex] < 1500. or yloc[oindex] > 196500. or zloc[oindex] < 0. or zloc[oindex] >14000:
                   #if xloc[oindex] < 3000. or xloc[oindex] > 194500. or yloc[oindex] < 3000. or yloc[oindex] > 194500. or zloc[oindex] < 0. or zloc[oindex] >14000:
                   if xloc[oindex] < 3000. or xloc[oindex] > 597000. or yloc[oindex] < 3000. or yloc[oindex] > 297000. or zloc[oindex] < 0. or zloc[oindex] > 14000:
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
   def __init__(self,filepath,ntilt=None,ny=None,nx=None,radar=False,conv=False):
      self.read_header(filepath)
      self.obs = {}
      self.missing = 99999.
      if radar: #--- Read in Radar Observations
         self.read_rad_obs(filepath)
         if ntilt is None or ny is None or nx is None:
            pass
         else:
            self.reshape_rad_obs(ntilt,ny,nx)
      elif conv:
        print('Reading Profilers')
        self.read_conv_obs(filepath)
        print('Reading Surface')
        self.read_conv_obs(filepath,sfc=True) #--- Read in Surface Types
      else:
        print('You are reading in neither radar nor conventional observations')


   def read_header(self,filepath):
      """
      The goal of this function is to read an ascii output file and save the relevant information into an observation class

      Required inputs: Filepath
      """
      lines = open(filepath)
      self.obcode = []
      self.obname = []

      self.copy_name = []
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
            #--- Read in all the observation copy types (e.g., mem0X_bckd, mem0X_analy)
            self.copy_name.append(line)            

         elif lindex > qc_start and lindex <= qc_stop:
            #--- Read in the QC flags
            self.qc.append(line)
         
         elif lindex > qc_stop:
            #---Read in observation indiced
            if 'OBS' in line:
               self.obindex.append(lindex)

   def read_rad_obs(self,filepath):
      """
      A function used to read radar observations from obs_seq file

      Required Arguments:
         filepath (string) - The path to the obs_sequence file

      """

      #--- Read in the sequence file
      #--- Define the important indices for each observation
      lines = open(filepath)
      rad_obs = lines.readlines()
      plat_index   = len(self.copy_name) + 1    #--- Observation Platform
      flag_index   = len(self.copy_name) + 2    #--- JDL Flag   
      dartqc_index = len(self.copy_name) + 3    #--- DART QC Flag
      loc_index    = len(self.qc) + len(self.copy_name) + 4 #--- Location of the Observation locations
      kindex       = len(self.qc) + len(self.copy_name) + 6 #--- Observation Kind


      nobs = 0
      nplat = 1
      nobs_radar = []
      dbz_ob_index = []
      vr_ob_index = []

      #--- Step 1: Get The Number of Observations For Each Radar
      #--- Save the location where the dbz/vr observations can
      #--- be found
   
      #--- nobs_radar = The number of either dbz / vr obs per radar
      #--- dbz_ob_index/vr_ob_index = Indices where radar obs start
      for index in self.obindex:
         platform = float(rad_obs[index+plat_index].split()[0])
         if int(rad_obs[kindex+index]) == 46:
            dbz_ob_index.append(index)
            if platform > nplat:
               nobs_radar.append(nobs) 
               nobs = 1
               nplat = platform
            else:
               nobs += 1

         elif int(rad_obs[kindex+index]) == 45:
            vr_ob_index.append(index)

      #--- Add Final Radar Platform Info
      nobs_radar.append(nobs)


      #--- Step 2: Create the Radar Dictionary
      #--- First Loop over Each Radar
      for nrdr,nob_per_rdr in enumerate(nobs_radar):
         plat_name = 'radar_%03d'%(nrdr+1)
         self.obs[plat_name]={}

         #--------------------------#
         #--- Radar Reflectivity ---#
         #--------------------------#
         self.obs[plat_name]['RADAR_REFLECTIVITY'] = {}
         self.obs[plat_name]['RADAR_REFLECTIVITY']['xloc']     = np.zeros((nob_per_rdr))
         self.obs[plat_name]['RADAR_REFLECTIVITY']['yloc']     = np.zeros((nob_per_rdr))
         self.obs[plat_name]['RADAR_REFLECTIVITY']['zloc']     = np.zeros((nob_per_rdr))
         self.obs[plat_name]['RADAR_REFLECTIVITY']['flag']    = np.zeros((nob_per_rdr))  
         self.obs[plat_name]['RADAR_REFLECTIVITY']['DART_QC'] = np.zeros((nob_per_rdr))
         for cindex,copy in enumerate(self.copy_name):
            #--- Loop Over the Different Copy Types
            self.obs[plat_name]['RADAR_REFLECTIVITY'][copy] = np.zeros((nob_per_rdr)) 
            obcount = 0 
            for oindex in dbz_ob_index:
               #--- Make sure platform number is the same as radar number
               if float(rad_obs[oindex+plat_index]) == nrdr+1: 
                  ob_loc = rad_obs[oindex+loc_index].split()
                  flag  = float(rad_obs[oindex+flag_index])
                  if flag >= self.missing:
                     self.obs[plat_name]['RADAR_REFLECTIVITY'][copy][obcount]      = np.nan
                  else:
                     self.obs[plat_name]['RADAR_REFLECTIVITY'][copy][obcount]      = float(rad_obs[oindex+cindex+1])
                  if cindex == 0:
                     if flag >= self.missing:
                        self.obs[plat_name]['RADAR_REFLECTIVITY']['xloc'][obcount]    = np.nan
                        self.obs[plat_name]['RADAR_REFLECTIVITY']['yloc'][obcount]    = np.nan
                        self.obs[plat_name]['RADAR_REFLECTIVITY']['zloc'][obcount]    = np.nan 
                        self.obs[plat_name]['RADAR_REFLECTIVITY']['flag'][obcount]    = float(rad_obs[oindex+flag_index])
                        self.obs[plat_name]['RADAR_REFLECTIVITY']['DART_QC'][obcount] = np.nan
                     else:
                        self.obs[plat_name]['RADAR_REFLECTIVITY']['xloc'][obcount]    = float(ob_loc[0])
                        self.obs[plat_name]['RADAR_REFLECTIVITY']['yloc'][obcount]    = float(ob_loc[1])
                        self.obs[plat_name]['RADAR_REFLECTIVITY']['zloc'][obcount]    = float(ob_loc[2]) 
                        self.obs[plat_name]['RADAR_REFLECTIVITY']['flag'][obcount]    = float(rad_obs[oindex+flag_index])
                        #--- DART QC is not created until after DA so make sure the statement is option
                        try:
                           self.obs[plat_name]['RADAR_REFLECTIVITY']['DART_QC'][obcount] = float(rad_obs[oindex+dartqc_index])
                        except:
                           self.obs[plat_name]['RADAR_REFLECTIVITY']['DART_QC'][obcount] = np.nan
                  obcount += 1
         #--------------------------#
         #---  Radial Velocity   ---#
         #--------------------------#
         self.obs[plat_name]['DOPPLER_RADIAL_VELOCITY'] = {}
         self.obs[plat_name]['DOPPLER_RADIAL_VELOCITY']['xloc']     = np.zeros((nob_per_rdr))
         self.obs[plat_name]['DOPPLER_RADIAL_VELOCITY']['yloc']     = np.zeros((nob_per_rdr))
         self.obs[plat_name]['DOPPLER_RADIAL_VELOCITY']['zloc']     = np.zeros((nob_per_rdr))
         self.obs[plat_name]['DOPPLER_RADIAL_VELOCITY']['flag']    = np.zeros((nob_per_rdr))  
         self.obs[plat_name]['DOPPLER_RADIAL_VELOCITY']['DART_QC'] = np.zeros((nob_per_rdr))
         for cindex,copy in enumerate(self.copy_name):
            #--- Loop Over the Different Copies
            self.obs[plat_name]['DOPPLER_RADIAL_VELOCITY'][copy] = np.zeros((nob_per_rdr)) 
            obcount = 0 
            for oindex in vr_ob_index:
               #--- Make sure platform number is the same as radar number
               if float(rad_obs[oindex+plat_index]) == nrdr+1:       
                  ob_loc = rad_obs[oindex+loc_index].split()
                  flag  = float(rad_obs[oindex+flag_index])
                  if flag >= self.missing:
                     self.obs[plat_name]['DOPPLER_RADIAL_VELOCITY'][copy][obcount]      = np.nan
                  else:
                     self.obs[plat_name]['DOPPLER_RADIAL_VELOCITY'][copy][obcount]      = float(rad_obs[oindex+cindex+1])
 
                  if cindex == 0:
                     if flag >= self.missing:
                        self.obs[plat_name]['DOPPLER_RADIAL_VELOCITY']['xloc'][obcount]     = np.nan
                        self.obs[plat_name]['DOPPLER_RADIAL_VELOCITY']['yloc'][obcount]     = np.nan
                        self.obs[plat_name]['DOPPLER_RADIAL_VELOCITY']['zloc'][obcount]     = np.nan
                        self.obs[plat_name]['DOPPLER_RADIAL_VELOCITY']['flag'][obcount]     = float(rad_obs[oindex+flag_index])
                        #--- DART QC is not created until after DA so make sure the statement is option
                        self.obs[plat_name]['DOPPLER_RADIAL_VELOCITY']['DART_QC'][obcount]  = np.nan 
 
                     else:
                        self.obs[plat_name]['DOPPLER_RADIAL_VELOCITY']['xloc'][obcount]     = float(ob_loc[0])
                        self.obs[plat_name]['DOPPLER_RADIAL_VELOCITY']['yloc'][obcount]     = float(ob_loc[1])
                        self.obs[plat_name]['DOPPLER_RADIAL_VELOCITY']['zloc'][obcount]     = float(ob_loc[2]) 
                        self.obs[plat_name]['DOPPLER_RADIAL_VELOCITY']['flag'][obcount]    = float(rad_obs[oindex+flag_index])
                        #--- DART QC is not created until after DA so make sure the statement is option
                        try:
                           self.obs[plat_name]['DOPPLER_RADIAL_VELOCITY']['DART_QC'][obcount] = float(rad_obs[oindex+dartqc_index])
                        except:
                           self.obs[plat_name]['DOPPLER_RADIAL_VELOCITY']['DART_QC'][obcount] = np.nan 
            
                  obcount += 1



   def reshape_rad_obs(self,ntilt,ny,nx):
      """
      A function that is used to reshape all radar observation
      to fill a three-dimensional volume

      Required Arguments:
         ntilt (integer) - The number of radar tilts
         ny    (integer) - The number of grid points in the y direction
         nx    (interge) - The number of grid points in the x direction

      """ 
      for platform in self.obs.keys():
         if 'radar' in platform:
            for obtype in  self.obs[platform].keys():
               for reshape_key in self.obs[platform][obtype].keys():
                  self.obs[platform][obtype][reshape_key] = np.reshape(self.obs[platform][obtype][reshape_key],(ntilt,ny,nx),order='F')

   def read_conv_obs(self,filepath,sfc=False):
      """
      A function that is used to read in conventional observations saved at different heights

      Required arguments:
         filepath - The pathways to the sequence file that is to be read in
         sfc - Read in Surface Obs Instead of Profilers
      """  
      #--- Read in the sequence file
      #--- Define the important indices for each observation
      lines = open(filepath)
      conv_obs = lines.readlines()
      plat_index   = len(self.copy_name) + 1    #--- Observation Platform
      flag_index   = len(self.copy_name) + 2    #--- JDL Flag   
      dartqc_index = len(self.copy_name) + 3    #--- DART QC Flag
      loc_index    = len(self.qc) + len(self.copy_name) + 4 #--- Location of the Observation locations
      kindex       = len(self.qc) + len(self.copy_name) + 6 #--- Observation Kind


      nobs = 0
      nplat = 1
      
      nobs_conv     = []
      if sfc:
         conv_obnames = ['U_WIND_10M','V_WIND_10M','TEMPERATURE_2M','DEWPOINT_2_METER']
      else:
         conv_obnames = ['RADIOSONDE_U_WIND_COMPONENT','RADIOSONDE_V_WIND_COMPONENT',
                    'RADIOSONDE_TEMPERATURE','RADIOSONDE_SPECIFIC_HUMIDITY','RADIOSONDE_DEWPOINT']


      conv_uindex   = []
      conv_vindex   = []
      conv_qvindex  = []
      conv_tindex   = []
      #conv_u10index = [] #--- sfc
      #conv_v10index = [] #--- sfc
      #conv_t2index  = [] #--- sfc
      #conv_td2index = [] #--- sfc

      #conv_obnames = ['RADIOSONDE_U_WIND_COMPONENT','RADIOSONDE_V_WIND_COMPONENT',
      #           'RADIOSONDE_TEMPERATURE','RADIOSONDE_SPECIFIC_HUMIDITY','RADIOSONDE_DEWPOINT',
      #            'U_WIND_10M','V_WIND_10M','TEMPERATURE_2M','DEWPOINT_2_METER']
   
     
      obcode = pydart.readvar.obcode() #--- Listed Observations
      for key in obcode.keys():
         if sfc:
            if 'U_WIND_10M'                   in key: conv_ucode = obcode[key]       
            if 'V_WIND_10M'                   in key: conv_vcode = obcode[key]       
            if 'TEMPERATURE_2M'               in key: conv_tcode = obcode[key]       
            if 'DEWPOINT_2_METER'             in key: conv_qvcode = obcode[key]       
         else:
            if 'RADIOSONDE_U_WIND_COMPONENT'  in key: conv_ucode = obcode[key]
            if 'RADIOSONDE_V_WIND_COMPONENT'  in key: conv_vcode = obcode[key]
            if 'RADIOSONDE_TEMPERATURE'       in key: conv_tcode = obcode[key]
            if 'RADIOSONDE_SPECIFIC_HUMIDITY' in key: conv_qvcode = obcode[key]
            if 'RADIOSONDE_DEWPOINT'          in key: conv_qvcode = obcode[key]


      print('JDL UCODE = ',conv_ucode)
      print('JDL VCODE = ',conv_vcode)
      print('JDL TCODE = ',conv_tcode)
      print('JDL QVCODE = ',conv_qvcode)
      #--- Step 1: Get The Number of Observations For Each Profiler
      #--- Save the location where the u,v,T,qv observations can
      #--- be found
   
      #--- nobs_conv = The number of either dbz / vr obs per radar
      #--- conv_(u,v,qv,t)index = Indices where radar obs start

      #--- Loop through observations for the respective type
      for index in self.obindex:
         platform = float(conv_obs[index+plat_index].split()[0])
         if int(conv_obs[kindex+index]) == conv_ucode:
            conv_uindex.append(index)
            if platform > nplat:
               nobs_conv.append(nobs) 
               nobs = 1
               nplat = platform
            else:
               nobs += 1
         elif int(conv_obs[kindex+index]) == conv_vcode:
            conv_vindex.append(index)
         elif int(conv_obs[kindex+index]) == conv_qvcode:
            conv_qvindex.append(index)
         elif int(conv_obs[kindex+index]) == conv_tcode:
            conv_tindex.append(index)
   


      #--- Add Final Radar Platform Info
      nobs_conv.append(nobs)


      #--- Step 2: Create the Radar Dictionary
      #--- First Loop over Each Radar
      for nconv,nob_per_conv in enumerate(nobs_conv):
         if sfc:
            plat_name = 'sfc_%03d'%(nconv+1)
         else:
            plat_name = 'pro_%03d'%(nconv+1)
         self.obs[plat_name]={}
         #--------------------------#
         #--- LOOP OVER CONV OBS ---#
         #--------------------------#
         for obname in conv_obnames:
            self.obs[plat_name][obname] = {}
            self.obs[plat_name][obname]['xloc']     = np.zeros((nob_per_conv))
            self.obs[plat_name][obname]['yloc']     = np.zeros((nob_per_conv))
            self.obs[plat_name][obname]['zloc']     = np.zeros((nob_per_conv))
            self.obs[plat_name][obname]['flag']    = np.zeros((nob_per_conv))  
            self.obs[plat_name][obname]['DART_QC'] = np.zeros((nob_per_conv))
            for cindex,copy in enumerate(self.copy_name):
               #--- Loop Over the Different Copy Types
               self.obs[plat_name][obname][copy] = np.zeros((nob_per_conv)) 
               obcount = 0
             
               #--- U - Wind
               if 'RADIOSONDE_U_WIND_COMPONENT'    in obname: ob_index = conv_uindex
               elif 'U_WIND_10M'                   in obname: ob_index = conv_uindex
               #--- V - Wind
               elif 'RADIOSONDE_V_WIND_COMPONENT'  in obname: ob_index = conv_vindex 
               elif 'V_WIND_10M'                   in obname: ob_index = conv_vindex
               #--- Temperature
               elif 'RADIOSONDE_TEMPERATURE'       in obname: ob_index = conv_tindex
               elif 'TEMPERATURE_2M'               in obname: ob_index = conv_tindex
               #--- Moisture
               elif 'RADIOSONDE_SPECIFIC_HUMIDITY' in obname: ob_index = conv_qvindex
               elif 'RADIOSONDE_DEWPOINT' in obname: ob_index = conv_qvindex
               elif 'DEWPOINT_2_METER' in obname: ob_index = conv_qvindex

               else: print("Unspecified Observation Type")
 
               for oindex in ob_index: #--- Loop over respective indicies for each observation type
                  #--- Make sure platform number is the same as conventional number
                  if float(conv_obs[oindex+plat_index]) == nconv+1: 
                     ob_loc = conv_obs[oindex+loc_index].split()
                     flag  = float(conv_obs[oindex+flag_index])
                     if flag >= self.missing:
                        self.obs[plat_name][obname][copy][obcount]      = np.nan
                     else:
                        self.obs[plat_name][obname][copy][obcount]      = float(conv_obs[oindex+cindex+1]) #--- JDL ORIG (Probably gets truth)
                        # JDL THIS IS NOT THE LAST
                        #self.obs[plat_name][obname][copy][obcount]      = float(conv_obs[oindex+cindex])


                     if cindex == 0:
                        if flag >= self.missing:
                           self.obs[plat_name][obname]['xloc'][obcount]    = np.nan
                           self.obs[plat_name][obname]['yloc'][obcount]    = np.nan
                           self.obs[plat_name][obname]['zloc'][obcount]    = np.nan 
                           self.obs[plat_name][obname]['flag'][obcount]    = float(conv_obs[oindex+flag_index])
                           self.obs[plat_name][obname]['DART_QC'][obcount] = np.nan
                        else:
                           self.obs[plat_name][obname]['xloc'][obcount]    = float(ob_loc[0])
                           self.obs[plat_name][obname]['yloc'][obcount]    = float(ob_loc[1])
                           self.obs[plat_name][obname]['zloc'][obcount]    = float(ob_loc[2]) 
                           self.obs[plat_name][obname]['flag'][obcount]    = float(conv_obs[oindex+flag_index])
                           #--- DART QC is not created until after DA so make sure the statement is option
                           try:
                              self.obs[plat_name][obname]['DART_QC'][obcount] = float(conv_obs[oindex+dartqc_index])
                           except:
                              self.obs[plat_name][obname]['DART_QC'][obcount] = np.nan
                     obcount += 1



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
