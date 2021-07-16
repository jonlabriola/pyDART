def add_header(txtfile,obname,obcode,nob):
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

   txtfile.write("obs_sequence\n")
   txtfile.write("obs_kind_definitions\n")
   txtfile.write("           %s\n"%len(obname))
   
   print('The number of obs =',nobtype)

   #--- Writing out the assimilated ob types   
   for index in range(0,nobtype):
      txtfile.write("          %s %s\n"%(obcode[index],obname[index]))

   #--- Bottom half 
   txtfile.write("  num_copies:            %d  num_qc:            %d\n"%(2,1))  #--- Don't know exactly the purpose of this line
   txtfile.write("  num_obs:            %d  max_num_obs:            %d\n"%(nob,nob))
   txtfile.write("observations\n")
   txtfile.write("truth\n")
   txtfile.write("Quality Control\n")
   txtfile.write("  first:            %d  last:       %d\n"%(1,nob))

#def add_entry():
#   """ 
#   Add an entry to the observation Sequence File
#   """ 
#   if ascii == None:
#         fi = open(self.ascii[:-4]+".tmp.out", "w")
#        else:
#            fi = open(ascii, "w")
#
## Write out header information
#
#        fi.write(" obs_sequence\n")
#        fi.write("obs_kind_definitions\n")
#
## Observation types are in the h5file.root.obs.kinds directory
#
#        kinds = h5file.root.obs.kinds
#
#        fi.write("       %d\n" % N.size(kinds))
#
#        for r in kinds.iterrows():
#            fi.write("    %d          %s   \n" % (r['index'], r['name']) )
#            if self.debug:  print(('pyDart/hdf2ascii:  Written observational types:  ', r))
#
#        attr = h5file.root.header.attributes
#        nobs = attr.col('num_obs')[0]
#
#        fi.write("  num_copies:            %d  num_qc:            %d\n" % (attr.col('num_copies')[0], 0 ))
#
#        if len(self.index) > 0:
#            fi.write(" num_obs:       %d  max_num_obs:       %d\n" % (len(self.index), len(self.index)) )
#
#            fi.write("observations\n")
#            if attr.col('num_copies')[0] == 2:
#                fi.write("truth\n")
#
#            fi.write("  first:            %d  last:       %d\n" % (1, len(self.index)) )
#            if self.debug:
#                print(("pyDart/hdf2ascii:  Max number of observations:    ", len(self.index)))
#
#        else:
#            fi.write(" num_obs:       %d  max_num_obs:       %d\n" % (attr.col('num_obs')[0], attr.col('max_num_obs')[0]))
#
#            fi.write("observations\n")
#            if attr.col('num_copies')[0] == 2:
#                fi.write("truth\n")
#
#            fi.write("  first:            %d  last:       %d\n" % (attr.col('first')[0], attr.col('last')[0]))
#            if self.debug:
#                print(("pyDart/hdf2ascii:  Max number of observations:    ", attr.col('max_num_obs')[0]))
#




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
