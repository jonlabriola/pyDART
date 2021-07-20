import sys, os
import string
import re
import glob
import pylab as P
import numpy as N
import time
import netCDF4 as ncdf
#from pyproj import Proj
from optparse import OptionParser
from tables import *
#from netcdftime import utime
#from datetime import datetime as py_datetime

#=========================================================================================
# PARAMETERS FOR MRMS

dbz_zeros         = False           # True: missing < dBZ < dbz_min=>0, False: set to missing (creates halo)
dbz_missing_zeros = True            # True => fill _missing with zeros, False => OFF
dbz_min           = 20.             # Reflectivity values below this are _missing or set to dbz_clear_air
dbz_max           = 70.             # Set reflectivity values above this value to dbz_max
dbz_thin          = 4               # thin data by skipping this in x&y (1 ==> no thinning, usually 3 or 6)
dbz_thin_z        = 2               # Vertically thin data by this factor (usually 2)
dbz_clear_thin    = 2               # Thin zeros additionally by skipping this in x&y (1 ==> no thinning)
dbz_clear_air     = 0.0             # value for clear air reflectivity (usually 0 or -10)
dbz_clear_hgts    = (True, 3000., 7000.)
dbz_stdev         = 7.5
lat_bound_box     = (33.0,38.0)     # max and min latitudes to output
lon_bound_box     = (-100.,-95.0)   # max and min longitudes to output
hgt_bound_box     = (0.0,10000.0)   # max and min heights to output
dbz_verts         = [(-2000.,-2000.),(-2000.,2000.),(2000.,2000.),(2000.,-2000.)]
dbz_missing_fudge = -99.

# read MRMS created superobbed sweep (on a rectangular grid) and create a
#      DART data base file.  Code also thins zero reflecitivty obs, thresholds
#      a max and min, and removes data close to a radar.

def mrms(self, ncdf_file, filename=None, lat_bbox=None, lon_bbox=None, hgt_bbox = hgt_bound_box):

        if filename != None:
            self.file(filename=filename)

# Open netcdf file(s)

        if ncdf_file.find(".netcdf") != -1:
            fncdf = [ncdf_file]
        else:
            fncdf = glob.glob(ncdf_file+"/*.netcdf")

        print(("\npyDart/MRMS->  Number of files to be read in:  %d" % (len(fncdf))))
        print(("\npyDart/MRMS->  First file:  %s " % fncdf[0]))
        if len(fncdf) > 1:  print(("\npyDart/MRMS->  Last file:  %s " % fncdf[-1]))

        if dbz_zeros:
            print(("\npyDart/MRMS:  Values of reflectivity < %d are set to clear_air values of: %d" % (dbz_min, dbz_clear_air)))
        else:
            print("\npyDart/MRMS:  dbz_zero flag is False")

        if dbz_missing_zeros:
            print(("\npyDart/MRMS:  Missing reflectivity are set to clear_air values of  %d" % dbz_clear_air))
        else:
            print("\npyDart/MRMS:  dbz_missing zero flag is False")

        if dbz_max:
            print(("\npyDart/MRMS: Values of reflectivity greater than %d will be set to %d" % (dbz_max, dbz_max)))

        print(("\npyDart/MRMS:  dBZ values will be thinned horizontally by a factor of %d" % dbz_thin))

        print(("\npyDart/MRMS:  dBZ values will be thinned vertically by a factor of %d" % dbz_thin_z))

        if dbz_clear_thin > 1:
             print(("\npyDart/MRMS:  Clear air reflectivity values thinned by an additional factor of %d" % dbz_clear_thin))
        else:
             print("\npyDart/MRMS:  No extra thning of clear air data requested")

        print(("\npyDart/MRMS->  Lat Bounding Box on:  %f to %f" % (lat_bbox[0], lat_bbox[1])))
        print(("\npyDart/MRMS->  Lon Bounding Box on:  %f to %f" % (lon_bbox[0], lon_bbox[1])))

# Create PyTables file

        filter_spec = Filters(complevel=5, complib="zlib", shuffle=1, fletcher32=0)

        h5file = open_file(self.hdf5, mode = "w", title = version_string, filters=filter_spec)

        group_ob_kinds = h5file.create_group("/", 'obs', 'Obs for DART file')

        group_header = h5file.create_group("/", 'header', 'Header Information for DART file')

# Find the obs group to create table in

        root         = h5file.root
        group_obs    = root.obs
        group_header = root.header

# create table that will hold the observation information

        table_obs = h5file.create_table(group_obs, 'observations', DART_obs, 'Observations from DART file')

        row = table_obs.row

# set counters

        count       = 0
        count_dbz   = 0
        count_zeros = 0

#  Open netcdf files, read data, flatten, and append

        for file in fncdf:

            print(("\npyDart/MRMS->  Processing file:  %s" % file))

            f = ncdf.Dataset(file, "r")

            lats  = f.variables['Lat'][::dbz_thin]
            lons  = f.variables['Lon'][::dbz_thin]

            try:
                hgts  = f.variables['Hgt'][0:-1:dbz_thin_z]
            except KeyError:
                try:
                    hgts  = f.variables['Ht'][0:-1:dbz_thin_z]
                except KeyError:
                    hgts = N.array((f.Height,))

  # For radar reflectivity - do some pre-processing based on flags at top of file
  # if dbz_zeros - for all MISSING values, set to zero - except we use the dbz_thin to thin zero reflectivity
  # if dbz_min   - then set all values of dbz < dbz_min to MISSING - this helps create a halo of missing values
  #                around the actual storms.
  # if dbz_max   - finally, as a last step, clip the reflectivity to be less than dbz_max

            dbz = None


            try:
                dbz = f.variables['MergedReflectivityQCComposite'][0,::dbz_thin,::dbz_thin]
                print("\n==>MRMS2HDF:  Found MergedReflectivityQC variable in netCDF file!\n")
                dbz = N.expand_dims(dbz, axis=0)
            except KeyError:
                print("\n==>MRMS2HDF:  Cannot find MergedReflectivityQCComposite variable in netCDF file!")

            try:
                dbz = f.variables['MergedReflectivityQC'][0,::dbz_thin,::dbz_thin]
                print("\n==>MRMS2HDF:  Found MergedReflectivityQC variable in netCDF file!\n")
                dbz = N.expand_dims(dbz, axis=0)
            except KeyError:
                print("\n==>MRMS2HDF:  Cannot find MergedReflectivityQC variable in netCDF file!")

            try:
                dbz  = f.variables['mrefl_mosaic'][0,::dbz_thin_z,::-dbz_thin,::dbz_thin]
                print("\n==>MRMS2HDF:  Found mrefl_mosaic variable in netCDF file!\n")
                print("\n==>MRMS2HDF:  Fixing the missing values!\n")
                dbz = N.where(dbz < dbz_missing_fudge+0.1, f.MissingData, dbz)
            except KeyError:
                print("\n==>MRMS2HDF:  Cannot find mrefl_mosaic variable in netCDF file!\n")

            if dbz == None:
                print("\n==>MRMS2HDF:  Cannot find any of the specified reflectivity variables in netCDF file, exiting....\n")
                raise SystemExit

# Figure out the bounding box indices to make the output go a LOT faster

            index  = (lats > lat_bbox[0] ) & (lats < lat_bbox[1])
            lindex = N.arange(lats.size)[index]
            index  = (lons > lon_bbox[0] ) & (lons < lon_bbox[1])
            mindex = N.arange(lons.size)[index]
            index  = (hgts >= hgt_bbox[0] ) & (hgts < hgt_bbox[1])
            kindex = N.arange(hgts.size)[index]

# Now process the data (note, 2D arrays still have a 'k'-index, the k-loop executes once)

            for k in kindex:

                if dbz_clear_hgts[0] == True:
                    if hgts[k] == dbz_clear_hgts[1] or hgts[k] == dbz_clear_hgts[2]:
                        hgt_flag = 0
                    else:
                        hgt_flag = 1
                else:
                    hgt_flag = 0

                hgt = hgts[k]

                for l in lindex:

                    lat = lats[l]

                    for m in mindex:

                        lon = lons[m]

                        thin_zeros = hgt_flag + (l + m) % dbz_clear_thin

                        data = min(dbz[k,l,m], dbz_max)

                        my_mask = False

                        if data == f.MissingData and dbz_missing_zeros and thin_zeros == 0:
                              data = dbz_clear_air
                              data_kind = ObType_LookUp("CLEARAIR_REFLECTIVITY")
                              my_mask = True
                              count_zeros += 1

                        if (data > f.MissingData and data < dbz_min) and dbz_zeros and thin_zeros == 0:
                              data = dbz_clear_air
                              data_kind = ObType_LookUp("CLEARAIR_REFLECTIVITY")
                              my_mask = True
                              count_zeros += 1

                        if data >= dbz_min:
                              my_mask = True
                              data_kind = ObType_LookUp("REFLECTIVITY")
                              count_dbz += 1

                        if my_mask:

                            row['number'] = int(count+1)
                            row['value']  = data

                            if  count == 0:
                                row['previous'] = -1
                            else:
                                row['previous'] = count-1

                            row['next']       = count+1
                            row['cov_group']  = -1

            # Set lat and lon in the codes....

                            row['lon']        = lon
                            row['lat']        = lat
                            row['height']     = hgt
                            row['z']          = hgt
                            row['vert_coord'] = 3
                            row['kind']       = data_kind

            # Add in a full date and time string to data set, as well as create a UTIME in seconds for searching

                            dt = ncdf.num2date(f.Time, units="seconds since 1970-01-01 00:00:00")
                            gc = ncdf.date2num(dt, units = "days since 1601-01-01 00:00:00")

                            row['days']    = N.int(gc)
                            row['seconds'] = (gc - N.int(gc)) * 86400

                            row['date']    = dt.strftime(time_format)
                            #row['utime']   = f.Time
                            row['error_var'] = dbz_stdev

                            row['index'] = count

                            count = count + 1

                            if count % 5000 == 0:
                                print(("PyDart.MRMS:  Processed observation # %d" %  (count+1)))
                                #print((" Date,sec_utime,utime = ", row['date'], row['utime']))

                            row.append()
            f.close()

        table_obs.flush()

        print(("\n >%s<" % ("=" * 79)))
        print(("\n  PyDart.MRMS->  Total number of observations processed observation %d" % int(count)))
        print(("\n  PyDart.MRMS->  Total number of dBZ > dbz_min: %d" % int(count_dbz)))
        print(("\n  PyDart.MRMS->  Total number of ZERO dBZs:     %d" % int(count_zeros)))
        print(("\n >%s<" % ("=" * 79)))

# Create header information for the table......

        table_ob_kinds = h5file.create_table(group_ob_kinds, 'kinds', DART_ob_kinds, 'Observation Descriptions')
        row = table_ob_kinds.row

        if count_dbz > 0:
            row['index'] = ObType_LookUp("DBZ")
            row['name']  = "RADAR_REFLECTIVITY" #"DBZ"
            row.append()

        if count_zeros > 0:
            row['index'] = ObType_LookUp("CLEARAIR_REFLECTIVITY")
            row['name']  = "CLEARAIR_REFLECTIVITY" #"DBZ"
            row.append()

        table_ob_kinds.flush()

        if self.debug:  print('Completed reading obs_kind_definitions')

# Create file header information

        table_header = h5file.create_table(group_header, 'attributes', DART_header, 'Attributes of the observational file')

        row = table_header.row

        row['num_copies']  = 1
        row['num_qc']      = 0
        row['num_obs']     = count
        row['max_num_obs'] = count
        row['first']       = 1
        row['last']        = count

        row.append()

        table_header.flush()

        if self.debug:
                print(("Number of observation copies:  ", 1))
                print(("Number of QC'd observations:   ", 0))
                print(("Number of observations:        ", count))
                print(("Max number of observations:    ", count))

        h5file.close()

        print("\npyDART.MRMS->  Converted MRMS file to HDF5 DART file")

        return



#-------------------------------------------------------------------------------
# Main function defined to return correct sys.exit() calls

def main(argv=None):
    if argv is None:
           argv = sys.argv

# Init some local variables

    start = None
    end   = None

# Initialize class object...

    #myDART = pyDART(verbose=_verbose, debug=_debug)

# Command line interface for PyDart

    parser = OptionParser()
    parser.add_option("-f", "--file",        dest="file",      type="string", help = "Filename of ascii or PyDART file for conversion/search")
    #parser.add_option(      "--lat_box",     dest="lat_box",   default=None,  type = "float",  nargs=2, help = "Search for MRMS within these lat limits. Usage:  --lat_box lat_south lat_north")
    #parser.add_option(      "--lon_box",     dest="lon_box",   default=None,  type = "float",  nargs=2, help = "Search for MRMS within these lon limits. Usage:  --lon_box lon_west lon_east")

    (options, args) = parser.parse_args()

    #if options.lat_box:
    #   lat_bbox = options.lat_box
    #else:
    lat_bbox = lat_bound_box

    #if options.lon_box:
    #   lon_bbox = options.lon_box
    #else:
    lon_bbox = lon_bound_box

    mrms(options.mrms, filename = options.file, lat_bbox = lat_bbox, lon_bbox = lon_bbox )


#-------------------------------------------------------------------------------
# Main program for testing...
#
if __name__ == "__main__":
    sys.exit(main())

# End of file
