# pyDART

There are four primary goals supported by this package:

1.) Grab existing observations or simulated observations (all saved in the pyDART class) and save within an NetCDF File 
   -Initially do by observation type (MRMS, NEXRAD, Idealized Radar,  Conventional, Idealized Conventional)

   a.) idealRadObs.py - Generates simulated radar observations (vr/dBZ) at the specified locations

       generates radar observation file that contains full volume scans from different radars
   b.) idealConvObs.py - Generates conventional observations at specified locatiosn
   			 includes surface stations, profilers, soundings  
       generates conventional observation file that contains soundings, profilers, and surface obs

2.) Combine all unique observations into a single netcdf file
   

2.) Convert the combined NetCDF file into an observation sequence file (hdf2ascii)
   -hdf2ascii program used to do this

3.) Convert the observation sequence file back into 
   -ascii2hdf program used to do this
