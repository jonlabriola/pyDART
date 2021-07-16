# pyDART

There are four primary goals supported by this package:

1.) Grab existing observations or simulated observations (all saved in the pyDART class) and save within an NetCDF File 
   -Initially do by observation type (MRMS, NEXRAD, Idealized Radar,  Conventional, Idealized Conventional)

   a.) idealRadObs.py - Generates simulated radar observations (vr/dBZ) at the specified locations

       generates radar observation class that contains full volume scans from different radars
   b.) idealConvObs.py - Generates conventional observations at specified locatiosn
   			 includes surface stations, profilers, soundings  
       generates conventional observation class that contains soundings, profilers, and surface obs

   Both programs save pickle files of the class instance, which are then used to produce the subsequent radar observations

2.) Convert the combined radar and conventional observations classes into an observation sequence file (hdf2ascii)
   -hdf2ascii program used to do this

3.) Convert the observation sequence file back and place into observation class for appropriate time 
   -ascii2hdf program used to do this
