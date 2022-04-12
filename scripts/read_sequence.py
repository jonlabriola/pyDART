import pydart,pickle,argparse

parser = argparse.ArgumentParser()
parser.add_argument("path",type=str,help='The full path to the obs_sequence file you want to read in')
parser.add_argument("--time",type=str,default='20090516002300',help='The time of the obs_sequence file (YYYYMMDDHHMMSS)')
#parser.add_argument("varname",type=str,help='The variables used to plot RADAR_REFLECTIVITY, DOPPLER_RADIAL_VELOCITY, fine_z')
parser.add_argument("--radar",action='store_true',help="When you want to read in the radar variables")
parser.add_argument("--conv",action='store_true',help="When you want to read in the profiler observations")
#parser.add_argument("platform",type=int,help='Work only with radar obs')
parser.add_argument("--pickle",action='store_true',help='When working with a pickle file')

arguments = parser.parse_args()
stored_pickle = False 

if arguments.radar:
   if arguments.pickle:
      radar = pickle.load(open('radar_%s.pickle'%arguments.time,"rb"))
   else:
      ntilt = 14
      #ny = 99
      #nx = 99
      #ny = 40
      #nx = 40
      (14, 59, 119)
      nx = 119 #50
      ny = 59 #119
      radar = pydart.sequence.read_sequence(arguments.path,ntilt=ntilt,ny=ny,nx=nx,radar=True)
      pickle.dump(radar,open('radar_%s.pickle'%arguments.time, "wb" ))

   #radar_name = 'radar_%03d'%arguments.platform
   radar_name = radar.obs.keys()[0] 
   #plotvar = radar.obs[radar_name][arguments.varname]
   varnames = ['RADAR_REFLECTIVITY','DOPPLER_RADIAL_VELOCITY']
   for varname in varnames:
      plotvar = radar.obs[radar_name][arguments.varname]
      for key in plotvar.keys():
         print('Working on key %s'%key)
         if key in ['xloc','yloc','zloc','flag','DART_QC']:
            pass
         else:
            outname = 'seq_'
            for key_split in key.split():
               outname += key_split+'_'
            pydart.plotting.twod_plot(plotvar,arguments.varname,tilt=1,outname='%s%s_%s.png'%(outname,arguments.varname,radar_name),copy_name=key)

if arguments.conv:
   if arguments.pickle:
      conv = pickle.load(open('conventional_%s.pickle'%arguments.time,"rb"))
   else:
      conv = pydart.sequence.read_sequence(arguments.path,conv=True)
      pickle.dump(conv,open('conventional_%s.pickle'%arguments.time, "wb" ))

for plat_name in conv.obs.keys():
   print('The platform is ...',plat_name)
   if 'pro' in plat_name:
      varnames = ['RADIOSONDE_U_WIND_COMPONENT','RADIOSONDE_V_WIND_COMPONENT',
                        'RADIOSONDE_TEMPERATURE','RADIOSONDE_DEWPOINT']
   elif 'sfc' in plat_name:
      varnames = ['TEMPERATURE_2M','DEWPOINT_2_METER','U_WIND_10M','V_WIND_10M']
   for varname in varnames:
      plotvar = conv.obs[plat_name][varname]
      if 'sfc' in plat_name: print('varname = ',varname)
   #   for key in plotvar.keys():
   #      print('Working on key %s'%key)
