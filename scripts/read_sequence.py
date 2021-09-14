import pydart,pickle,argparse

parser = argparse.ArgumentParser()
parser.add_argument("path",type=str,default='200905152323',help='The last output time (if different from first)')
parser.add_argument("varname",type=str,help='The variables used to plot RADAR_REFLECTIVITY, DOPPLER_RADIAL_VELOCITY, fine_z')
parser.add_argument("platform",type=int,help='Work only with radar obs')
arguments = parser.parse_args()

#filepath = '/scratch/jonathan.labriola/osse/both/cm1_filter/archive/20090516002300/obs_seq.final'
stored_pickle = True

if stored_pickle:
  radar = pickle.load(open('TEST.pickle',"rb"))
else:
   ntilt = 14
   ny = 100
   nx = 100
   radar = pydart.sequence.read_sequence(arguments.path,ntilt=ntilt,ny=ny,nx=nx)
   pickle.dump(radar,open('TEST.pickle', "wb" ))

radar_name = 'radar_%03d'%arguments.platform
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
