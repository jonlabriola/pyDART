#import pickle, numpy, pydart, os
#import numpy as np
import pydart,pickle



#filepath = '/scratch/jonathan.labriola/osse/cm1_test/cm1_filter/archive/20090515232300/obs_seq.final' 
#--- Create a Sequence Files
#filepath = '../output/20090515232300_obs_seq.prior'
#filepath = '../output/20090515232300_obs_seq.prior'
#filepath = '/scratch/jonathan.labriola/osse/obs/20090516001300_obs_seq.prior'
#filepath = '../output/20090516001300_obs_seq.prior'
#filepath = '/scratch/jonathan.labriola/osse/cm1_test/cm1_filter/archive/20090516001300/obs_seq.final'
filepath = '/scratch/jonathan.labriola/osse/cm1_test/cm1_filter/archive/20090515232300/obs_seq.final'
stored_pickle = False
#varname = 'RADAR_REFLECTIVITY'
varname = 'DOPPLER_RADIAL_VELOCITY'
#outname = 'obs_sequence_pic.png'

if stored_pickle:
  radar = pickle.load(open('TEST.pickle',"rb"))
else:
   ntilt = 14
   ny = 100
   nx = 100
   radar = pydart.sequence.read_sequence(filepath,ntilt=ntilt,ny=ny,nx=nx)
   pickle.dump(radar,open('TEST.pickle', "wb" ))

plotvar = radar.obs['radar_001'][varname]

for key in plotvar.keys():
   print('Working on key %s'%key)
   if key in ['xloc','yloc','zloc','flag','DART_QC']:
      pass
   else:
      outname = 'seq_'
      for key_split in key.split():
         outname += key_split+'_'
      pydart.plotting.twod_plot(plotvar,varname,tilt=1,outname='%s%s.png'%(outname,varname),copy_name=key)
