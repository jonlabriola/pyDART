import pickle, numpy, pydart, os, argparse
import numpy as np
#--- Convert either observation class to obs_sequence file 
parser = argparse.ArgumentParser()
parser.add_argument("file_loc",type=str,default='./',help='The pathway the pickle files are located')
parser.add_argument("--time",type=str,default='20090515232300',help='Evaluated Output Time (YYYYMMDDHHMMSS)')
parser.add_argument("--conv",action='store_true',help='Work only with conventional obs')
parser.add_argument("--radar",action='store_true',help='Work only with radar obs')
parser.add_argument("--outpath",type=str,default='N/A',help='Output path, if desired')
arguments = parser.parse_args()

#conv_class = '/work/jonathan.labriola/python_scripts/pyDART/output/conv_obs_%s.pickle'%arguments.time
#radar_class = '/work/jonathan.labriola/python_scripts/pyDART/output/radar_obs_%s.pickle'%arguments.time
conv_class = '%s/conv_obs_%s.pickle'%(arguments.file_loc,arguments.time)
radar_class = '%s/radar_obs_%s.pickle'%(arguments.file_loc,arguments.time)

if arguments.conv:
   conv = pickle.load(open(conv_class,"rb"))
   observations = conv.obs  

elif arguments.radar:
   radar = pickle.load(open(radar_class,"rb"))
   observations = radar.obs  

else:
   conv = pickle.load(open(conv_class,"rb"))
   radar = pickle.load(open(radar_class,"rb"))
   observations = {**radar.obs, **conv.obs} #--- Merge Observation Dictionaries
   
#--- Create a Sequence Files
if 'N/A' in arguments.outpath:
   filepath = '/work/jonathan.labriola/python_scripts/pyDART/output/%s00_obs_seq.prior'%arguments.time
else:
   filepath = '%s/%s00_obs_seq.prior'%(arguments.outpath,arguments.time)
   
pydart.sequence.create_sequence(filepath,observations)
