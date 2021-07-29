import pickle, numpy, pydart, os
import numpy as np

#--- Convert either observation class to obs_sequence file 

#--- or obs_sequence file to observation class

#--- flatten out radar observations into obs_sequence (allows the observations to be reintegrated)

time = 10
create_sequence = True
create_class = False 
types = ['radar','conventional']
time = '200905152313'
#radar_class = '../output/radar_obs_200905160138.pickle'
#conv_class = '../output/conv_obs_200905160143.pickle' 
#conv_class = '../output/conv_obs_200905160208.pickle'
conv_class = '../output/conv_obs_%s.pickle'%time
#--- Merge Together Radar/Conventional Observations
conv = pickle.load(open(conv_class,"rb"))
#radar = pickle.load(open(radar_class,"rb"))
#observations = {**radar.obs, **conv.obs} #--- Merge Observation Dictionaries
observations = conv.obs

#--- Create a Sequence Files
filepath = '../output/ob_seq_%s.out'%time
obsseq = pydart.sequence.create_sequence(filepath)
obsseq.add_header(observations)   #--- Make Sequence Header
obsseq.add_entries(observations)  #--- Make Observation Entries
