import pickle, numpy, pydart, os
import numpy as np

#--- Convert either observation class to obs_sequence file 

#--- or obs_sequence file to observation class

#--- flatten out radar observations into obs_sequence (allows the observations to be reintegrated)

time = 10
create_sequence = True
create_class = False 
types = ['radar','conventional']

radar_class = 'output/radar_obs_200905160138.pickle'
conv_class = './output/conv_obs_200905152228.pickle'

obname = []
obcode = []
obs_codes = pydart.readvar.obcode() 
dart_obs = obs_codes.keys() #--- A List of the different observation names
#--- Observation Sequence File Should Not Exist

#--- JDL Conducting a Test
conv = pickle.load(open(conv_class,"rb"))
radar = pickle.load(open(radar_class,"rb"))
observations = {**radar.obs, **conv.obs} #--- Merge Observation Dictionaries


nobs = 0 
for platform in observations.keys():          #--- Loop over platforms
   for key in observations[platform]:         #--- Loop over saved values
      if key in dart_obs:
         nobs += np.size(observations.obs[platform][key]['obs'])
         if key in obname:
            pass
         else:
            obname.append(key)
            obcode.append(obs_codes[key])

print('Total Assimilated Obser = ',nobs)
print('Observation name = ',obname)


#--- Sequence Files
try:
  filepath = "./output/obsequence.txt"
  os.system('rm %s'%filepath)
except:
  print('File does not exist')
file1 = open(filepath,"w")
pydart.sequence.add_header(file1,obname,obcode,nobs)
file1.close()
