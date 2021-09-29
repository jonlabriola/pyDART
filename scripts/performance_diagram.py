import pydart,argparse,pickle
import numpy as np
import matplotlib.pyplot as plt

#--- This Code Loops Through Each Ensemble
#--- Calculates the FOH/POD of each Ensemble Member
#--- Plots Individual Ensemble Members and Then the Ensemble Aeverage

parser = argparse.ArgumentParser()
parser.add_argument("var",type=str,help = 'Var to plot (RADAR_REFLECTIVITY, fine_z, etc.)')
parser.add_argument("--time",type=str,default='200905160033',help='The location of the dmax file')
parser.add_argument("--platform",type=int,default=1,help='The platform number (default = 1)')
parser.add_argument("--tilt",type=int,default=1,help='Radar tilt (default = 1)')
arguments = parser.parse_args()

#figure = pydart.plotting.gen_performance() 
platform_name = 'radar_%03d'%arguments.platform
cols = ['r','b','g','c','y','k','m','lime','cornflowerblue','tan']
for ens in range(1,11):
   print('Ensemble = ',ens)
   figure = pydart.plotting.gen_performance() #--- Generate the Initial Performance Diagram
   #obs_path = '/work/jonathan.labriola/OSSE/QLCS/mem%03d/radar_obs_%s.pickle'%(ens,arguments.time)
   obs_path = '/work/jonathan.labriola/OSSE/QLCS/WRONG_NAMING/mem%03d/radar_obs_%s.pickle'%(ens,arguments.time)
   radar = pickle.load(open(obs_path,"rb"))
   obs = radar.obs[platform_name][arguments.var]['obs'][arguments.tilt]
   FOH_ENS = []
   POD_ENS = []
   CSI_ENS = []
   BIAS_ENS = []

   for mem in range(1,41):
      fcst_path = '/scratch/jonathan.labriola/osse/ensemble/forecasts/RADAR/ens%03d/mem%03d/radar_obs_%s.pickle'%(ens,mem,arguments.time)
      radar = pickle.load(open(fcst_path,"rb"))
      fcst = radar.obs[platform_name][arguments.var]['obs'][arguments.tilt]

      #--- Calculate Statistics For Each Forecast Member
      FOH,POD,CSI,BIAS = pydart.verification.performance_scores(fcst,obs,45,kernel=3.)
      plt.scatter(FOH,POD,color=cols[ens-1],alpha=0.5,s=200)
      FOH_ENS.append(FOH) #--- Saving to save an ensemble average 
      POD_ENS.append(POD)
      CSI_ENS.append(CSI)
      BIAS_ENS.append(BIAS)

   plt.scatter(np.mean(FOH_ENS),np.mean(POD_ENS),color='k',s=500)
   #txt = plt.text(0.975,0.95,'MEAN CSI = %4.3f\nMEAN BIAS = %4.3f'%(np.mean(CSI),np.mean(BIAS_ENS)),fontsize=15,horizontalalignment='right', verticalalignment='center')
   txt = plt.text(0.025,0.95,'MEAN CSI = %4.3f\nMEAN BIA = %4.3f'%(np.mean(CSI),np.mean(BIAS_ENS)),fontsize=15,horizontalalignment='left', verticalalignment='center')
   txt.set_bbox(dict(facecolor='white', alpha=0.9, edgecolor='k'))
   plt.savefig('Performance_%s_ens%03d.png'%(arguments.time,ens),dpi=300)
   plt.clf()
