import pydart,argparse,pickle
import numpy as np
import matplotlib.pyplot as plt

#--- This Code Loops Through Each Ensemble
#--- Calculates the FOH/POD of each Ensemble Member
#--- Plots Individual Ensemble Members and Then the Ensemble Aeverage

parser = argparse.ArgumentParser()
parser.add_argument("var",type=str,help = 'Var to plot (RADAR_REFLECTIVITY, fine_z, etc.)')
parser.add_argument("--kernel",type=float,default=2.5,help='The kernel radius applied to forecasts/observations')
parser.add_argument("--time",type=str,default='20090516005300',help='The location of the dmax file')
parser.add_argument("--platform",type=int,default=1,help='The platform number (default = 1)')
parser.add_argument("--tilt",type=int,default=1,help='Radar tilt (default = 1)')
parser.add_argument("--nen",type=int,default=10,help='The number of ensembles (default=10)')
parser.add_argument("--nmem",type=int,default=40,help='The number of members in an ensemble (default=40)')
parser.add_argument("--no_prof",action='store_true',help='Select experiment with no profiler during DA')
parser.add_argument("--mor",action='store_true',help='Select experiments run with Morrison Scheme')
arguments = parser.parse_args()


if arguments.no_prof:
   da_tag = 'ens_rdr'
else:
   da_tag = 'ens_all'
if arguments.mor:
   mp_tag = 'Morrison'
else:
   mp_tag = 'NSSL'


platform_name = 'radar_%03d'%arguments.platform
cols = ['r','b','g','c','y','orangered','m','lime','cornflowerblue','tan']
FOH_ENS = np.zeros((arguments.nen))
POD_ENS = np.zeros((arguments.nen))
CSI_ENS = np.zeros((arguments.nen))
BIA_ENS = np.zeros((arguments.nen))

#--- Loop through each ensemble, plot individual members and ensemble mean
for eindex,ens in enumerate(range(1,arguments.nen+1)):
   print('Ensemble = ',ens)
   figure = pydart.plotting.gen_performance() #--- Generate the Initial Performance Diagram
   obs_path = '/work/jonathan.labriola/OSSE/QLCS/ens%03d/radar_obs_%s.pickle'%(ens,arguments.time)
   radar = pickle.load(open(obs_path,"rb"))
   obs = radar.obs[platform_name][arguments.var]['obs'][arguments.tilt]
   FOH = np.zeros((arguments.nmem))
   POD = np.zeros((arguments.nmem))
   CSI = np.zeros((arguments.nmem))
   BIA = np.zeros((arguments.nmem))

   #--- Loop through each member, grab statistics and Plot
   for mindex,mem in enumerate(range(1,41)):
      #if arguments.no_prof:
      #   fcst_path = '/scratch/jonathan.labriola/osse/ens_rdr/ens%03d/mem%03d/radar_obs_%s.pickle'%(ens,mem,arguments.time)
      #else:
      fcst_path = '/scratch/jonathan.labriola/osse/%s/%s/ens%03d/mem%03d/radar_obs_%s.pickle'%(mp_tag,da_tag,ens,mem,arguments.time)
      print('Opening ...',fcst_path)
  
      radar = pickle.load(open(fcst_path,"rb"))
      fcst = radar.obs[platform_name][arguments.var]['obs'][arguments.tilt]
      print('Max forecast Z = ',np.nanmax(fcst))
      FOH[mindex],POD[mindex],CSI[mindex],BIA[mindex] = pydart.verification.performance_scores(fcst,obs,45,kernel=arguments.kernel)
      plt.scatter(FOH[mindex],POD[mindex],color=cols[ens-1],alpha=0.5,s=200)

   plt.scatter(np.mean(FOH),np.mean(POD),color='k',s=500)
   txt = plt.text(0.025,0.95,'MEAN CSI = %4.3f\nMEAN BIA = %4.3f'%(np.mean(CSI),np.mean(BIA)),fontsize=15,horizontalalignment='left', verticalalignment='center')
   txt.set_bbox(dict(facecolor='white', alpha=0.9, edgecolor='k'))
   plt.savefig('Performance_%s_ens%03d_%3.2f_%s_%s.png'%(arguments.time,ens,arguments.kernel,mp_tag,da_tag),dpi=300)
   plt.clf()

   #--- Save the Ensemble Mean Statistics For Final Plot 
   FOH_ENS[eindex] = np.mean(FOH)
   POD_ENS[eindex] = np.mean(POD)
   CSI_ENS[eindex] = np.mean(CSI)
   BIA_ENS[eindex] = np.mean(BIA)

figure = pydart.plotting.gen_performance()
for eindex,ens in enumerate(range(1,arguments.nen+1)): 
   plt.scatter(FOH_ENS[eindex],POD_ENS[eindex],color=cols[eindex],s=500,label='ENS%03d: CSI(%2.2f) BIAS(%2.2f)'%(ens,CSI_ENS[eindex],BIA_ENS[eindex]))
plt.scatter(np.mean(FOH_ENS),np.mean(POD_ENS),color='k',s=500,label='MEAN: CSI(%2.2f) BIAS(%2.2f)'%(np.mean(CSI_ENS),np.mean(BIA_ENS)))
plt.legend(loc=2)
plt.savefig('Performance_%s_mean_%3.2f_%s_%s.png'%(arguments.time,arguments.kernel,mp_tag,da_tag),dpi=300)
plt.clf()
