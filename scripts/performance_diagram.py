import pydart,argparse,pickle
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
matplotlib.rcParams['axes.linewidth'] = 1.5
matplotlib.rcParams.update({'font.size': 22})
#--- This Code Loops Through Each Ensemble
#--- Calculates the FOH/POD of each Ensemble Member
#--- Plots Individual Ensemble Members and Then the Ensemble Aeverage

parser = argparse.ArgumentParser()
parser.add_argument("var",type=str,help = 'Var to plot (RADAR_REFLECTIVITY, fine_z, etc.)')
parser.add_argument("--kernel",type=float,default=2.5,help='The kernel radius applied to forecasts/observations')
parser.add_argument("--time",type=str,default='20090516005300',help='The location of the dmax file')
parser.add_argument("--platform",type=int,default=1,help='The platform number (default = 1)')
parser.add_argument("--tilt",type=int,default=1,help='Radar tilt (default = 1)')
parser.add_argument("--threshold",type=float,default=40,help='Threshold')
parser.add_argument("--nen",type=int,default=10,help='The number of ensembles (default=10)')
parser.add_argument("--nmem",type=int,default=40,help='The number of members in an ensemble (default=40)')
parser.add_argument("--no_prof",action='store_true',help='Select experiment with no profiler during DA')
parser.add_argument("--new_err",action='store_true',help='New Observation Errors')
parser.add_argument("--multi_prof",action='store_true',help='Select experiments run with Morrison Scheme')
parser.add_argument("--coarse_wind_shift",action='store_true',help='3 km experiments with a 7m/s wind shift in each direction')
parser.add_argument("--coarse",action='store_true',help='3 km experiments with a 7m/s wind shift in each direction')
parser.add_argument("--mor",action='store_true',help='Select experiments run with Morrison Scheme')
parser.add_argument("--coarse_perts",action='store_true',help='Experiments run on 3 km grid with initial perts')
parser.add_argument("--coarse_dryperts",action='store_true',help='Experiments run on 3 km grid with initial perts')
arguments = parser.parse_args()


if arguments.no_prof:
   da_tag = 'ens_rdr'
elif arguments.new_err:
   da_tag = 'ens_all_new_error'
else:
   da_tag = 'ens_all'
if arguments.mor:
   mp_tag = 'Morrison'
elif arguments.multi_prof:
   mp_tag = 'Multi_Sounding'
elif arguments.coarse_wind_shift:
   mp_tag = '3km_NSSL_u7v7'
elif arguments.coarse:
   mp_tag = '3km_NSSL'
elif arguments.coarse_perts:
   mp_tag = '3km_NSSL_pert'
elif arguments.coarse_dryperts:
   mp_tag = '3km_NSSL_dry'
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
   if arguments.coarse_wind_shift or arguments.coarse or arguments.coarse_perts or arguments.coarse_dryperts:
      obs_path = '/work/jonathan.labriola/OSSE/obs/3km/ens%03d/radar_obs_%s.pickle'%(ens,arguments.time)
      radar = pickle.load(open(obs_path,"rb"))
      obs = radar.obs[platform_name][arguments.var]['obs'][:,:-1,:-1]
   else:
      obs_path = '/work/jonathan.labriola/OSSE/QLCS/ens%03d/radar_obs_%s.pickle'%(ens,arguments.time)
      radar = pickle.load(open(obs_path,"rb"))
      obs = radar.obs[platform_name][arguments.var]['obs'][:,:,:]
   if arguments.tilt < 0:
      obs = np.nanmax(obs,axis=0)
   else:
      obs = obs[arguments.tilt]

  
   
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
      if arguments.tilt < 0:
         fcst = np.nanmax(radar.obs[platform_name][arguments.var]['obs'][:,:,:],axis=0)
      else:
         fcst = radar.obs[platform_name][arguments.var]['obs'][arguments.tilt]
      print('Max forecast Z = ',np.nanmax(fcst))
      FOH[mindex],POD[mindex],CSI[mindex],BIA[mindex] = pydart.verification.performance_scores(fcst,obs,arguments.threshold,kernel=arguments.kernel)
      plt.scatter(FOH[mindex],POD[mindex],color=cols[ens-1],alpha=0.5,s=200)

   plt.scatter(np.median(FOH),np.median(POD),color='k',s=500)
   txt = plt.text(0.025,0.95,'MEAN CSI = %4.3f\nMEAN BIA = %4.3f'%(np.median(CSI),np.median(BIA)),fontsize=15,horizontalalignment='left', verticalalignment='center')
   txt.set_bbox(dict(facecolor='white', alpha=0.9, edgecolor='k'))
   plt.savefig('Performance_%s_ens%03d_%3.2f_%s_%s.png'%(arguments.time,ens,arguments.kernel,mp_tag,da_tag),dpi=300)
   plt.clf()

   #--- Save the Ensemble Mean Statistics For Final Plot 
   FOH_ENS[eindex] = np.median(FOH)
   POD_ENS[eindex] = np.median(POD)
   CSI_ENS[eindex] = np.median(CSI)
   BIA_ENS[eindex] = np.median(BIA)

figure = pydart.plotting.gen_performance()
for eindex,ens in enumerate(range(1,arguments.nen+1)): 
   plt.scatter(FOH_ENS[eindex],POD_ENS[eindex],color=cols[eindex],s=500,label='ENS%03d: CSI(%2.2f) BIAS(%2.2f)'%(ens,CSI_ENS[eindex],BIA_ENS[eindex]))
plt.scatter(np.median(FOH_ENS),np.median(POD_ENS),color='k',s=500,label='MEDIAN: CSI(%2.2f) BIAS(%2.2f)'%(np.median(CSI_ENS),np.median(BIA_ENS)))
plt.xlabel('Success Ratio (1-FAR)')
plt.ylabel('Probability of Detection')
plt.legend(loc=2,prop={"size":10})
plt.savefig('Performance_%s_mean_%3.2f_%s_%s_%3.2fthresh.png'%(arguments.time,arguments.kernel,mp_tag,da_tag,arguments.threshold),dpi=300)
plt.clf()
