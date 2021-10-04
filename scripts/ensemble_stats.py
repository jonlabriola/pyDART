import pydart,argparse,pickle
import numpy as np
import matplotlib.pyplot as plt

#--- This Code Loops Through Each Ensemble
#--- Calculates the FOH/POD of each Ensemble Member
#--- Plots Individual Ensemble Members and Then the Ensemble Aeverage

parser = argparse.ArgumentParser()
parser.add_argument("var",type=str,help = 'Var to plot (RADAR_REFLECTIVITY, fine_z, etc.)')
parser.add_argument("--mor",action='store_true',help='Load experiments run with Morrison Scheme')
parser.add_argument("--no_prof",action='store_true',help='Load experiment run with no profilers during DA')
parser.add_argument("--kernel",type=float,default=2.5,help='The kernel radius applied to forecasts/observations')
parser.add_argument("--time",type=str,default='20090516003300',help='The location of the dmax file')
parser.add_argument("--platform",type=int,default=1,help='The platform number (default = 1)')
parser.add_argument("--threshold",type=float,default=45.,help='The verified threshold (default=45)')
parser.add_argument("--tilt",type=int,default=1,help='Radar tilt (default = 1)')
parser.add_argument("--nen",type=int,default=10,help='The number of ensembles (default=10)')
parser.add_argument("--nmem",type=int,default=40,help='The number of members in an ensemble (default=40)')
parser.add_argument("--nmep",action='store_true',help='Plot the reliability')
parser.add_argument("--reliability",action='store_true',help='Plot the reliability')
parser.add_argument("--pmmean",action='store_true',help='Plot the reliability')
parser.add_argument("--mean",action='store_true',help='Plot the reliability')
parser.add_argument("--auc",action='store_true',help='Print AUC value')
arguments = parser.parse_args()

platform_name = 'radar_%03d'%arguments.platform
#--- Selecting all the specified tags
if arguments.no_prof:
   exp_tag = 'ens_rdr'
else:
   exp_tag = 'ens_all'
if arguments.mor:
   mp_tag = 'Morrison'
else:
   mp_tag = 'NSSL'

#--- Loop through each ensemble, plot individual members and ensemble mean
for eindex,ens in enumerate(range(1,arguments.nen+1)):
   print('Ensemble = ',ens)
   obs_path = '/work/jonathan.labriola/OSSE/QLCS/ens%03d/radar_obs_%s.pickle'%(ens,arguments.time)
   radar = pickle.load(open(obs_path,"rb"))
   obs = radar.obs[platform_name][arguments.var]['obs'][arguments.tilt]
   for mindex,mem in enumerate(range(1,arguments.nmem+1)):
      if arguments.no_prof:
         fcst_path = '/scratch/jonathan.labriola/osse/%s/%s/ens%03d/mem%03d/radar_obs_%s.pickle'%(mp_tag,exp_tag,ens,mem,arguments.time)
      else:
         fcst_path = '/scratch/jonathan.labriola/osse/%s/%s/ens%03d/mem%03d/radar_obs_%s.pickle'%(mp_tag,exp_tag,ens,mem,arguments.time)
      radar = pickle.load(open(fcst_path,"rb"))
      if mindex == 0:
         tmp = radar.obs[platform_name][arguments.var]['obs'][arguments.tilt]
         [nyy,nxx] = tmp.shape
         fcst = np.zeros((arguments.nmem,nyy,nxx))
         fcst[mindex] = tmp
      else:
         fcst[mindex] = radar.obs[platform_name][arguments.var]['obs'][arguments.tilt]

   if arguments.pmmean:  #---Probabilist matched mean
      var = pydart.verification.pmmean(fcst)
      varname = arguments.var

   elif arguments.mean:  #--- Probability matched mean
      var = np.mean(fcst,axis=0)
      varname = arguments.var

   else:  #--- NMEP/Reliability/AUC
      var = pydart.verification.nmep(fcst,arguments.kernel,arguments.threshold)      
      obs = pydart.verification.gen_buffer(obs,arguments.kernel,arguments.threshold,ones=False) 
      varname = 'NEP'

      if arguments.reliability:
         sample_climo,no_skill,obs_frequency,bin_centers,bin_climo,reliability_bins,bin_count = pydart.verification.createReliability(var,obs,arguments.threshold)
         figure = plt.figure
         BSS = pydart.verification.BSS(var,obs,arguments.threshold)
         pydart.plotting.plotReliability(sample_climo,no_skill,obs_frequency,bin_centers,bin_climo,linecolor = 'r',label='BSS = (%0.3f)'%(BSS))
         plt.legend()
         plt.savefig('Reliability_%s_ens%03d_%3.2f_%3.2f_%s.png'%(arguments.time,ens,arguments.kernel,arguments.threshold,exp_tag),dpi=300)
         plt.clf()
      elif arguments.auc:
         var_auc,pod,pofd = pydart.verification.auc(var,obs,arguments.threshold) 
         print('The AUC for ens%03d is .. %3.3f'%(ens,var_auc))

   if arguments.reliability or arguments.auc:
      pass
   else:
      outname ='%s_%s_ens%03d_%3.2f_%3.2f_%s_%s.png'%(varname,arguments.time,ens,arguments.kernel,arguments.threshold,exp_tag,mp_tag) 
      if arguments.pmmean or arguments.mean:
         pydart.plotting.rough_plot(var,varname,outname=outname)
      else:
         pydart.plotting.rough_plot(var,varname,outname=outname,contour_var=obs,threshold=arguments.threshold)

      plt.clf()
