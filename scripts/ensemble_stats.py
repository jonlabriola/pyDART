import pydart,argparse,pickle
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['axes.linewidth'] = 1.5
matplotlib.rcParams.update({'font.size': 22})
#--- This Code Loops Through Each Ensemble
#--- Calculates the FOH/POD of each Ensemble Member
#--- Plots Individual Ensemble Members and Then the Ensemble Aeverage

parser = argparse.ArgumentParser()
parser.add_argument("var",type=str,help = 'Var to plot (RADAR_REFLECTIVITY, fine_z, etc.)')
parser.add_argument("--mor",action='store_true',help='Load experiments run with Morrison Scheme')
parser.add_argument("--multi_sound",action='store_true',help='Load experiment run with no profilers during DA')
parser.add_argument("--no_prof",action='store_true',help='Load experiment run with no profilers during DA')
parser.add_argument("--limit",action='store_true',help='Load experiment run with no profilers during DA')
parser.add_argument("--coarse_wind_shift",action='store_true',help='Load experiment run on 3 km grid with u7/v7 wind shift')
parser.add_argument("--coarse",action='store_true',help='Load experiment run on 3 km grid')
parser.add_argument("--coarse_pert",action='store_true',help='Load experiment run on 3 km grid with random perturbations')
parser.add_argument("--coarse_drypert",action='store_true',help='Load experiment run on 3 km grid with dry random perturbations')
parser.add_argument("--kernel",type=float,default=2.5,help='The kernel radius applied to forecasts/observations')
parser.add_argument("--time",type=str,default='20090516003300',help='The location of the dmax file')
parser.add_argument("--platform",type=int,default=1,help='The platform number (default = 1)')
parser.add_argument("--threshold",type=float,default=45.,help='The verified threshold (default=45)')
parser.add_argument("--tilt",type=int,default=1,help='Radar tilt (default = 1)')
parser.add_argument("--nen",type=int,default=10,help='The number of ensembles (default=10)')
parser.add_argument("--nmem",type=int,default=40,help='The number of members in an ensemble (default=40)')
parser.add_argument("--nplat",type=int,default=1,help='Include number of platforms')
parser.add_argument("--nmep",action='store_true',help='Plot the reliability')
parser.add_argument("--reliability",action='store_true',help='Plot the reliability')
parser.add_argument("--pmmean",action='store_true',help='Plot the reliability')
parser.add_argument("--mean",action='store_true',help='Plot the reliability')
parser.add_argument("--auc",action='store_true',help='Print AUC value')
parser.add_argument("--save",action='store_true',help='Save a file with the ensemble information')
parser.add_argument("--xmin",default=0,type=int,help='X minimum')
parser.add_argument("--xmax",default=-1,type=int,help='X minimum')
parser.add_argument("--ymin",default=0,type=int,help='X minimum')
parser.add_argument("--ymax",default=-1,type=int,help='X minimum')
arguments = parser.parse_args()

platform_name = 'radar_%03d'%arguments.platform
#--- Selecting all the specified tags
if arguments.no_prof:
   exp_tag = 'ens_rdr'
else:
   exp_tag = 'ens_all'
if arguments.mor:
   mp_tag = 'Morrison'
elif arguments.multi_sound:
   mp_tag = 'Multi_Sounding'
elif arguments.coarse_wind_shift:
   mp_tag = '3km_NSSL_u7v7'
elif arguments.coarse:
   mp_tag = '3km_NSSL'
elif arguments.coarse_pert:
   mp_tag = '3km_NSSL_pert'
elif arguments.coarse_drypert:
   mp_tag = '3km_NSSL_dry'
else:
   mp_tag = 'NSSL'

#--- Loop through each ensemble, plot individual members and ensemble mean
for eindex,ens in enumerate(range(1,arguments.nen+1)):
   print('Ensemble = ',ens)
   if arguments.coarse_wind_shift or arguments.coarse or arguments.coarse_pert or arguments.coarse_drypert:
      #obs_path = '/work/jonathan.labriola/OSSE/obs/3km/ens%03d/radar_obs_%s.pickle'%(ens,arguments.time)
      obs_path = '/work/jonathan.labriola/OSSE/QLCS/Nature_Runs/simobs/ens%03d/radar_obs_%s.pickle'%(ens,arguments.time)
      radar = pickle.load(open(obs_path,"rb"))
      for platforms in range(1,arguments.nplat+1):
         platform_name = 'radar_%03d'%platforms
         print('Reading in ...',platform_name)
         obs_tmp = radar.obs[platform_name][arguments.var]['obs'][:,arguments.ymin:arguments.ymax,arguments.xmin:arguments.xmax]
         obs_tmp = np.where(np.isnan(obs_tmp),0,obs_tmp)
         if arguments.tilt < 0:
            obs_tmp = np.nanmax(obs_tmp,axis=0)
         else:
            obs_tmp = obs_tmp[arguments.tilt]
         if platforms == 1:
           obs = np.zeros(obs_tmp.shape)
         obs = np.where(obs>obs_tmp,obs,obs_tmp)

   else:
      print('Find a new platform')      
      #obs_path = '/work/jonathan.labriola/OSSE/QLCS/ens%03d/radar_obs_%s.pickle'%(ens,arguments.time)
      #radar = pickle.load(open(obs_path,"rb"))
      #obs = np.amax(radar.obs[platform_name][arguments.var]['obs'][:,:,:],axis=0)


   for mindex,mem in enumerate(range(1,arguments.nmem+1)):
      if arguments.no_prof:
         #fcst_path = '/scratch/jonathan.labriola/osse/%s/%s/ens%03d/mem%03d/radar_obs_%s.pickle'%(mp_tag,exp_tag,ens,mem,arguments.time)
         if arguments.limit:
            fcst_path = '/scratch/jonathan.labriola/osse/500m_NR/NSSL_fast4ms/ens_rdr/ens%03d_limit_poorcov/mem%03d/radar_obs_%s.pickle'%(ens,mem,arguments.time)
         else:
            fcst_path = '/scratch/jonathan.labriola/osse/500m_NR/NSSL_fast4ms/ens_rdr/ens%03d/mem%03d/radar_obs_%s.pickle'%(ens,mem,arguments.time)
      else:
         if arguments.limit:
            fcst_path = '/scratch/jonathan.labriola/osse/500m_NR/NSSL_fast4ms/ens_all/ens%03d_limit_poorcov/mem%03d/radar_obs_%s.pickle'%(ens,mem,arguments.time)
         else:
            fcst_path = '/scratch/jonathan.labriola/osse/500m_NR/NSSL_fast4ms/ens_all/ens%03d/mem%03d/radar_obs_%s.pickle'%(ens,mem,arguments.time)
         #fcst_path = '/scratch/jonathan.labriola/osse/%s/%s/ens%03d/mem%03d/radar_obs_%s.pickle'%(mp_tag,exp_tag,ens,mem,arguments.time)

      print('Opening ...',fcst_path)
      radar = pickle.load(open(fcst_path,"rb"))

      #--- Get Grid Information (Observation and Forecasts Should be on the same grid
      #--- radar_001 will always exist so just use this site
      xh,yh = pydart.plotting.calc_grid_dims(radar.obs['radar_001'][arguments.var]['xloc'][0],radar.obs['radar_001'][arguments.var]['yloc'][0],oned=True,return_km=True)
      xh = xh[arguments.xmin:arguments.xmax]
      yh = yh[arguments.ymin:arguments.ymax]


      for platforms in range(1,arguments.nplat+1):
         platform_name = 'radar_%03d'%platforms
         tmp_fcst_sngl = radar.obs[platform_name][arguments.var]['obs'][:,arguments.ymin:arguments.ymax,arguments.xmin:arguments.xmax]
         tmp_fcst_sngl = np.where(np.isnan(tmp_fcst_sngl),0.,tmp_fcst_sngl)
         #print('Reading in ...',platform_name)
         if arguments.tilt < 0:
            tmp_fcst_sngl = np.nanmax(tmp_fcst_sngl,axis=0)
         else:
            tmp_fcst_sngl = tmp_fcst_sngl[arguments.tilt]

         if platforms == 1:
            tmp_fcst = np.zeros(tmp_fcst_sngl.shape)
         tmp_fcst = np.where(tmp_fcst>tmp_fcst_sngl,tmp_fcst,tmp_fcst_sngl)


      if mindex == 0 :
          [nyy,nxx] = tmp_fcst.shape
          fcst = np.zeros((arguments.nmem,nyy,nxx))
          fcst[mindex] = tmp_fcst
      else:
          fcst[mindex] = tmp_fcst #radar.obs[platform_name][arguments.var]['obs'][arguments.tilt]





   if arguments.pmmean:  #---Probabilist matched mean
      var = pydart.verification.pmmean(fcst)
      #var = np.percentile(fcst,90,axis=0)
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
         if eindex == 0:
            sample_climo_ens = np.zeros((arguments.nen,sample_climo.shape[0]))
            
            no_skill_ens = np.zeros((arguments.nen,no_skill.shape[0]))
            obs_frequency_ens = np.zeros((arguments.nen,obs_frequency.shape[0]))
            bin_centers_ens = np.zeros((arguments.nen,bin_centers.shape[0]))
            bin_climo_ens = np.zeros((arguments.nen,bin_climo.shape[0]))
            reliability_bins_ens = np.zeros((arguments.nen,reliability_bins.shape[0]))
            bin_count_ens = np.zeros((arguments.nen,bin_count.shape[0]))
            bss_ens = np.zeros((arguments.nen))
         sample_climo_ens[eindex] = sample_climo
         no_skill_ens[eindex] = no_skill
         obs_frequency_ens[eindex] = obs_frequency
         bin_centers_ens[eindex] = bin_centers
         bin_climo_ens[eindex] = bin_climo
         reliability_bins_ens[eindex] = reliability_bins
         bin_count_ens[eindex] = bin_count
 
         figure = plt.figure
         BSS = pydart.verification.BSS(var,obs,arguments.threshold)
         #BSS,pod,pofd = pydart.verification.auc(var,obs,arguments.threshold)
         bss_ens[eindex] = BSS
         pydart.plotting.plotReliability(sample_climo,no_skill,obs_frequency,bin_centers,bin_climo,linecolor = 'r',label='BSS = (%0.3f)'%(BSS))
         plt.legend()
         plt.savefig('Reliability_%s_ens%03d_%3.2f_%3.2f_%s_%02dtilt.png'%(arguments.time,ens,arguments.kernel,arguments.threshold,exp_tag,arguments.tilt),dpi=300)
         plt.clf()
      elif arguments.auc:
         var_auc,pod,pofd = pydart.verification.auc(var,obs,arguments.threshold) 
         print('The AUC for ens%03d is .. %3.3f'%(ens,var_auc))

   if arguments.reliability or arguments.auc:
      pass
   else:
      if arguments.limit:
         outname ='%s_%s_ens%03d_%3.2f_%3.2f_%s_%s_%02dtilt_limit_poorcov.png'%(varname,arguments.time,ens,arguments.kernel,arguments.threshold,exp_tag,mp_tag,arguments.tilt)
      else:
         outname ='%s_%s_ens%03d_%3.2f_%3.2f_%s_%s_%02dtilt.png'%(varname,arguments.time,ens,arguments.kernel,arguments.threshold,exp_tag,mp_tag,arguments.tilt) 
      if arguments.pmmean or arguments.mean:
         pydart.plotting.rough_plot(var,varname,outname=outname)
         if arguments.save:
            np.savez('mean_%02d_%02s.png'%(arguments.tstart,arguments.tstop),field=var,varname=varname,xh=xh,yh=yh)
      else:
         pydart.plotting.rough_plot(var,varname,outname=outname,contour_var=obs,threshold=arguments.threshold)
         if arguments.save:
            np.savez('NMEP_%s_thresh%3.2f_kern%3.2f.npz'%(arguments.time,arguments.threshold,arguments.kernel),field=var,obs=obs,xh=xh,yh=yh)

      plt.clf()


if arguments.reliability:
   #figure = plt.figure()
   figure = plt.figure(figsize=(14,11))
   cols = ['r','b','g','c','y','orangered','m','lime','cornflowerblue','tan']
   for eindex,ens in enumerate(range(1,arguments.nen+1)):
      pydart.plotting.plotReliability(np.mean(sample_climo_ens,axis=0),np.mean(no_skill_ens,axis=0),obs_frequency_ens[eindex],
                                      np.mean(bin_centers_ens,axis=0),np.mean(bin_climo_ens,axis=0),linecolor = cols[eindex],label='ENS%02d (%0.3f)'%(ens,bss_ens[eindex]),linewidth=2.5)
      #pydart.plotting.plotReliability(np.mean(sample_climo_ens[eindex],axis=0),np.mean(no_skill_ens[eindex],axis=0),obs_frequency_ens[eindex],
      #                                np.mean(bin_centers_ens[eindex],axis=0),np.mean(bin_climo_ens[eindex],axis=0),linecolor = cols[eindex],label='BSS = (%0.3f)'%(bss_ens[eindex]))
   pydart.plotting.plotReliability(np.mean(sample_climo_ens,axis=0),np.mean(no_skill_ens,axis=0),np.mean(obs_frequency_ens,axis=0),
                                np.mean(bin_centers_ens,axis=0),np.mean(bin_climo_ens,axis=0),linecolor = 'k',label='MEAN (%0.3f)'%(np.mean(bss_ens)),linewsith=5)

   plt.xlabel('Forecast Probability')
   plt.ylabel('Observed Frequency')
   plt.legend(prop={"size":10})
   plt.savefig('Reliability_%s_ensembles_%3.2f_%3.2f_%s_%s_%02d.png'%(arguments.time,arguments.kernel,arguments.threshold,mp_tag,exp_tag,arguments.tilt),dpi=300)
