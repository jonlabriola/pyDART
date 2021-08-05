import pydart,argparse,pickle
parser = argparse.ArgumentParser()
parser.add_argument("var",type=str,help = 'Var to plot (RADAR_REFLECTIVITY, fine_z, etc.)')
parser.add_argument("file_path",type=str,help='The location of the dmax file')
parser.add_argument("--platform",type=int,default=1,help='The platform number (default = 1)')
parser.add_argument("--tilt",type=int,default=1,help='Radar tilt (default = 1)')
arguments = parser.parse_args()

radar = pickle.load(open(arguments.file_path,"rb"))

platform_name = 'radar_%03d'%arguments.platform

plotvar = radar.obs[platform_name][arguments.var]
outname = 'radar_%s_%03d_%03dtilt.png'%(arguments.var,arguments.platform,arguments.tilt)
pydart.plotting.twod_plot(plotvar,arguments.var,tilt=1,outname=outname)

