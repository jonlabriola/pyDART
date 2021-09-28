import matplotlib
import matplotlib.pyplot as plt
import numpy as np
#--- Code Used to generate rough plots of EnKF Output
def get_colors(varname):
   cb_ticks = None
   if varname.lower() in ['dbz','radar_reflectivity', 'fine_z']:
      #contours = np.arange(0,70.,1.)
      #cb_ticks = np.arange(0,75,5.)
      #colormap = plt.get_cmap("jet")
      contours = np.arange(0.,75.,5.)
      colormap = ['#FFFFFF', '#00FFFF', '#6495ED', '#0000CC', '#00FF00', '#00BB00',
                   '#008800', '#FFFF00', '#FFDD00', '#DAA520', '#FF0000', '#B21111',
                   '#990000', '#FF00FF', '#BB55DD']
      cb_ticks = [10., 20., 30., 40., 50., 60., 70.]
      colormap = matplotlib.colors.ListedColormap(colormap,cb_ticks)
   elif varname.lower() in ['vr','doppler_radial_velocity']:
      contours = np.arange(-25,26,1.0)
      colormap = plt.get_cmap("BrBG")
      cb_ticks = np.arange(-25,26,2) 
   elif varname == 'az':
      contours = np.arange(0,360,10.)
      colormap = plt.get_cmap("gist_stern")
   else:
      contours = np.arange(0,5000.,500.)
      colormap = plt.get_cmap("gist_stern")
   if cb_ticks is None:
      cb_ticks = contours[::2]
   cmap = plt.get_cmap(colormap)
   norm = matplotlib.colors.BoundaryNorm(contours, ncolors=cmap.N, clip=True)

   return cmap,norm,cb_ticks

def rough_plot(var,varname,outname=None):
   """
   JDL - A Rough and dirty plot to show to check results
   """
   colormap = get_colors(varname)
   CS = plt.pcolormesh(var,cmap=cmap,norm=norm,alpha=1.0,shading='auto',edgecolors='none')
   plt.colorbar(CS)
   if outname is None:
      plt.show()
   else:
      plt.savefig(outname)
      plt.clf()

def twod_plot(var_class,varname,tilt=1,outname=None,copy_name=None):
   cmap,norm,cb_ticks = get_colors(varname)
   if copy_name is None:
      copy_name = 'obs'


   if np.isnan(np.amax(var_class['xloc'])):
      dx = np.nan
      j = 1
      while np.isnan(dx):
        try:
           dx = np.nanmin(var_class['xloc'][tilt,j,1:] - var_class['xloc'][tilt,j,0:-1])/1000.
           j+=1
        except:
           print('end of the line no more grid points')
      #dx = np.nanmin(var_class['xloc'][tilt,0,1:] - var_class['xloc'][tilt,1:,0:-1])
      [nz,ny,nx] = var_class['xloc'].shape 
      xh = np.arange(0,nx*dx,dx)
      yh = np.arange(0,ny*dx,dx)
      xx,yy = np.meshgrid(xh,yh)  
   else:
      xx,yy = np.meshgrid(var_class['xloc'][tilt,0,:]/1000.,var_class['yloc'][tilt,:,0]/1000.)
 
   zh = var_class['zloc'][tilt,:,:]
   var_class[copy_name][tilt] = np.where(zh < 150, np.nan, var_class[copy_name][tilt])  
    
   print('Var max = ',np.nanmax(var_class[copy_name][tilt]))
   print('Var min = ',np.nanmin(var_class[copy_name][tilt]))
   CS = plt.pcolormesh(xx,yy,var_class[copy_name][tilt],cmap=cmap,norm=norm,alpha=1.0,shading='auto',edgecolors='none')
   CB = plt.colorbar(CS, shrink=0.8, ticks = cb_ticks) 
   if outname is None:
      plt.show()
   else:
      plt.savefig(outname)
      plt.clf()

