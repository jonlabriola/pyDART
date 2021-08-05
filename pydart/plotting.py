import matplotlib
import matplotlib.pyplot as plt
import numpy as np
#--- Code Used to generate rough plots of EnKF Output
def get_colors(varname):
   cb_ticks = None
   if varname.lower() in ['dbz','radar_reflectivity', 'fine_z']:
      contours = np.arange(0,70.,1.)
      cb_ticks = np.arange(0,75,5.)
      colormap = plt.get_cmap("jet")
   elif varname == 'vr':
      contours = np.arange(-15,15,0.5)
      colormap = plt.get_cmap("BrBG")
      cb_ticks = np.arange(-15,15,2) 
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
   #cmap = plt.get_cmap(colormap)
   #norm = matplotlib.colors.BoundaryNorm(contours, ncolors=cmap.N, clip=True)
   CS = plt.pcolormesh(var,cmap=cmap,norm=norm,alpha=1.0,shading='auto',edgecolors='none')
   plt.colorbar(CS)
   if outname is None:
      plt.show()
   else:
      plt.savefig(outname)
      plt.clf()

def twod_plot(var_class,varname,tilt=1,outname=None):
   cmap,norm,cb_ticks = get_colors(varname)
   #print(var_class['xloc'].shape)
   #xx,yy  = np.meshgrid(var_class['xloc'][tilt],var_class['yloc'][tilt])
   #CS = plt.pcolormesh(var_class['xloc'][tilt],var_class['yloc'][tilt],var_class['obs'][tilt],cmap=cmap,norm=norm,alpha=1.0,shading='auto',edgecolors='none')

   if np.isnan(np.amax(var_class['xloc'])):
      dx = np.nanmin(var_class['xloc'][tilt,0,1:] - var_class['xloc'][tilt,1:,0:-1])
      [nz,ny,nx] = var_class['xloc'].shape 
      xh = np.arange(0,nx*dx,dx)
      yh = np.arange(0,ny*dx,dx)
      xx,yy = np.meshgrid(xh,yh)  
   else:
      xx,yy = np.meshgrid(var_class['xloc'][tilt,0,:],var_class['yloc'][tilt,:,0])
  
   #print(xx[4:6,4:6])
   #print(var_class['xloc'][tilt,4:6,4:6])
   #print(np.amax(xx))
   CS = plt.pcolormesh(xx,yy,var_class['obs'][tilt],cmap=cmap,norm=norm,alpha=1.0,shading='auto',edgecolors='none')
   CB = plt.colorbar(CS, shrink=0.8, ticks = cb_ticks) 
   if outname is None:
      plt.show()
   else:
      plt.savefig(outname)
      plt.clf()

