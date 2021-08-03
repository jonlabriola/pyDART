import matplotlib
import matplotlib.pyplot as plt
import numpy as np
#--- Code Used to generate rough plots of EnKF Output
def rough_plot(var,varname,outname=None):
   """
   JDL - A Rough and dirty plot to show to check results
   """

   if varname =='dbz':
      contours = np.arange(0,70.,1.)
      #cb_ticks = np.arange(0,5000.,5.)
      colormap = plt.get_cmap("jet")
   elif varname == 'vr':
      contours = np.arange(-15,15,0.5)
      colormap = plt.get_cmap("BrBG")
   elif varname == 'az':
      contours = np.arange(0,360,10.)
      #cb_ticks = np.arange(0,5000.,500.)
      colormap = plt.get_cmap("gist_stern")
   else:
      contours = np.arange(0,5000.,500.)
      #cb_ticks = np.arange(0,5000.,500.)
      colormap = plt.get_cmap("gist_stern")
   cmap = plt.get_cmap(colormap)
   norm = matplotlib.colors.BoundaryNorm(contours, ncolors=cmap.N, clip=True)
   CS = plt.pcolormesh(var,cmap=cmap,norm=norm,alpha=1.0,shading='auto',edgecolors='none')
   plt.colorbar(CS)
   if outname is None:
      plt.show()
   else:
      plt.savefig(outname)
      plt.clf()
