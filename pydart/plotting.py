import matplotlib
import matplotlib.pyplot as plt
import numpy as np
#--- Code Used to generate rough plots of EnKF Output
def get_colors(varname):
   cb_ticks = None
   if varname.lower() in ['dbz','radar_reflectivity', 'fine_z']:
      contours = np.arange(0,70.,0.5)
      cb_ticks = np.arange(0,75,5.)
      colormap = plt.get_cmap("jet")
      #contours = np.arange(0.,75.,5.)
      #colormap = ['#FFFFFF', '#00FFFF', '#6495ED', '#0000CC', '#00FF00', '#00BB00',
      #             '#008800', '#FFFF00', '#FFDD00', '#DAA520', '#FF0000', '#B21111',
      #             '#990000', '#FF00FF', '#BB55DD']
      #cb_ticks = [10., 20., 30., 40., 50., 60., 70.]
      #colormap = matplotlib.colors.ListedColormap(colormap,cb_ticks)
   elif varname.lower() in ['vr','doppler_radial_velocity']:
      contours = np.arange(-25,26,1.0)
      colormap = plt.get_cmap("BrBG")
      cb_ticks = np.arange(-25,26,2) 
   elif varname == 'az':
      contours = np.arange(0,360,10.)
      colormap = plt.get_cmap("gist_stern")
   elif varname.lower() in ['nep','nmep']:
      contours = np.array([0.0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1.0])
      colormap = ['#FFFFFF', '#DDF0FF', '#CCE0FF', '#88AAFF', '#4477FF', '#FFFF99', '#F0F000',
                        '#C0C000', '#FF7777', '#FF2222', '#CC0000', '#AA0000']
      cb_ticks = np.array([0.0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1.0])
      colormap = matplotlib.colors.ListedColormap(colormap,cb_ticks)
   else:
      contours = np.arange(0,5000.,500.)
      colormap = plt.get_cmap("gist_stern")
   if cb_ticks is None:
      cb_ticks = contours[::2]
   cmap = plt.get_cmap(colormap)
   norm = matplotlib.colors.BoundaryNorm(contours, ncolors=cmap.N, clip=True)

   return cmap,norm,cb_ticks

def rough_plot(var,varname,outname=None,**kwargs):
   """
   JDL - A Rough and dirty plot to show to check results
   """
   cmap,norm,cb_ticks = get_colors(varname)
   CS = plt.pcolormesh(var,cmap=cmap,norm=norm,alpha=1.0,shading='auto',edgecolors='none')
   plt.colorbar(CS)
  
   #--- Add on Contours (If Desired)
   if 'contour_var' in kwargs:
      if 'threshold' in kwargs:
         threshold = [0,kwargs['threshold']]
      else:
         threshold = [0,1]
      plt.contour(kwargs['contour_var'],threshold,linewidth=5.0)

   if outname is None:
      plt.show()
   else:
      print('Saving ...',outname)
      plt.savefig(outname)
      plt.clf()

def twod_plot(var_class,varname,tilt=1,outname=None,copy_name=None):
   cmap,norm,cb_ticks = get_colors(varname)
   print('JDL the tilt is ...',tilt)
   if copy_name is None:
      copy_name = 'obs'

   [nz,ny,nx] = var_class['xloc'].shape 
   if np.isnan(np.amax(var_class['xloc'])):
      dx = np.nan
      j = 1
      while np.isnan(dx):
        print('J = ',j)
        if j > ny: break
        try:
           dx = np.nanmin(var_class['xloc'][tilt,j,1:] - var_class['xloc'][tilt,j,0:-1])/1000.
           j+=1
        except:
           j+=1
           print('J = ',j)
           #pass
           #print('end of the line no more grid points')
      #dx = np.nanmin(var_class['xloc'][tilt,0,1:] - var_class['xloc'][tilt,1:,0:-1])
      print('dx = ',dx)
      if np.isnan(dx):
         xh = np.arange(0,nx,1)
         yh = np.arange(0,ny,1)
         xx,yy = np.meshgrid(xh,yh)
      else:
         [nz,ny,nx] = var_class['xloc'].shape 
         xh = np.arange(0,nx*dx,dx)
         yh = np.arange(0,ny*dx,dx)
         xx,yy = np.meshgrid(xh,yh)  
   else:
      xx,yy = np.meshgrid(var_class['xloc'][tilt,0,:]/1000.,var_class['yloc'][tilt,:,0]/1000.)

   for tindex in range(0,14):
      print('Var min for tindex = ',np.nanmin(var_class[copy_name][tindex]))
      print('Var max for tindex = ',var_class[copy_name][tindex].shape) 
   print('Var max = ',np.nanmax(var_class[copy_name][tilt]))
   print('Var min = ',np.nanmin(var_class[copy_name][tilt]))
   CS = plt.pcolormesh(xx,yy,var_class[copy_name][tilt],cmap=cmap,norm=norm,alpha=1.0,shading='auto',edgecolors='none')
   CB = plt.colorbar(CS, shrink=0.8, ticks = cb_ticks) 
   if outname is None:
      plt.show()
   else:
      plt.savefig(outname)
      plt.clf()


def gen_performance():
   """
   This function can be used to generate a simple performance diagram.
   Once the diaggram is created you can add your own FOH/POD pairs
  
   Required Inputs:
      None

   Returns: Figure Object to be used for plotting
   """
   figure = plt.figure(figsize=(14,11))
   x=y= np.arange(0,1.1,.01)
   X,Y = np.meshgrid(x,y)
   CSI = ((1/X)+(1/Y)-1)**-1
   #CS = plt.contourf(X,Y,CSI,np.arange(0.1,1.2,.2),colors=['#ffe6e6','#ff9999','#ff6666','#ff3333','#cc0000'])
   CS = plt.contourf(X,Y,CSI,np.arange(0.1,1.2,.2),colors=['#bfbfbf','#a6a6a6','#8c8c8c','#737373','#595959'])
   biases = [0.25,0.7,1,1.5,4]
   for i in biases:
      plt.annotate('%.2f'%i,(.5,.5*i),fontsize=10)
   Bias = Y/X
   BS = plt.contour(X,Y,Bias,biases,colors='black',linestyles='--')
   CB = plt.colorbar(CS, shrink=0.8, ticks = [0.1,0.3,0.5,0.7,0.9])
   CB.set_label('CSI')
   plt.xlim([0,1])
   plt.ylim([0,1])
   return figure


def plotReliability(sample_climo,no_skill,obs_frequency,bin_centers,bin_climo,**kwargs):
   """
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Goal:  Plot standard variables Reliability Curve
   !
   ! Need:  All of the data obtained from Reliability
   !
   ! Produced June 9, 2015
   ! Author: Jon Labriola
   !
   ! Modifications: NONE
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   """

   plt.xticks(np.arange(0.0, 1.2, 0.2))
   plt.yticks(np.arange(0.0, 1.2, 0.2))

   plt.axis([0., 1., 0., 1.])
   if 'linecolor' in kwargs:
     linecolor = kwargs['linecolor']
   else:
     linecolor = 'black'

   if 'alpha' in kwargs:
      alpha = kwargs['alpha']
   else:
      alpha = 1.0

   if 'linewidth' in kwargs:
     linewidth = kwargs['linewidth']
   else:
     linewidth=5.

   #Need sample_climo, no_skill
   plt.fill_between(np.arange(0, 1.001, 0.001), no_skill, np.ones((1001)), where=no_skill>sample_climo, edgecolor='#CCCCCC', facecolor='#CCCCCC', interpolate=True)
   plt.fill_between(np.arange(0, 1.001, 0.001), np.zeros((1001)), no_skill, where=no_skill<sample_climo, edgecolor='#CCCCCC', facecolor='#CCCCCC', interpolate=True)

   #Need obs_frequency
   frequency_masked = np.ma.masked_where(obs_frequency < 0, obs_frequency)

   #--- Plotting
   if 'label' in kwargs:
      plt.plot(bin_centers, frequency_masked, color = linecolor, linewidth = linewidth,alpha=alpha,label=kwargs['label'])
   else:
      plt.plot(bin_centers, frequency_masked, color = linecolor, linewidth = linewidth,alpha=alpha) #label = 'ROC= '+'%.3f'%ROC)

   #Perfect reliability line
   #Need bin_centers
   plt.plot(bin_centers, bin_centers, color='#000000', linewidth=1.5, linestyle='--')

   #No resolution line (sample climatology)
   #bin_climo and bin_centers
   plt.plot(bin_centers, bin_climo, color='#000000', linewidth=1.5, linestyle='--')


   if 'outpath' in kwargs:
     print("Plot saved to " + str(kwargs['outpath']))
     plt.savefig(kwargs['outpath'], figsize = (13, 13), dpi=300)
     plt.clf()


def calc_grid_dims(xloc,yloc,oned=False,return_km=False):
   """
  Calculate the locations of the horizontal grid points

  Inputs:
     xloc [ny,nx]  The x-grid information (can be computed from obs sequence)
     yloc [ny,nx]  The y-grid information (can be computed from obs sequence)
  

  Optional:
     oned: A flag to return one-dimension grid information similar to xh,yh
     return_km: A flag to convert units to km

  Returns

     xx,yy [ny,nx] the horizontal grid point locations
   """
   if np.isnan(np.amax(xloc)) or np.isnan(np.amax(yloc)):
      dx = np.nan
      j = 1
      while np.isnan(dx):
         try:
            dx = np.nanmin(xloc[j,1:] - xloc[j,0:-1])
            j+=1
         except:
           print('end of the line no more grid points')
      [ny,nx] = xloc.shape
      xh = np.arange(0,nx*dx,dx)
      yh = np.arange(0,ny*dx,dx)

      if oned:
         xx = xh
         yy = yh
      else:
         xx,yy = np.meshgrid(xh,yh)
   else:
      if oned:
         xx = xloc[0,:]
         yy = yloc[:,0]
      else:
         xx,yy = np.meshgrid(xloc[0,:],yloc[:,0])

   if return_km:
      return xx/1000.,yy/1000.
   else:
      return xx,yy
