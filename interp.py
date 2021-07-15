#---- Scripts are used to interpolate location and distance
#---- and converted between lat/lon and cartesian coordinates
import numpy as np
map_projection = 'latlon'

def dll_2_dxy(lat1, lat2, lon1, lon2, degrees=False, azimuth=False, proj = map_projection):

  """dll_2_dxy returns the approximate distance in meters between two lat/lon pairs

     Valid projections: Lambert Conformal
                        Lat - Lon

     INPUTS: Two (lat,lon) pairs in radians, or if degrees==True, degrees (default)

     if lon2 > lon1:  x > 0

     if lat2 > lat1:  y > 0

     OUTPUTS:  DX, DY in meters

     Azimuth formula from http://www.movable-type.co.uk/scripts/latlong.html
  """

  if degrees:
    rlon1 = np.deg2rad(lon1)
    rlon2 = np.deg2rad(lon2)
    rlat1 = np.deg2rad(lat1)
    rlat2 = np.deg2rad(lat2)
  else:
    rlon1 = lon1
    rlon2 = lon2
    rlat1 = lat1
    rlat2 = lat2

# Simple Lat-Lon grid

  if proj == 'latlon':
    rearth  = 1000.0 * 6367.0
    x       = rearth * np.cos(0.5*(rlat1+rlat2)) * (rlon2-rlon1)
    y       = rearth * (rlat2-rlat1)

# Lambert Conformal

  if proj == 'lcc':
    p1 = Proj(proj='lcc', ellps='WGS84', datum='WGS84', lat_1=truelat1, lat_2=truelat2, lat_0=lat1, lon_0=lon1)
    x, y = p1(lon2, lat2, errchk = True)

  if azimuth:
    ay = np.sin(rlon2-rlon1)*np.cos(rlat2)
    ax = np.cos(rlat1)*np.sin(rlat2)-np.sin(rlat1)*np.cos(rlat2)*np.cos(rlon2-rlon1)
    az = np.degrees(np.arctan2(ay,ax))
    return x, y, az

  return x, y


#===============================================================================
def dxy_2_dll(x, y, lat1, lon1, degrees=True, proj = map_projection):

  """dxy_2_dll returns the approximate lat/lon between an x,y coordinate and
     a reference lat/lon point.

     Valid projections: Lambert Conformal
                        Lat - Lon

     INPUTS:  x,y in meters, lat1, lon1 in radians, or if degrees == True,
              then degrees (default value)

     if x > 0, lon > lon1

     if y > 0, lat > lat1

     OUTPUTS:  lat, lon in radians, or if degrees == True, degrees
               i.e., the input and output units for lat/lon are held the same.
  """

  if degrees:
    rlon1 = np.deg2rad(lon1)
    rlat1 = np.deg2rad(lat1)
  else:
    rlon1 = lon1
    rlat1 = lat1

# Simple Lat-Lon grid

  if proj == 'latlon':
    rearth = 1000.0 * 6367.0
    rlat2  = rlat1 + y / rearth
    lon    = np.rad2deg(rlon1 + x / ( rearth * np.cos(0.5*(rlat1+rlat1)) ) )
    lat    = np.rad2deg(rlat2)

# Lambert Conformal

  if proj == 'lcc':
    p1 = Proj(proj='lcc', ellps='WGS84', datum='WGS84', lat_1=truelat1, lat_2=truelat2, lat_0=lat1, lon_0=lon1)
    lon, lat = p1(x, y, inverse = True)

  if degrees == False:
    return np.deg2rad(lat), np.deg2rad(lon)
  else:
    return lat, lon
#=================================================================
def beam_elv(sfc_range, z):
########################################################################
#
#     PURPOSE:
#
#     Calculate the elevation angle (elvang) and the along
#     ray-path distance (range) of a radar beam
#     crossing through the given height and along-ground
#     distance.
#
#     This method assumes dn/dh is constant such that the
#     beam curves with a radius of 4/3 of the earth's radius.
#     This is dervied from Eq. 2.28 of Doviak and Zrnic',
#     Doppler Radar and Weather Observations, 1st Ed.
#
########################################################################
#
#     AUTHOR: Keith Brewster
#     10/10/95
#
#     MODIFICATION HISTORY: adapted to python by Lou Wicker (thanks Keith)
#
########################################################################
#
#     INPUT:
#       sfc_range:    Distance (meters) along ground from radar
#       z        :    Height above radar
#
#     OUTPUT
#       elvang   Elevation angle (degrees) of radar beam
#
########################################################################
    eradius=6371000.
    frthrde=(4.*eradius/3.)
    eighthre=(8.*eradius/3.)
    fthsq=(frthrde*frthrde)

    if sfc_range > 0.0:
        hgtdb = frthrde + z
        rngdb = sfc_range/frthrde

        elvrad = np.arctan((hgtdb*np.cos(rngdb) - frthrde)/(hgtdb * np.sin(rngdb)))

        return np.rad2deg(elvrad)

    else:

        return -999.
#===============================================================================
def interp_wghts(x, xc, extrapolate=False):
  """interp_weights returs the linear interpolation weights for a given
     ascending array of coordinates.

     x = location to be interpolated to
     xc = locations of grid
     extrapolate:  what to do at the boundary?
                   for now, if at edge, set return values to missing

     OUTPUTS:  i0, i1, dx0, dx1 locations and weights for the interpolation
  """

  indices = np.where(xc <= x)

  #print('indicies = ',indices)
  #print('x =',x)
  #print('xc = ',xc)

  if np.size(indices[0]) == 0:
    return -1, -1, None, None, None
  else:
    i0 = indices[0][-1]
    if i0 == np.size(xc)-1:
      return -1, -1, None, None, None
    else:
      i1  = i0 + 1
      dx  = xc[i1] - xc[i0]
      dx0 = xc[i1] - x
      dx1 = dx-dx0
      return i0, i1, dx0, dx1, dx
#===============================================================================

def unstagger_grid(var,axis):
   """
   Program to unstagger grids
   Required Inputs:
   var:  array with an undetermined number of dimensions
   axis: The axis to perform "unstaggering"

   """

   if axis == 0:
      var = (var[1:] + var[:-1]) / 2.
   elif axis == 1:
      var = (var[:,1:] + var[:,:-1]) / 2.
   elif axis == 2:
      var = (var[:,:,1:] + var[:,:,:-1]) / 2.
   elif axis == 3: #--- Unlikely to ever happen
      var = (var[:,:,:,1:] + var[:,:,:,:-1]) / 2.
   return var

#==============================================================================

def rad_obs_loc(model,xloc,yloc,zloc,radtilt,min_range=3000.,max_range=150000.):
#def rad_obs_loc(radobs,min_range=3000.,max_range=150000.):
   """
   Define the location of different radar locations for a radar volume
  
   reads in the radobs class, which contains:
     model:  Model Information including xx,yy [ny,nx] (m)
     xloc: The x-location of the radar  float   (m)
     yloc: The y-location of the radar  float   (m)
     zloc: The height of the radar      float   (m)
     tilts:   The radar tilts          [ntilt]   (radians)

   optional:
      min_range: Minimum Distance of radar observations
      max_range: Maximum Distance of radar observations 

   returns:
      obx:  The x-location of the radar observation   [ntilt,ny,nx] (m)
      oby:  The y-location of the radar observation   [ntilt,ny,nx] (m)
      obz:  The height of the radar observation       [ntilt,ny,nx] (m)
      radtilt: The radar tilt angle                      [ntilt,ny,nx] (radians)
      az:   Observation azimuth angle                 [ntilt,ny,nx] (radians)
   
   """
   import numpy as np
   #--- Earth Information
   ke = 4./3.
   erad = 6378.14 * 1000.

   #--- Defining Forecast Grid Information 
   irange = np.arange(0,max_range+100.,100.) # Radar Range
   xh = model['xh2d'] #---2D Plane of X-Distance Values
   nx = model['nx']
   yh = model['yh2d'] #---2D Plane of Y-Distance Values
   ny = model['ny']
   ntilt = len(radtilt)

   #--- Calculating the Distance from the Radar
   rad_dis = np.sqrt(((xh-xloc)**2) + ((yh-yloc)**2))

   #--- Calculating the Radar Azimuth Angle
   azimuth = np.arctan((xh-xloc)/(yh-yloc))
   azimuth = np.where((yh-yloc)<0,azimuth+np.pi,azimuth)
   azimuth = azimuth * (180./np.pi)
   azimuth = np.where(azimuth < 0, 360 + azimuth,azimuth)
   azimuth = np.radians(azimuth)

   #--- Calculating Beam Height Based Upon Radar Tilt
   obz = np.zeros((ntilt,ny,nx))
   obx = np.zeros(obz.shape)
   oby = np.zeros(obz.shape)
   rtilt = np.zeros(obz.shape)
   az = np.zeros(obz.shape)

   #--- Loop Over the different radar tilts, calculate height
   for ntilt, tilt in enumerate(radtilt):
      beam_hgt =  ((irange**2) + (ke*erad)**2 + (2*irange*ke*erad*np.sin(tilt)))**0.5 - (ke * erad)+zloc
      min_hgt = ((min_range**2) + (ke*erad)**2 + (2*min_range*ke*erad*np.sin(tilt)))**0.5 - (ke * erad)+zloc
      max_hgt = ((max_range**2) + (ke*erad)**2 + (2*max_range*ke*erad*np.sin(tilt)))**0.5 - (ke * erad)+zloc
      circle_rad = ke*erad*np.arcsin((irange*np.cos(tilt))/((ke*erad)+beam_hgt)[ntilt])

      #--- Loop over the different radii to obtain approximate beam height
      for hindex in range(0,circle_rad.shape[0]-1):
         obz[ntilt] = np.where((rad_dis >= circle_rad[hindex]) & (rad_dis < circle_rad[hindex+1]),beam_hgt[hindex],obz[ntilt])
      obz[ntilt]    = np.where((obz[ntilt] > max_hgt) | (obz[ntilt] < min_hgt),np.nan,obz[ntilt])
      az[ntilt]     = np.where(np.isnan(obz[ntilt]),np.nan,azimuth)
      obx[ntilt]    = np.where(np.isnan(obz[ntilt]),np.nan,xh)
      oby[ntilt]    = np.where(np.isnan(obz[ntilt]),np.nan,yh)
      rtilt[ntilt]  = np.where(np.isnan(obz[ntilt]),np.nan,tilt)


   return obx,oby,obz,rtilt,az

#=====================================================================================================
def cressman(x, y, obs, x0, y0, roi, missing=-999999999.):
   """ 
   Returns a data value for the point
   Arguments: x, y, obs    : [nz,ny, nx]
              x0, y0       : [nz,ny_coarse,nx_coarse]
              roi          : radius of influence
              missing:  value to assign if no data, default = 0.0

   This routine can also be implemented in FORTRAN much quicker...
   """
   # Create distance array
   [nz,ny,nx] = x0.shape

   #nobs = np.size(y0)
   new_ob = np.zeros((nz,ny,nx))
   new_ob[:,:,:] = np.nan
   print('new_obs shape = ', new_ob.shape)
   R2    = roi**2.0
   print(np.amax(x0))
   for k in range(0,nz):
      print('Performing cressman on the %03d level'%k)
      for j in range(0,ny):
         for i in range(0,nx):
            
            if np.isnan(x0[k,j,i]): continue #--- Skip NaN's
           

            if np.isnan(x0[k,j,i]):
               print(x0[k,j,i]) 
            dis = np.sqrt( (x[k]-x0[k,j,i])**2 + (y[k]-y0[k,j,i])**2 )
            indices = np.where(dis<= roi)
            size = np.size(indices)
            if size != 0:
               #print('indices = ',indices)
               w_sum = 0.0
               top   = 0.0
               rk2 = dis[indices]**2.
               wk = (R2-rk2)/(R2+rk2)
               w_sum = np.sum(wk)
               top = np.sum(wk*obs[k][indices])
               if w_sum<0.01:
                  new_ob[k,j,i] = missing
               else:
                  new_ob[k,j,i] = top/w_sum

            else:
               new_ob[k,j,i] = np.nan
   return new_ob

#===============================================================================


#=====================================================================================================
def cressman_orig(x, y, obs, x0, y0, roi, missing=-999999999.):
   """
   Returns a data value for the point
   Arguments: x, y    : [nz,ny, nx]
              obs:      [nz,ny,nx]
              x0, y0:   [nz,ny_coarse,nx_coarse]
              roi:      radius of influence
              missing:  value to assign if no data, default = 0.0

   This routine can also be implemented in FORTRAN much quicker...
   """
   # Create distance array

   nobs = np.size(y0)
   new_ob = np.zeros((nobs))
   R2    = roi**2.0
   for ob in  range(0,nobs):
     if np.isnan(x0[ob]): continue
     dis = np.sqrt( (x-x0[ob])**2 + (y-y0[ob])**2 )
     # Cut the size o n by choosing a smart threshold (2 km)
     indices = np.where(dis <= roi)
     # Check to see if there are any data ponts that are within search radius
     size = np.size(indices)
     # go thru values w/in radius: calculate weight, mult by value and sum also sum weights
     if size != 0:
       w_sum = 0.0
       top   = 0.0
       rk2 = dis[indices]**2.
       wk = (R2-rk2)/(R2+rk2)
       w_sum = np.sum(wk)
       top = np.sum(wk*obs[indices])
       if w_sum < 0.01:
          new_ob[ob] = missing
       else:
          new_ob[ob] = top/w_sum
     else:  #if there are no values assign the data point a NaN
       new_ob[ob] = np.nan
   return new_ob

#===============================================================================
