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
    rlon1 = N.deg2rad(lon1)
    rlon2 = N.deg2rad(lon2)
    rlat1 = N.deg2rad(lat1)
    rlat2 = N.deg2rad(lat2)
  else:
    rlon1 = lon1
    rlon2 = lon2
    rlat1 = lat1
    rlat2 = lat2

# Simple Lat-Lon grid

  if proj == 'latlon':
    rearth  = 1000.0 * 6367.0
    x       = rearth * N.cos(0.5*(rlat1+rlat2)) * (rlon2-rlon1)
    y       = rearth * (rlat2-rlat1)

# Lambert Conformal

  if proj == 'lcc':
    p1 = Proj(proj='lcc', ellps='WGS84', datum='WGS84', lat_1=truelat1, lat_2=truelat2, lat_0=lat1, lon_0=lon1)
    x, y = p1(lon2, lat2, errchk = True)

  if azimuth:
    ay = N.sin(rlon2-rlon1)*N.cos(rlat2)
    ax = N.cos(rlat1)*N.sin(rlat2)-N.sin(rlat1)*N.cos(rlat2)*N.cos(rlon2-rlon1)
    az = N.degrees(N.arctan2(ay,ax))
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
    rlon1 = N.deg2rad(lon1)
    rlat1 = N.deg2rad(lat1)
  else:
    rlon1 = lon1
    rlat1 = lat1

# Simple Lat-Lon grid

  if proj == 'latlon':
    rearth = 1000.0 * 6367.0
    rlat2  = rlat1 + y / rearth
    lon    = N.rad2deg(rlon1 + x / ( rearth * N.cos(0.5*(rlat1+rlat1)) ) )
    lat    = N.rad2deg(rlat2)

# Lambert Conformal

  if proj == 'lcc':
    p1 = Proj(proj='lcc', ellps='WGS84', datum='WGS84', lat_1=truelat1, lat_2=truelat2, lat_0=lat1, lon_0=lon1)
    lon, lat = p1(x, y, inverse = True)

  if degrees == False:
    return N.deg2rad(lat), N.deg2rad(lon)
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

        elvrad = N.arctan((hgtdb*N.cos(rngdb) - frthrde)/(hgtdb * N.sin(rngdb)))

        return N.rad2deg(elvrad)

    else:

        return -999.
#===============================================================================

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

def rad_obs_loc(fcst,locx,locy,rdr_hgt,tilt,min_range=3000.,max_range=150000.):
   """
   Define the location of different radar locations for a single tilt
  
   locx: The x-location of the radar [float]   (m)
   locy: The y-location of the radar [float]   (m)
   rdr_hgt: The height of the radar  [float]   (m)
   tilts:   The radar tilts          [float] (radians)

   optional:
      min_range: Minimum Distance of radar observations
      max_range: Maximum Distance of radar observations 

   returns:
      obx:  The x-location of the radar observation   [ntilt,ny,nx] (m)
      oby:  The y-location of the radar observation   [ntilt,ny,nx] (m)
      obz:  The height of the radar observation       [ntilt,ny,nx] (m)
      tilt: The radar tilt angle                      [ntilt,ny,nx] (radians)
      az:   Observation azimuth angle                 [ntilt,ny,nx] (radians)
   
   """
   import numpy as np
   #--- Earth Information
   ke = 4./3.
   erad = 6378.14 * 1000.

   #--- Defining Forecast Grid Information 
   irange = np.arange(0,max_range+100.,100.) # Radar Range
   nx = fcst['xh'].shape[0]
   ny = fcst['yh'].shape[0]
   xh = np.tile(fcst['xh'], (ny,1))
   yh = np.transpose(np.tile(fcst['yh'], (nx,1)))

   #--- Calculating the Distance from the Radar
   rad_dis = np.sqrt(((xh-locx)**2) + ((yh-locy)**2))

   #--- Calculating the Radar Azimuth Angle
   azimuth = np.arctan((xh-locx)/(yh-locy))
   azimuth = np.where((yh-locy)<0,azimuth+np.pi,azimuth)
   azimuth = azimuth * (180./np.pi)
   azimuth = np.where(azimuth < 0, 360 + azimuth,azimuth)

   #--- Calculating Radar Height Based Upon Distance from Radar
   #obz = np.zeros((nlev,ny,nx))
   #obz[:,:,:] = np.nan
   #az = np.zeros(obz.shape)
   #obx = np.zeros(az.shape)
   #oby = np.zeros(az.shape)

   #obz = np.zeros((ny,nx))
   #obz[:,:] = np.nan
   #az = np.zeros(obz.shape)
   #obx = np.zeros(az.shape)
   #oby = np.zeros(az.shape)

   #for ntilt,radtilt in enumerate(tilts):
   #   print("Evalulated tilt = ",radtilt)
   #   beam_hgt =  ((irange**2) + (ke*erad)**2 + (2*irange*ke*erad*np.sin(radtilt)))**0.5 - (ke * erad)+rdr_hgt
   #   min_hgt = ((min_range**2) + (ke*erad)**2 + (2*min_range*ke*erad*np.sin(radtilt)))**0.5 - (ke * erad)+rdr_hgt
   #   max_hgt = ((max_range**2) + (ke*erad)**2 + (2*max_range*ke*erad*np.sin(radtilt)))**0.5 - (ke * erad)+rdr_hgt
   #   circle_rad = ke*erad*np.arcsin((irange*np.cos(radtilt))/((ke*erad)+beam_hgt))

   #   #--- Loop over the different radii to obtain radar height
   #   for hindex in range(0,circle_rad.shape[0]-1):
   #      obz[ntilt] = np.where((rad_dis >= circle_rad[hindex]) & (rad_dis < circle_rad[hindex+1]),beam_hgt[hindex],obz[ntilt])
   #   obz[ntilt] = np.where((obz[ntilt] > max_hgt) | (obz[ntilt] < min_hgt),np.nan,obz[ntilt])
   #   az  = np.where(np.isnan(obz[ntilt]),np.nan,azimuth)
   #   obx    = np.where(np.isnan(obz[ntilt]),np.nan,xh)
   #   oby    = np.where(np.isnan(obz[ntilt]),np.nan,yh)
   #   tilt    = np.where(np.isnan(obz[ntilt]),np.nan,radtilt)


   beam_hgt =  ((irange**2) + (ke*erad)**2 + (2*irange*ke*erad*np.sin(tilt)))**0.5 - (ke * erad)+rdr_hgt
   min_hgt = ((min_range**2) + (ke*erad)**2 + (2*min_range*ke*erad*np.sin(tilt)))**0.5 - (ke * erad)+rdr_hgt
   max_hgt = ((max_range**2) + (ke*erad)**2 + (2*max_range*ke*erad*np.sin(tilt)))**0.5 - (ke * erad)+rdr_hgt
   circle_rad = ke*erad*np.arcsin((irange*np.cos(tilt))/((ke*erad)+beam_hgt))

      #--- Loop over the different radii to obtain radar height
   obz = np.zeros((ny,nx))
   for hindex in range(0,circle_rad.shape[0]-1):
      obz = np.where((rad_dis >= circle_rad[hindex]) & (rad_dis < circle_rad[hindex+1]),beam_hgt[hindex],obz)
   obz    = np.where((obz > max_hgt) | (obz < min_hgt),np.nan,obz)
   az     = np.where(np.isnan(obz),np.nan,azimuth).flatten()
   obx    = np.where(np.isnan(obz),np.nan,xh).flatten()
   oby    = np.where(np.isnan(obz),np.nan,yh).flatten()
   tilt   = np.where(np.isnan(obz),np.nan,tilt).flatten()
   obz    = obz.flatten()

      #f = interpolate.interp1d(s, h)
      #
      #for ii  in range(0,nx):
      #   for jj in range(0,ny):
      #      if r[jj,ii] <= np.amax(s):
      #         hgt = f(r[jj,ii])
      #         if hgt > min_hgt and hgt < max_hgt:
      #            obz[ntilt,jj,ii] = f(r[jj,ii])
      #            az[ntilt,jj,ii] = azimuth[jj,ii]
      #         else:
      #            obz[ntilt,jj,ii] = np.nan
      #            az[ntilt,jj,ii] = np.nan

   return obx,oby,obz,tilt,az

#=====================================================================================================
def cressman(x, y, obs, x0, y0, roi, missing=-999999999.):
   """ 
   Returns a data value for the point
   Arguments: x/y/obs:  1D arrays of location
              x0, y0:   point to analyze to
              roi:      radius of influence
              missing:  value to assign if no data, default = 0.0

   This routine can also be implemented in FORTRAN much quicker...
   """
   import numpy as np
   # Create distance array

   nobs = np.size(y0)
   new_ob = np.zeros((nobs))
   for ob in  range(0,nobs):
     #if ob%1000 == 0:
     #   print('Working on observation =',ob)
     if np.isnan(x0[ob]): continue
     dis = np.sqrt( (x-x0[ob])**2 + (y-y0[ob])**2 )
     # Cut the size o n by choosing a smart threshold (2 km)
     indices = np.where(dis <= roi)
     # Check to see if there are any data ponts that are within search radius
     size = np.size(indices)
     # go thru values w/in radius: calculate weight, mult by value and sum also sum weights
     if size != 0:
       R2    = roi**2.0
       w_sum = 0.0
       top   = 0.0
       rk2 = dis[indices]**2.
       wk = (R2-rk2)/(R2+rk2)
       w_sum = np.sum(wk)
       #--- First Observation Called 
       top = np.sum(wk*obs[indices])
       if w_sum < 0.01:
          new_ob[ob] = missing      
       else:
          new_ob[ob] = top/w_sum


     else:  #if there are no values assign the data point a NaN
       new_ob[ob] = np.nan
   return new_ob
#     for n, value in enumerate(dis[indices]):
#        rk2 = value**2.0
#        wk = (R2-rk2) / (R2+rk2)
#        top = top + wk*obs[n]
#        w_sum = w_sum + wk
#       print n, value, data[n], R2, w_sum, top
#
#    if w_sum < 0.01:
#      return missing
#    else:
#      return top/w_sum
#
#  else:
#    return missing
#===============================================================================

