#--- Code for Forward Operators Used to Calculate Observations
import pydart
import numpy as np
import math

#ens_default_cmap = P.get_cmap('seismic')

_missing = -9900000000.
_dbz_min = 0.0
_dbz_max = 75.

_clevels = {'DBZ':[5,75,5], 'W':[-10,11,1], 'TH':[295.,325.,1.], 'THP':[-5,6,0.5], 'QV':[0.0,0.022,0.002], \
            'WZ': [-100.,120.,20.], 'U':[-20,22,2], 'V':[-20,22,2]}

def calcHx(fcst, xloc, yloc, height, elv, azimuth,
           cartesian=False,missing=None,clear_dbz=None,dbzOnly=False):
  """
  The Radar Forward Operator Used to Calculate Z and Radial Wind
  Required Arguments
      1.) fcst    - The  model state variables                [nvars] [nk,ny,nz] - object contains vars and grid info
      3.) yloc    - Observation Y-Location (m or Latitude)    [ntilt,ny,nx] (m, degrees lat)
      4.) xloc    - Observation X-Location (m or Longitude)   [ntilt,ny,nx] (m, degrees lon)  
      5.) height  - Observation Height above Surface (m)      [ntilt,ny,nx] (m)
      6.) elv    - Radar elevation angle                      [ntilt,ny,nx] (radians)
      7.) azimuth - radar azimuth angle                       [ntilt,ny,nx] (radians)


  Optional Arguments:
     cartesian    -  Are we working with a cartesian grid    [bool]
     missing      -  Value of a missing observation          [float]
     clear_dbz    -  What is defined as clear air (remove 
                     vr obs accordingly                      [float] 
     dbzOnly      -  Only calculate reflectivity, no Vr      [bool] 
  """
  nobs = np.size(yloc) 
  [ntilt,ny,nx] = yloc.shape
  if missing == None:
    missing = _missing

  if dbzOnly:
     Hx_Z  = missing * np.ones((ntilt,ny,nx))
  else:
     Hx_vr = missing * np.ones((ntilt,ny,nx))
     Hx_Z  = missing * np.ones((ntilt,ny,nx))

  if cartesian: #--- x/y/z coordinates
     yloc = yloc
     xloc = xloc

     ixloc = 0.0
     iyloc = 0.0
     ihgt = 0.0

     mxlocs = fcst['xh']
     mylocs = fcst['yh']
     mhgts = fcst['zh']

     print(' calcHx:  OBS   Y-DIS MIN/MAX:  ', np.nanmin(yloc),   np.nanmax(yloc))
     print(' calcHx:  MODEL Y-DIS MIN/MAX:  ', np.nanmin(mylocs), np.nanmax(mylocs))
     print(' calcHx:  OBS   X-DIS MIN/MAX:  ', np.nanmin(xloc),   np.nanmax(xloc))
     print(' calcHx:  MODEL X-DIS MIN/MAX:  ', np.nanmin(mxlocs), np.nanmax(mxlocs))

     print(' calcHx:  OBS   HGT MIN/MAX:  ', np.nanmin(height), np.nanmax(height))
     print(' calcHx:  MODEL HGT MIN/MAX:  ', np.nanmin(mhgts),  np.nanmax(mhgts))

  else: #--- Lat/Lon Coordinates
     lat = yloc
     lon = xloc

     ilon = 0.0
     ilat = 0.0
     iihgt = 0.0

     mlons = fcst['LON']  #--- JDL Will Need to UPDATE 
     mlats = fcst['LAT']  #--- JDL Will Need To UPDATE
     mhgts = fcst['zh']

     print(' calcHx:  OBS   LAT MIN/MAX:  ', np.nanmin(lat),    np.nanmax(lat))
     print(' calcHx:  MODEL LAT MIN/MAX:  ', np.nanmin(mlats),  np.nanmax(mlats))
     print(' calcHx:  OBS   LON MIN/MAX:  ', np.nanmin(lon),    np.nanmax(lon))
     print(' calcHx:  MODEL LON MIN/MAX:  ', np.nanmin(mlons),  np.nanmax(mlons))

     print(' calcHx:  OBS   HGT MIN/MAX:  ', np.nanmin(height), np.nanmax(height))
     print(' calcHx:  MODEL HGT MIN/MAX:  ', np.nanmin(mhgts),  np.nanmax(mhgts))

  b = np.zeros([5])

  for k in range(0,ntilt):
    for j in range(0,ny):
      for i in range(0,nx):
         #--- Skip all gridpoints with nan's
         if np.isnan(xloc[k,j,i]):
            if dbzOnly:
               Hx_Z[k,j,i] = np.nan
            else: 
               Hx_vr[k,j,i] = np.nan
               Hx_Z[k,j,i] = np.nan
            continue

         # need to acess lat lon alt info of ob loc from table and pass it to tlint
         # check to see if ob position is the same then less work

         #--- JDL If Statements are used to limit repeitition (the same i1,i2 ... values preserved)
         if cartesian: #--- Use x/y/z Coordinates
            if xloc[k,j,i] != ixloc:
               ixloc = xloc[k,j,i]
               i1, i2, dx1, dx2, dx = pydart.interp.interp_wghts(ixloc, mxlocs)
            if yloc[k,j,i] != iyloc:
               iyloc = yloc[k,j,i]
               j1, j2, dy1, dy2, dy = pydart.interp.interp_wghts(iyloc, mylocs)
            if height[k,j,i] != ihgt:
              ihgt = height[k,j,i]
              k1, k2, dz1, dz2, dz = pydart.interp.interp_wghts(ihgt, mhgts)

         else: #--- Use Lat/Lon Coordinates
            if lon[k,j,i] != ilon:
              ilon = lon[k,j,i]
              i1, i2, dx1, dx2, dx = pydart.interp.interp_wghts(ilon, mlons)
            if lat[k,j,i] != ilat:
              ilat = lat[k,j,i]
              j1, j2, dy1, dy2, dy = pydart.interp.interp_wghts(ilat, mlats)
            if height[n] != ihgt:
              ihgt = height[k,j,i]
              k1, k2, dz1, dz2, dz = pydart.interp.interp_wghts(ihgt, mhgts)

         if i1 < 0 or j1 < 0 or k1 < 0:  continue

         b[:] = missing

         if dbzOnly:  # DBZ
            for m, key in enumerate( ["dbz"]):
               q1     = dx1*fcst[key][k1,j1,i1] + dx2*fcst[key][k1,j1,i2]
               q2     = dx1*fcst[key][k1,j2,i1] + dx2*fcst[key][k1,j2,i2]
               vb     = (dy1*q1 + dy2*q2) / ( dx*dy )
               q1     = dx1*fcst[key][k2,j1,i1] + dx2*fcst[key][k2,j1,i2]
               q2     = dx1*fcst[key][k2,j2,i1] + dx2*fcst[key][k2,j2,i2]
               vt     = (dy1*q1 + dy2*q2) / ( dx*dy )
               b[m] = (dz1*vb + dz2*vt) / dz
            Hx_Z[k,j,i] = np.clip(b[0],_dbz_min,_dbz_max)


         else:  #DBZ/VR
            for m, key in enumerate( ["u", "v", "w", "dbz", "rho"]):
               #print(fcst[key].shape)
               #print('j1 ',j1)
               #print('i1 ',i1)
               #print('i2 ',i2)
               #print('mylocs ',mylocs.shape)
               #print('mxlocs ',mxlocs.shape)
               #print('WORK = ',fcst[key].shape)
               q1     = dx1*fcst[key][k1,j1,i1] + dx2*fcst[key][k1,j1,i2]
               q2     = dx1*fcst[key][k1,j2,i1] + dx2*fcst[key][k1,j2,i2]
               vb     = (dy1*q1 + dy2*q2) / ( dx*dy )
               q1     = dx1*fcst[key][k2,j1,i1] + dx2*fcst[key][k2,j1,i2]
               q2     = dx1*fcst[key][k2,j2,i1] + dx2*fcst[key][k2,j2,i2]
               vt     = (dy1*q1 + dy2*q2) / ( dx*dy )
               b[m] = (dz1*vb + dz2*vt) / dz
            # In CM1, dont have fall speed from microphysics, so use typical DBZ power law here
            refl   = 10.0**(0.1*np.clip(b[3],_dbz_min,_dbz_max))
            vfall  =  2.6 * refl**0.107 * (1.2/b[4])**0.4
            Hx_Z[k,j,i] = np.clip(b[3],_dbz_min,_dbz_max)  
            if Hx_Z[k,j,i] >= clear_dbz: # --- Remove Radial Velocity Observations from Regions of Clear Air
              Hx_vr[k,j,i] = b[0]*math.sin(azimuth[k,j,i])*math.cos(elv[k,j,i]) + b[1]*math.cos(azimuth[k,j,i])*math.cos(elv[k,j,i]) + (b[2]-vfall)*math.sin(elv[k,j,i])
            else:
               Hx_vr[k,j,i] = np.nan  #--- Don't consider Vr in regions of low Z

# END OBS_OP
  if dbzOnly:
     return Hx_Z
  else:
     return Hx_vr,Hx_Z  #,yloc,xloc

def theta_to_temp(pt,p):
   """
   Calculate temperature given potential temperature and pressure

   input: Both values can be an array or integer
      pt:   Potential Temperature (K
       p:   Air Pressure (Pa)

   return
       T:   Air Temperature (C)
   """

   R = 287.            # J[kg * K]^-1
   cp = 1004.          # J[kg * K]^-1
   R_over_cp = R / cp  # J[kg * K]^-1
   p_0 = 100000.
   T = pt / ( (p_0 / p) ** R_over_cp )

   return T - 273.15

def qv_to_spechum(qv) : #--- JDL Does Forward Operator Calculatoe Specific Humidity or Just qv
   """
   Calculate temperature given potential temperature and pressure

   input: Both values can be an array or integer
      qv:   water vapor mixing ratio (kg/kg)
   return
       w:   specific himidity
   """

   return qv/(1+qv)