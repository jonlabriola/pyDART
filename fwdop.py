#--- Code for Forward Operators Used to Calculate Observations
import interp
import numpy as np
import math

#ens_default_cmap = P.get_cmap('seismic')

_missing = -9900000000.
_dbz_min = 0.0
_dbz_max = 75.

_clevels = {'DBZ':[5,75,5], 'W':[-10,11,1], 'TH':[295.,325.,1.], 'THP':[-5,6,0.5], 'QV':[0.0,0.022,0.002], \
            'WZ': [-100.,120.,20.], 'U':[-20,22,2], 'V':[-20,22,2]}

def calcHx(fcst, kind, yloc, xloc, height, elev, azimuth,
           cartesian=False,rads=False,missing=None,clear_dbz=None):
  """
  The Radar Forward Operator Used to Calculate Z and Radial Wind
  Required Arguments
      1.) fcst    - The  model state variables                [nvars] [nk,ny,nz] - object contains vars and grid info
      2.) kind    - The kind of variables you work with       [nobs] - Observation Type (11 = Vr, 12 = dBZ)
      3.) yloc    - Observation Y-Location (m or Latitude)    [nobs]
      4.) xloc    - Observation X-Location (m or Longitude)   [nobs]   
      5.) height  - Observation Height above Surface (m)      [nobs]
      6.) elev    - Radar elevation angle                     [nobs]
      7.) azimuth - radar azimuth angle                       [nobs]


  Optional Arguments:
     cartesian    -  Are we working with a cartesian grid    [bool]
     rads         -  Are angle units radians?                [bool]
     missing      -  Value of a missing observation          [float]
     clear_dbz    -  What is defined as clear air (remove 
                     vr obs accordingly                      [float] 
    
  """
  nobs = np.size(kind) 

  if missing == None:
    missing = _missing

  Hx_vr = missing * np.ones([nobs])
  Hx_Z  = missing * np.ones([nobs])
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

  for n in np.arange(nobs):
    #if N.isnan(xloc[n]) == False: continue
    if np.isnan(xloc[n]): 
       Hx_vr[n] = np.nan
       Hx_Z[n] = np.nan
       continue

    if n%20000. == 0:
       print('Working on ob =',n)
  # need to acess lat lon alt info of ob loc from table and pass it to tlint
  # check to see if ob position is the same then less work

  #--- JDL If Statements are used to limit repeitition (the same i1,i2 ... values preserved)
    if cartesian: #--- Use x/y/z Coordinates
      if xloc[n] != ixloc:
         ixloc = xloc[n]
         i1, i2, dx1, dx2, dx = interp.interp_wghts(ixloc, mxlocs)
      if yloc[n] != iyloc:
         iyloc = yloc[n]
         j1, j2, dy1, dy2, dy = interp.interp_wghts(iyloc, mylocs)
      if height[n] != ihgt:
        ihgt = height[n]
        k1, k2, dz1, dz2, dz = interp.interp_wghts(ihgt, mhgts)

    else: #--- Use Lat/Lon Coordinates
      if lon[n] != ilon:
        ilon = lon[n]
        i1, i2, dx1, dx2, dx = interp.interp_wghts(ilon, mlons)
      if lat[n] != ilat:
        ilat = lat[n]
        j1, j2, dy1, dy2, dy = interp.interp_wghts(ilat, mlats)
      if height[n] != ihgt:
        ihgt = height[n]
        k1, k2, dz1, dz2, dz = interp.interp_wghts(ihgt, mhgts)

    if i1 < 0 or j1 < 0 or k1 < 0:  continue
    if rads:
       el = elev[n]
    else:
       el = np.deg2rad(elev[n])

  # The azimuth (in my PAR file) have north as 0 degreees, but I think this is what is needed
    if rads:
       az = azimuth[n]
    else:
       az = np.deg2rad(azimuth[n])

  # have to avoid using last member since it the mean

    b[:] = missing

    if kind[n] == 11:  # VR (AND DBZ Since Already Calculated)!
      #--- Bi-Linear Interpolation To Get Vr
      for m, key in enumerate( ["ua", "va", "wa", "dbz", "rho"]):
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

    # model's VR


      Hx_Z[n] = np.clip(b[3],_dbz_min,_dbz_max)  

 
      if Hx_Z[n] >= clear_dbz: # --- Remove Radial Velocity Observations from Regions of Clear Air
         Hx_vr[n] = b[0]*math.sin(az)*math.cos(el) + b[1]*math.cos(az)*math.cos(el) + (b[2]-vfall)*math.sin(el)
      else:
         Hx_vr[n] = np.nan  #--- Don't consider Vr in regions of low Z
    if kind[n] == 12:  # DBZ

      for m, key in enumerate( ["dbz"]):
        q1     = dx1*fcst[key][k1,j1,i1] + dx2*fcst[key][k1,j1,i2]
        q2     = dx1*fcst[key][k1,j2,i1] + dx2*fcst[key][k1,j2,i2]
        vb     = (dy1*q1 + dy2*q2) / ( dx*dy )
        q1     = dx1*fcst[key][k2,j1,i1] + dx2*fcst[key][k2,j1,i2]
        q2     = dx1*fcst[key][k2,j2,i1] + dx2*fcst[key][k2,j2,i2]
        vt     = (dy1*q1 + dy2*q2) / ( dx*dy )
        b[m] = (dz1*vb + dz2*vt) / dz

    # model's DBZ

      Hx_Z[n] = np.clip(b[0],_dbz_min,_dbz_max)
      
# END OBS_OP

  return Hx_vr,Hx_Z  #,yloc,xloc
