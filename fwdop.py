#--- Code for Forward Operators Used to Calculate Observations
import location
import numpy as N
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
  nobs = N.size(kind) 

  if missing == None:
    missing = _missing

  Hx_vr = missing * N.ones([nobs])
  Hx_Z  = missing * N.ones([nobs])
  if cartesian: #--- x/y/z coordinates
     yloc = yloc
     xloc = xloc

     ixloc = 0.0
     iyloc = 0.0
     ihgt = 0.0

     mxlocs = fcst['xh']
     mylocs = fcst['yh']
     mhgts = fcst['zh']

     print(' calcHx:  OBS   LAT MIN/MAX:  ', yloc.min(), yloc.max())
     print(' calcHx:  MODEL LAT MIN/MAX:  ', mylocs.min(), mylocs.max())
     print(' calcHx:  OBS   LON MIN/MAX:  ', xloc.min(), xloc.max())
     print(' calcHx:  MODEL LON MIN/MAX:  ', mxlocs.min(), mxlocs.max())

     print(' calcHx:  OBS   HGT MIN/MAX:  ', height.min(), height.max())
     print(' calcHx:  MODEL HGT MIN/MAX:  ', mhgts.min(), mhgts.max())

  else: #--- Lat/Lon Coordinates
     lat = yloc
     lon = xloc

     ilon = 0.0
     ilat = 0.0
     iihgt = 0.0

     mlons = fcst['LON']  #--- JDL Will Need to UPDATE 
     mlats = fcst['LAT']  #--- JDL Will Need To UPDATE
     mhgts = fcst['zh']

     print(' calcHx:  OBS   LAT MIN/MAX:  ', lat.min(), lat.max())
     print(' calcHx:  MODEL LAT MIN/MAX:  ', mlats.min(), mlats.max())
     print(' calcHx:  OBS   LON MIN/MAX:  ', lon.min(), lon.max())
     print(' calcHx:  MODEL LON MIN/MAX:  ', mlons.min(), mlons.max())

     print(' calcHx:  OBS   HGT MIN/MAX:  ', height.min(), height.max())
     print(' calcHx:  MODEL HGT MIN/MAX:  ', mhgts.min(), mhgts.max())

  b = N.zeros([5])

  for n in N.arange(nobs):
    #if N.isnan(xloc[n]) == False: continue
    if N.isnan(xloc[n]): 
       Hx_vr[n] = N.nan
       Hx_Z[n] = N.nan
       continue

    if n%20000. == 0:
       print('Working on ob =',n)
  # need to acess lat lon alt info of ob loc from table and pass it to tlint
  # check to see if ob position is the same then less work

  #--- JDL If Statements are used to limit repeitition (the same i1,i2 ... values preserved)
    if cartesian: #--- Use x/y/z Coordinates
      if xloc[n] != ixloc:
         ixloc = xloc[n]
         i1, i2, dx1, dx2, dx = location.interp_wghts(ixloc, mxlocs)
      if yloc[n] != iyloc:
         iyloc = yloc[n]
         j1, j2, dy1, dy2, dy = location.interp_wghts(iyloc, mylocs)
      if height[n] != ihgt:
        ihgt = height[n]
        k1, k2, dz1, dz2, dz = location.interp_wghts(ihgt, mhgts)

    else: #--- Use Lat/Lon Coordinates
      if lon[n] != ilon:
        ilon = lon[n]
        i1, i2, dx1, dx2, dx = location.interp_wghts(ilon, mlons)
      if lat[n] != ilat:
        ilat = lat[n]
        j1, j2, dy1, dy2, dy = location.interp_wghts(ilat, mlats)
      if height[n] != ihgt:
        ihgt = height[n]
        k1, k2, dz1, dz2, dz = location.interp_wghts(ihgt, mhgts)

    if i1 < 0 or j1 < 0 or k1 < 0:  continue
    if rads:
       el = elev[n]
    else:
       el = N.deg2rad(elev[n])

  # The azimuth (in my PAR file) have north as 0 degreees, but I think this is what is needed
    if rads:
       az = azimuth[n]
    else:
       az = N.deg2rad(azimuth[n])

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

      refl   = 10.0**(0.1*N.clip(b[3],_dbz_min,_dbz_max))
      vfall  =  2.6 * refl**0.107 * (1.2/b[4])**0.4

    # model's VR


      Hx_Z[n] = N.clip(b[3],_dbz_min,_dbz_max)  

 
      if Hx_Z[n] >= clear_dbz: # --- Remove Radial Velocity Observations from Regions of Clear Air
         Hx_vr[n] = b[0]*math.sin(az)*math.cos(el) + b[1]*math.cos(az)*math.cos(el) + (b[2]-vfall)*math.sin(el)
      else:
         Hx_vr[n] = N.nan  #--- Don't consider Vr in regions of low Z
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

      Hx_Z[n] = N.clip(b[0],_dbz_min,_dbz_max)
      
# END OBS_OP

# Remove missing Hx's and strip the input data of those points as well...
  #Hxf[:]=Hx[:]
  return Hx_vr,Hx_Z  #,yloc,xloc
  #idx = N.where(Hx[:] != missing)
  #nobs= N.size(idx)
  #Hxf = N.zeros([nobs])
  #Hxf[:]= Hx[idx]

  #if nobs > 0:
  #  if cartesian:
  #     return idx, Hxf, kind[idx], yloc[idx], xloc[idx], height[idx], elev[idx], azimuth[idx]
  #  else:
  #     return idx, Hxf, kind[idx], lat[idx], lon[idx], height[idx], elev[idx], azimuth[idx]
  #else:
  #  return None, 0, 0, 0, 0, 0, 0, 0
