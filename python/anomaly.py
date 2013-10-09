import numpy as np
import netCDF4
import glob
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, addcyclic

import conn
import output_data
import plotPython

gref = 980. #big scope so accessible

def calc_dGravity_height(lat, z):
  #for lat in radians, z in meters,
  #return gravity adjustment for height, cm/s^2
  #http://glossary.ametsoc.org/wiki/Acceleration_of_gravity
  dg = -(3.085462e-4+2.27e-7*np.cos(2*lat))*z
  return dg

def calc_geopotential(lat, z):
  cos2phi = np.cos(2*lat)
  g0 = 980.6160*(1-0.0026373*cos2phi+0.0000059*cos2phi*cos2phi)
  #gref = 980.
  
  #discretize height and integrate
  nSteps = 31
  zRange, dz = np.linspace(0., z, num=nSteps, endpoint=True, retstep=True)
  
  #else, perform integration
  sum_gdz = 0.
  for i in xrange(nSteps-1):
    zMid = .5*(zRange[i]+zRange[i+1])
    gLayer = g0+calc_dGravity_height(lat, zMid)
    sum_gdz = sum_gdz + gLayer*dz
  return sum_gdz/gref

def mpas_hgt2geopotential(nCells, htCells, latCells):
  geop = np.empty(nCells, dtype=float)
  for iCell in xrange(nCells):
    geop[iCell] = calc_geopotential(latCells[iCell], htCells[iCell])
    
  return geop
  
def makeMap_ll2mpas(latSrc, lonSrc, latDest, lonDest, nEdgesOnCell, cellsOnCell):
  #map the lat lon points in the src to the closest in the dest.
  #return the map
  
  nlat = len(latSrc)
  nlon = len(lonSrc)
  llMap = np.empty((nlat,nlon),dtype=int)
  
  pt_ll = np.empty(2,dtype=float)
  cellId = 0
  for iLat in xrange(nlat):
    for iLon in xrange(nlon):
      pt_ll[0] = latSrc[iLat]
      pt_ll[1] = lonSrc[iLon]
      
      #cellId=conn.findOwner_horizNbrs_latLon(pt_ll, cellId, latDest, lonDest, nEdgesOnCell, cellsOnCell)
      cellId=conn.findOwner_horizNbrs_latLon_storeDistance(pt_ll, cellId, latDest, lonDest, nEdgesOnCell, cellsOnCell)
      llMap[iLat,iLon] = cellId
      
  return llMap

def calc_diff_field(field_ll, field_mpas, llMap, nlat, nlon):
  #mpas-mean
  diffField = np.empty((nlat, nlon), dtype=float)
  
  for iLat in xrange(nlat):
    for iLon in xrange(nlon):
      iCell = llMap[iLat, iLon]
      diffField[iLat, iLon] = field_mpas[iCell]-field_ll[iLat, iLon]
  return diffField

def test_run_hgt():
  #long term mean
  fname_mean = '/data01/forONR/hgt.mon.1981-2010.ltm.nc'
  lev_500 = 5 #level of 500hPa
  month = 8 #august
  
  dataMean = netCDF4.Dataset(fname_mean,'r')
  hgtMean = dataMean.variables['hgt'][month,lev_500,:,:] #lat,lon
  #hgtMean += dataMean.variables['hgt'][month-1,lev_500,:,:]
  #hgtMean += dataMean.variables['hgt'][month-2,lev_500,:,:]
  #hgtMean /= 3.
  latMean = dataMean.variables['lat'][:]*np.pi/180. #stored in degrees!
  lonMean = dataMean.variables['lon'][:]*np.pi/180.
  dataMean.close()
  
  nlat = len(latMean)
  nlon = len(lonMean)
  var = np.zeros((nlat,nlon), dtype=float)
  
  #mpas files
  #fnames = []
  fnames = sorted(glob.glob('/arctic1/nick/cases/vduda/x4/x4.t.output.2006-08-*'))
  fnames.extend(sorted(glob.glob('/arctic1/nick/cases/vduda/x4.t.output.2006-08-15*')))
  #fnames = sorted(glob.glob('/arctic1/nick/cases/vduda/x4/*.output.2006-08*'))
  counter = 0
  for iFile, fpath in enumerate(fnames):
    data = output_data.open_netcdf_data(fpath)
    nEdgesOnCell = data.variables['nEdgesOnCell'][:];
    cellsOnCell = data.variables['cellsOnCell'][:]-1;
    latCell = data.variables['latCell'][:];
    lonCell = data.variables['lonCell'][:];
    nTimes = len(data.dimensions['Time'])
    #nTimes = 1 #3
    
    nCells = len(latCell)
    geop_mpas = np.zeros(nCells,dtype=float)
    
    for iTime in xrange(nTimes):
      #get average value for mesh over these times so can even anomaly different meshes together
      output_data.print_xtime(data, iTime)
      hgt_mpas = data.variables['height_500hPa'][iTime,:]
      #average here to avoid summing to large number over many times
      geop_mpas += mpas_hgt2geopotential(nCells, hgt_mpas, latCell)/nTimes
      
    llMap = makeMap_ll2mpas(latMean, lonMean, latCell, lonCell, nEdgesOnCell, cellsOnCell)
    var += calc_diff_field(hgtMean, geop_mpas, llMap, nlat, nlon)
    #print var
    counter = counter+1
    
    data.close()
    
  var /= counter
  #var = hgtMean
  # add wrap-around point in longitude.
  var, lonMean = addcyclic(var, lonMean)
  np.savez('tmp.dat',var,latMean,lonMean)
  lats,lons = np.meshgrid(latMean, lonMean)
  lats *= 180./np.pi; lons *= 180./np.pi
  plot_anomaly_ll(lats, lons, np.transpose(var))

def test_plot():
  fname_mean = '/data01/forONR/hgt.mon.1981-2010.ltm.nc'  
  dataMean = netCDF4.Dataset(fname_mean,'r')
  latMean = dataMean.variables['lat'][:]
  lonMean = dataMean.variables['lon'][:]
  dataMean.close()
  
  nlat = len(latMean)
  nlon = len(lonMean)
  var = np.zeros((nlat,nlon), dtype=float)
  for iLat in xrange(nlat):
    for iLon in xrange(nlon):
      var[iLat, iLon] = lonMean[iLon] #latMean[iLat]
  
  var, lonMean = addcyclic(var, lonMean)
  lats,lons = np.meshgrid(latMean, lonMean)
  plot_anomaly_ll(lats, lons, np.transpose(var))
  
def compare_mpas_pHgts():
  f = '/arctic1/nick/cases/vduda/x4/x4.t.output.2006-08-01_00.00.00.nc'
  data = netCDF4.Dataset(f,'r')
  nLevels = len(data.dimensions['nVertLevels'])
  nCells = len(data.dimensions['nCells'])
  iTime = 0
  hgt_mpas = data.variables['height_500hPa'][iTime,:]
  
  zgrid = data.variables['zgrid'][:]
  zMid = .5*(zgrid[:,0:-1]+zgrid[:,1:])
  p = data.variables['pressure_p'][iTime,:,:]+data.variables['pressure_base'][iTime,:,:]
  
  hgt_calc = np.empty(nCells,dtype=float)
  
  pLev = 50000.
  for iCell in xrange(nCells):
    (l,dl) = output_data.calcIndexOfValue(pLev, p[iCell,:], nLevels)
    z = output_data.calcValueOfIndex(l,dl,zMid[iCell,:])
    hgt_calc[iCell] = z
  
  diff = hgt_calc-hgt_mpas
  return diff

def compare_runs():
  f6 = '/arctic1/nick/cases/vduda/x4/x4.t.output.2006-08-01_00.00.00.nc'
  f7 = '/arctic1/nick/cases/2007/x4.t.output.2007-08-01_00.00.00.nc'
  
  data6 = netCDF4.Dataset(f6,'r')
  data7 = netCDF4.Dataset(f7,'r')
  nTimes6 = len(data6.dimensions['Time'])
  nTimes7 = len(data7.dimensions['Time'])
  nTimes = min((nTimes6, nTimes7))
  iTime = 0
  hDiff = data6.variables['height_500hPa'][iTime,:]-data7.variables['height_500hPa'][iTime,:]
  for iTime in xrange(nTimes):
    hDiff += data6.variables['height_500hPa'][iTime,:]-data7.variables['height_500hPa'][iTime,:]
  hDiff /= nTimes
  
  
  latCell = data6.variables['latCell'][:];
  lonCell = data6.variables['lonCell'][:];
  plotPython.plotSfc(latCell, lonCell, hDiff)
  
  data6.close()
  data7.close()

def compare_runs_vorticity():
  #time series of vorticity North of 60 North
  latThresh = 80.*np.pi/180.

  f6 = '/arctic1/nick/cases/vduda/x4/x4.t.output.2006-08-01_00.00.00.nc'
  f7 = '/arctic1/nick/cases/2007/x4.t.output.2007-08-01_00.00.00.nc'

  data6 = netCDF4.Dataset(f6,'r')
  data7 = netCDF4.Dataset(f7,'r')
  
  latVertex = data6.variables['latVertex'][:]
  vertN = latVertex>latThresh
  areaTri = data6.variables['areaTriangle'][:]
  nTimes6 = len(data6.dimensions['Time'])
  nTimes7 = len(data7.dimensions['Time'])
  nTimes = min((nTimes6, nTimes7))
  iTime = 0
  vort6 = []
  for iTime in xrange(nTimes):
    vort = data6.variables['vorticity_500hPa'][iTime,:]
    netVort = netVorticity_cells(vertN, vort, areaTri)
    vort6.append(netVort)
  vort7 = []
  for iTime in xrange(nTimes):
    vort = data7.variables['vorticity_500hPa'][iTime,:]
    netVort = netVorticity_cells(vertN, vort, areaTri)
    vort7.append(netVort)

  tRange = np.array(range(nTimes))*6 #in hours
  print len(tRange), len(vort6), len(vort7)
  plt.figure()
  plt.plot(tRange, vort6, 'b.-', label='2006')
  plt.plot(tRange, vort7, 'r--', label='2007')
  plt.ylabel("Vorticity, 1/s")
  plt.ylabel("Forecast hour, hour")
  plt.legend()
  plt.show()
  data6.close()
  data7.close()

def netVorticity_cells(verts, vortVertex, areaTri):
  netCirc = np.dot(vortVertex[verts],areaTri[verts])
  netArea = np.sum(areaTri[verts])
  return netCirc/netArea

def plot_anomaly_ll(lat, lon, var):
  #Input lat and lon in degrees!!!
  plt.figure()
  
  #m = Basemap(projection='ortho',lon_0=100,lat_0=60, resolution='l')
  m = Basemap(projection='ortho',lon_0=0,lat_0=90, resolution='l')
  x,y = m(lon, lat)
  
  m.drawcoastlines()
  m.drawmapboundary()
  
  #'''
  #set bounds on var
  maxVal = 100.; minVal = -100.
  var[var>maxVal] = maxVal
  var[var<minVal] = minVal
  #'''
  
  m.pcolor(x,y,var,shading='flat',edgecolors='none',cmap=plt.cm.jet) #cmap=plt.cm.hot_r) #,vmin=100,vmax=1000)
  #m.contour(x,y,var,10,tri=True,shading='flat',edgecolors='none',cmap=plt.cm.jet)
  #m.contourf(x,y,var,10,tri=True,shading='flat',edgecolors='none',cmap=plt.cm.jet)
  #cbar = m.colorbar(plt1,location='bottom',pad="5%")
  #cbar.set_label('\Delta m')
  plt.colorbar()
  plt.show()

if __name__=='__main__':
  test_run_hgt()
  #compare_runs_vorticity()
