#matplotlib basemap examples (all meteorology ones):
# http://matplotlib.org/basemap/users/examples.html

#python libs
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from mpl_toolkits.basemap import Basemap
import datetime as dt

#my libs
import output_data
import vars
import sys; sys.path.append("/data01/tracks/");
import cfsrStats
import mpasStats

def plot_dateRange_example():
  base = dt.datetime.today()
  numdays = 10
  dateList = [ base + dt.timedelta(days=x) for x in range(0,numdays) ]
  
  vals = np.arange(len(dateList))

  plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y/%m/%d/%H'))
  #matplotlib.dates.DayLocator(bymonthday=None, interval=1, tz=None)
  plt.gca().xaxis.set_major_locator(mdates.DayLocator())
  plt.plot(dateList,vals)
  plt.plot(dateList,-vals)
  plt.gcf().autofmt_xdate()
  
  plt.show()

def plotPoints_example():
  fpath = '/arctic1/nick/cases/cfsr/output/x4.cfsr.output.2006-08-01_00.00.00.nc'
  data = output_data.open_netcdf_data(fpath)
  latCell = data.variables['latCell'][:];
  lonCell = data.variables['lonCell'][:];
  data.close()
  
  plt.figure()
  map = Basemap(projection='ortho',lon_0=0,lat_0=90, resolution='l')
  '''
  cellInds_cfsr = [56355, 28662, 11263, 131726, 84718, 84714, 46582, 105453, 114078, 
  25622, 37141, 6428, 50655, 160926, 144268, 156040, 98770, 42664, 137842, 24214, 
  151512, 2568, 96729, 11246, 131490, 19940, 163620, 119379, 136917, 104557, 36439, 
  150649, 21865, 30589, 119446, 71817, 10022, 131487, 84705, 156231, 84691, 5142, 
  151475, 79538, 28657, 79527, 162909, 151551, 160897, 119355, 151474, 114095, 
  151473, 41958, 131519, 131551, 131556, 3564]
  
  cellInds_t = [56355, 14136, 62441, 105460, 62444, 131720, 84715, 84718, 46589, 58118, 103293, 
  119276, 119271, 38013, 131725, 20012, 131711, 96797, 1630, 20010, 79513, 79515, 
  151526, 131702, 41980, 71822, 144522, 137231, 14763]
  
  cellInds_kf = [56355, 14136, 62441, 105460, 62444, 50736, 84714, 84713, 46589, 
  58118, 10538, 119277, 15864, 84711, 131725, 20012, 50734, 50738, 131727, 119266, 
  38009, 71828, 10022, 33361, 8554, 137254, 66155, 4203, 137296]
  '''
  cellInds_cfsr = cfsrStats.x4_cfsr_july_cells
  
  lonp = lonCell[cellInds_cfsr]*180./np.pi
  latp = latCell[cellInds_cfsr]*180./np.pi
  xcfsr,ycfsr = map(lonp, latp)
  
  #cells_guess = [40419,26537,10552,1634,7338,20030,437,2896,11920,22053,19965,19961,28617]
  cells_guess = mpasStats.x4_t_0724_cells
  cells_f500 = [138657, 149403, 38051, 28687, 131782, 131766, 71977, 58328, 22007, 11933, 14811, 72391, 103406]
  cells_f1000 = [138657, 149403, 38051, 28687, 131782, 131766, 71977, 58328, 33409, 149420, 5194, 93157, 16019]

  lonp = lonCell[cells_guess]*180./np.pi
  latp = latCell[cells_guess]*180./np.pi
  xcfsr_g,ycfsr_g = map(lonp, latp)
  
  lonp = lonCell[cells_f500]*180./np.pi
  latp = latCell[cells_f500]*180./np.pi
  xcfsr_f500,ycfsr_f500 = map(lonp, latp)
  
  lonp = lonCell[cells_f1000]*180./np.pi
  latp = latCell[cells_f1000]*180./np.pi
  xcfsr_f1000,ycfsr_f1000 = map(lonp, latp)

  '''
  lonp = lonCell[cellInds_t]*180./np.pi
  latp = latCell[cellInds_t]*180./np.pi
  xtiedtke,ytiedtke = map(lonp, latp)
  
  lonp = lonCell[cellInds_kf]*180./np.pi
  latp = latCell[cellInds_kf]*180./np.pi
  xkf,ykf = map(lonp, latp)
  '''
  map.drawcoastlines()
  map.drawmapboundary()
  #map.scatter(x,y, picker=5)
  #map.scatter(x,y)
  map.plot(xcfsr,ycfsr,'b*--', label='back')
  map.plot(xcfsr_g,ycfsr_g, 'go-', label='07-24_12') # label='user')
  map.plot(xcfsr_f500,ycfsr_f500,'gs:', label='f,r500km')
  map.plot(xcfsr_f1000,ycfsr_f1000,'r*-.', label='f,r1000km')
  plt.legend()
  '''
  legend = plt.legend('back', 'user', 'forward')
  for label in legend.get_texts():
    label.set_fontsize('medium')
  '''
  #map.plot(xcfsr,ycfsr,'b*--', xtiedtke,ytiedtke, 'go-.', xkf,ykf, 'rs:') #, markersize=10)
  #plt.legend('cfsr', 'tiedtke', 'kf')
  plt.show()

def plotError_example():
  import netCDF4
  import datetime as dt
  #fFile = '/arctic1/mduda/60km/output.163842.2006-07-13_00.00.00.nc' #forecast
  fFile ='/arctic1/mduda/60km/restart.163842.2006-07-13_00.00.00.nc'
  anaFile = '/arctic1/nick/cases/163842/r2614/sfc_update.nc' #"analysis"
  
  fdata = netCDF4.Dataset(fFile,'r')
  adata = netCDF4.Dataset(anaFile,'r')
  
  #these files can start at different times but assume both are 6 hours.
  #we get these values from the attributes: print data.ncattrs()
  fStart = '2006-07-13_00:00:00' #fdata.config_start_time
  aStart = '2006-07-08_00:00:00' #adata.config_start_time
  tfStart = dt.datetime.strptime(fStart, '%Y-%m-%d_%H:%M:%S')
  taStart = dt.datetime.strptime(aStart, '%Y-%m-%d_%H:%M:%S')
  
  #aIndex = fIndex+f2aInd
  f2aInd = int((tfStart-taStart).total_seconds()/(60.*60.*6)) #assuming 6 hour steps per index
  print "For f:"+fStart+" and a:"+aStart+", dIndex="+str(f2aInd)+'\n'
  
  #field differences
  fIndex=0; aIndex=fIndex+f2aInd
  key = 'sst'
  avar = adata.variables[key][aIndex,:]
  fvar = fdata.variables[key][fIndex,:]
  diff = fvar-avar
  
  lat = fdata.variables['latCell'][:]*180./np.pi; #in degrees
  lon = fdata.variables['lonCell'][:]*180./np.pi;
  plotSfc(lat, lon, diff)
  
def example_structured_interp():
  #see http://www.scipy.org/Cookbook/Matplotlib/Gridding_irregularly_spaced_data
  #rather than work in unstructured space, we can also move everything to lat/lon world
  from scipy.interpolate import griddata #remember that distance is calculated in euclidean space!!!
  
  #since distance is calculated as euclidean, is nearest nbr interpolation "less wrong"?
  pass

def example_structured_weight():
  #we can also do a cressman or histogram procedure about each regular grid point
  #with the model output as observations.
  pass

def example_sfcUpdate_horiz(time,var):
  #The surfaceUpdate.nc files don't have all of the same fields/dimensions,...
  #so we get errors when attempting to access them
  
  import netCDF4
  
  ncfname = '/arctic1/nick/cases/163842/r2614/sfc_update.nc'
  ncfnameG = '/arctic1/nick/cases/163842/r2614/output.163842.2006-07-08_00.00.00.nc' #used for grid since, say, sfc_update.nc doesn't have this info
  dataG = netCDF4.Dataset(ncfnameG,'r')
  data = netCDF4.Dataset(ncfname,'r')
  nCells = len(data.dimensions['nCells'])
  
  lat = dataG.variables['latCell'][:]*180./np.pi; #in degrees
  lon = dataG.variables['lonCell'][:]*180./np.pi;
  
  level = 0; #time=20; var='sst';
  var = data.variables[var][time,:];# minVar = np.amin(var); maxVar = np.amax(var);
  
  #map = Basemap(projection='ortho',lon_0=-105,lat_0=40, resolution='l')
  map = Basemap(projection='ortho',lon_0=-100,lat_0=60, resolution='l')
  x,y = map(lon, lat)
  
  fig1 = plt.figure(1)
  map.drawcoastlines()
  map.drawmapboundary()
  map.pcolor(x,y,var,tri=True,shading='flat',edgecolors='none',cmap=plt.cm.jet) #cmap=plt.cm.hot_r) #,vmin=100,vmax=1000)
  #map.contour(x,y,var,10,tri=True,shading='flat',edgecolors='none',cmap=plt.cm.jet)
  #map.contourf(x,y,var,10,tri=True,shading='flat',edgecolors='none',cmap=plt.cm.jet)
  plt.colorbar()
  
  if(1):
    map = Basemap(projection='ortho',lon_0=100,lat_0=-60, resolution='l')
    x,y = map(lon, lat)

    fig2 = plt.figure(2)
    map.drawcoastlines()
    map.drawmapboundary()
    map.pcolor(x,y,var,tri=True,shading='flat',edgecolors='none',cmap=plt.cm.jet) #cmap=plt.cm.hot_r) #,vmin=100,vmax=1000)
    #map.contour(x,y,var,10,tri=True,shading='flat',edgecolors='none',cmap=plt.cm.jet)
    #map.contourf(x,y,var,10,tri=True,shading='flat',edgecolors='none',cmap=plt.cm.jet)
    plt.colorbar()
  
  #plotName = '/home/nickszap/Desktop/sst'+str(time)+'.png'
  #fig1.savefig(plotName)
  
  plt.show()

def plotSfc(lat, lon, var):  
  #map = Basemap(projection='ortho',lon_0=100,lat_0=60, resolution='l')
  map = Basemap(projection='ortho',lon_0=0,lat_0=90, resolution='l')
  x,y = map(lon, lat)
  
  fig1 = plt.figure(1)
  map.drawcoastlines()
  map.drawmapboundary()
  map.pcolor(x,y,var,tri=True,shading='flat',edgecolors='none',cmap=plt.cm.jet) #cmap=plt.cm.hot_r) #,vmin=100,vmax=1000)
  #map.contour(x,y,var,10,tri=True,shading='flat',edgecolors='none',cmap=plt.cm.jet)
  #map.contourf(x,y,var,10,tri=True,shading='flat',edgecolors='none',cmap=plt.cm.jet)
  plt.colorbar()
  plt.show()

def plot_2dll(lat, lon, var):  
  #map = Basemap(projection='ortho',lon_0=100,lat_0=60, resolution='l')
  map = Basemap(projection='ortho',lon_0=0,lat_0=90, resolution='l')
  x,y = map(lon, lat)
  
  fig1 = plt.figure(1)
  map.drawcoastlines()
  map.drawmapboundary()
  map.pcolor(x,y,var,shading='flat',edgecolors='none',cmap=plt.cm.jet) #cmap=plt.cm.hot_r) #,vmin=100,vmax=1000)
  #map.contour(x,y,var,10,tri=True,shading='flat',edgecolors='none',cmap=plt.cm.jet)
  #map.contourf(x,y,var,10,tri=True,shading='flat',edgecolors='none',cmap=plt.cm.jet)
  plt.colorbar()
  plt.show()

def example_mpas_horiz_tri():
  #plot triangulation on a map projection
  
  #ncfname = '/home/nickszap/research/mpas/output.2010-10-23_00:00:00.nc' #input file
  ncfname = '/arctic1/nick/cases/163842/r2614/output.163842.2006-07-08_00.00.00.nc'
  data = output_data.open_netcdf_data(ncfname)
  nCells = len(data.dimensions['nCells'])
  nVertices = len(data.dimensions['nVertices'])
  nLevels = len(data.dimensions['nVertLevels'])
  
  lat = data.variables['latCell'][:]*180./np.pi; #in degrees
  lon = data.variables['lonCell'][:]*180./np.pi;
  
  level = 0; time=15;
  var = data.variables['theta'][time,:,level];# minVar = np.amin(var); maxVar = np.amax(var);
  
  ''' #we can triangulate the map projection coords
  #although the needed connectivity information is in the mesh,
  #I don't feel like accumulating the edges into triangles so we'll compute a convex hull
  #that we "know" will be the same as the dual of the Voronoi
  xc = data.variables['xCell'][:]; yc = data.variables['yCell'][:]; zc = data.variables['zCell'][:];
  coords = np.array([xc,yc,zc]).transpose()
  triang = Delaunay(coords)
  '''
  
  '''
  #create the triangulation data structure since I've had issues with scipy's delaunay and such.
  #We can't triangulate all of the projected x,y since opposing sides of the earth will overlay each other
  tris = recoverTriangles(data, nCells); print "Triangulation finished\n"  
  
  nContours = 50
  #matplotlib.pyplot.tricontour
  #plt.tripcolor(x, y, tris, facecolors=zfaces, edgecolors='k')
  plt.tripcolor(lon, lat, tris, var)
  plt.colorbar()
  plt.show()
  '''
  
  #map = Basemap(projection='ortho',lon_0=-105,lat_0=40, resolution='l')
  map = Basemap(projection='ortho',lon_0=-100,lat_0=60, resolution='l')
  x,y = map(lon, lat)
  
  plt.figure(1)
  map.drawcoastlines()
  map.drawmapboundary()
  map.pcolor(x,y,var,tri=True,shading='flat',edgecolors='none',cmap=plt.cm.jet) #cmap=plt.cm.hot_r) #,vmin=100,vmax=1000)
  #map.contour(x,y,var,10,tri=True,shading='flat',edgecolors='none',cmap=plt.cm.jet)
  #map.contourf(x,y,var,10,tri=True,shading='flat',edgecolors='none',cmap=plt.cm.jet)
  plt.colorbar()
  
  if(0):
    plt.figure(3)
    map.drawcoastlines()
    map.drawmapboundary()
    map.contour(x,y,var,10,tri=True,shading='flat',edgecolors='none',cmap=plt.cm.jet)
    plt.colorbar()
  
  plt.show()
  '''
  map = Basemap(projection='ortho',lat_0=45,lon_0=-100,resolution='l')
  x, y = map(lon,lat)
  
  map.drawcoastlines(linewidth=0.25)
  map.drawcountries(linewidth=0.25)
  #cs = map.contour(lon,lat, 50, nContours, latlon=True, tri=True)
  #cs = map.contour(x,y, var, nContours, latlon=False, tri=True, triangulation=triang)
  #cs = map.pcolor(x,y, var, tri=True,triangulation=triang)
  
  plt.colorbar(cs)
  plt.title('Test plot')
  plt.show()
  '''
  
  data.close()

def example_mpas_vertical():
  #plot values of a radial column 
  pass
  
def recoverTriangles(data, nCells): #just use cellsOnVertex!!!
  #The cell centers form vertices of the triangulation.
  #From connectivity information in the mesh, we know edges between cell centers (cellsOnCell).
  #We need to construct elements out of these edges by:
  #edges->n2n ->cells are all triangles so n0 in n2n of nbr
  print "Why not use cellsOnVertex?\n"
  '''
  n2n = [[] for _ in xrange(nCells)]
  for i in xrange(nCells):
    #store cell nbrs
    n2n[i].extend(vars.getHorizNbrs(data, i).tolist()) #append would plop in the list as an entry
  #n2n finished (wow, that's easy!). probably can use a list comprehension too
  '''
  #1 liner? yay for list comprehensions
  n2n = [vars.getHorizNbrs(data, i).tolist() for i in xrange(nCells)]
  
  #We know # of triangles since each polygon vertex is center of triangle
  nTri = len(data.dimensions['nVertices'])
  tri = np.empty((nTri,3),np.dtype('int32'))
  
  triInd=0
  for i in xrange(nCells):
    for nbr in n2n[i]:
      #see if a nbr of this nbr goes back to original.
      #if it does, there will be 2 (no boundaries since we cover the earth)
      if (nbr<i): #we'll have already made this triangle
        continue
      for nbr2 in n2n[nbr]:
        if (nbr2<i or nbr2<nbr):
          continue
        if (i in n2n[nbr2]):
          #we've got a triangle
          tri[triInd,0]=i; tri[triInd,1]=nbr; tri[triInd,2]=nbr2;
          #print tri[triInd,:]
          triInd=triInd+1;
  #
  #some sanity checks
  if (triInd!=nTri):
    print "Uhoh. We didn't end up with the right number of triangles since %s is not %s\n"%(triInd,nTri)
  
  #do we need to wind the connectivities consistenly clockwise?
  return tri

def windTriangles(triInds, nTri, x,y,z):
  #winding should be such that normal points radially out, ie normal.radius>0
  for t in xrange(nTri):
    #normal is uxv
    n0 = triInds[t,0]; n1=triInds[t,1]; n2 = triInds[t,2];
    u = np.array([x[n1]-x[n0],y[n1]-y[n0],z[n1]-z[n0]])
    v = np.array([x[n2]-x[n0],y[n2]-y[n0],z[n2]-z[n0]])
    n = np.cross(u,v)
    r = np.array([x[n0], y[n0], z[n0]])
    if (np.dot(n,r)<0): #rewind
      triInds[t,1] = n2; triInds[t,2] = n1;
  return triInds

def form_triangulationStructure(data, nCells):
  #for now, it's simpler to have the plotting package compute triangulations.
  #however, apparently we can also pass in the triangulation data structure, but
  #i haven't managed to get it to work
  pass

if __name__=='__main__':
  plotPoints_example()