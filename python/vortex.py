# -*- coding: utf-8 -*-

# imports
import numpy as np
import glob

import output_data
import conn
import compare

import sys; sys.path.append("/home/nickszap/Dropbox/pythonScripts/mem/");
import varsM

def run_vortexLoc():
  #fnames = cfsrFiles2()
  #fnames = ['/arctic1/nick/cases/v1.0/x4/2week/tiedtke/x4.tiedtke.output.2006-07-24_12.00.00.nc',
  #          '/arctic1/nick/cases/v1.0/x4/2week/tiedtke/x4.tiedtke.output.2006-07-31_12.00.00.nc']
  #fnames = ['/arctic1/nick/cases/v1.0/x4/longer/x4.kf.output.2006-08-15_00.00.00.nc']
  #fnames = sorted(glob.glob('/arctic1/nick/cases/cfsr/output/x4.cfsr.output.2006-07-*'))[-1::-1] #track backwards
  #fnames = ['/arctic1/nick/cases/cfsr/output/x4.cfsr.output.2006-08-15_00.00.00.nc']
  fnames = sorted(glob.glob('/arctic1/nick/cases/vduda/x4/x4.t.output.2006-07-*'))
  #fnames = ['/arctic1/nick/cases/vduda/x7/x7.kf.output.2006-08-07_18.00.00.nc']
  #fnames = sorted(glob.glob('/arctic1/nick/cases/vduda/x4.t.output.2006-09-*'))
  
  pt_ll = np.empty(2,dtype='float')
  #pt_ll[0] = 90.*np.pi/180.; pt_ll[1] = 0.0
  pt_ll[0] = 87.05*np.pi/180.; pt_ll[1] = 18.86*np.pi/180. #for 20060801 start time
  #pt_ll[0] = 85.38*np.pi/180.; pt_ll[1] = 352.97*np.pi/180. #for 2006080312 start time
  #pt_ll[0] = 86.1*np.pi/180.; pt_ll[1] = 137.*np.pi/180. #for 2006080718
  #pt_ll[0] = 82.5*np.pi/180.; pt_ll[1] = 41.95*np.pi/180. # for 2006081500
  #pt_ll[0] = 71.55*np.pi/180.; pt_ll[1] = 153.1*np.pi/180. #for 2006072412
  #pt_ll[0] = 81.05*np.pi/180.; pt_ll[1] = 19.4*np.pi/180. #for 2006091900
  
  #pt_ll[0] = 1.41; pt_ll[1] = 0.27
  #give initial guess of cClose. guess is just seed to search, so can choose anything.
  cClose = 71791 #131735 #x4 mesh for 20060801 start time
  #cClose = 102439 #x7 mesh for 20060801 start time
  c0IsMin = True;
  
  #load in mesh data
  fpath = fnames[0]
  data = output_data.open_netcdf_data(fpath)
  nEdgesOnCell = data.variables['nEdgesOnCell'][:];
  cellsOnCell = data.variables['cellsOnCell'][:]-1;
  latCell = data.variables['latCell'][:];
  lonCell = data.variables['lonCell'][:];
  data.close()
  
  cellList = []
  minValList = []
  for iFile, fpath in enumerate(fnames):
    data = output_data.open_netcdf_data(fpath)

    #timeInds = xrange(0,28,4)
    #timeInds = [0]
    nTimes = len(data.dimensions['Time'])
    #timeInds = xrange(0,nTimes,1)
    timeInds = xrange(26,nTimes,1)
    #timeInds = xrange(nTimes-1,-1,-1)
    
    #search for cyclone center
    cClose = conn.findOwner_horizNbrs_latLon(pt_ll, cClose, latCell, lonCell, nEdgesOnCell, cellsOnCell)

    cCloseList, pt_ll, minVals = driver_vortexLoc(data, timeInds, pt_ll, cClose, cellsOnCell, nEdgesOnCell, latCell, lonCell, c0IsMin)
    
    cClose = cCloseList[-1]
    cellList.extend(cCloseList)
    minValList.extend(minVals)
    
    data.close()
    
  print cellList
  print minValList
'''
def run_vortexLoc_withSeeds():
  #track vortex with center of search disk already specified. use with interact.py
  #to get guesses if end up with wrong object.
  seeds = [40419,26537,10552,1634,7338,20030,437,2896,11920,22053,19965,19961,28617]
  
  print "Uhoh, not coded yet\n"
  exit(0)
  
  #fnames = cfsrFiles2()
  #fnames = ['/arctic1/nick/cases/v1.0/x4/2week/tiedtke/x4.tiedtke.output.2006-07-24_12.00.00.nc',
  #          '/arctic1/nick/cases/v1.0/x4/2week/tiedtke/x4.tiedtke.output.2006-07-31_12.00.00.nc']
  #fnames = ['/arctic1/nick/cases/v1.0/x4/longer/x4.kf.output.2006-08-15_00.00.00.nc']
  fnames = sorted(glob.glob('/arctic1/nick/cases/cfsr/output/x4.cfsr.output.2006-07-*'))
  #fnames = ['/arctic1/nick/cases/cfsr/output/x4.cfsr.output.2006-08-15_00.00.00.nc']
  #fnames = sorted(glob.glob('/arctic1/nick/cases/vduda/x4.t.output.2006-08-1*'))
  #fnames = ['/arctic1/nick/cases/vduda/x7/x7.kf.output.2006-08-07_18.00.00.nc']
  
  pt_ll = np.empty(2,dtype='float')
  #pt_ll[0] = 90.*np.pi/180.; pt_ll[1] = 0.0
  #pt_ll[0] = 87.05*np.pi/180.; pt_ll[1] = 18.86*np.pi/180. #for 20060801 start time
  #pt_ll[0] = 85.38*np.pi/180.; pt_ll[1] = 352.97*np.pi/180. #for 2006080312 start time
  #pt_ll[0] = 86.1*np.pi/180.; pt_ll[1] = 137.*np.pi/180. #for 2006080718
  #pt_ll[0] = 82.5*np.pi/180.; pt_ll[1] = 41.95*np.pi/180. # for 2006081500
  pt_ll[0] = 71.55*np.pi/180.; pt_ll[1] = 153.1*np.pi/180. #for 2006072412

  #pt_ll[0] = 1.41; pt_ll[1] = 0.27
  #give initial guess of cClose. guess is just seed to search, so can choose anything.
  cClose = 71791 #131735 #x4 mesh for 20060801 start time
  #cClose = 102439 #x7 mesh for 20060801 start time
  c0IsMin = True;
  
  #load in mesh data
  fpath = fnames[0]
  data = output_data.open_netcdf_data(fpath)
  nEdgesOnCell = data.variables['nEdgesOnCell'][:];
  cellsOnCell = data.variables['cellsOnCell'][:]-1;
  latCell = data.variables['latCell'][:];
  lonCell = data.variables['lonCell'][:];
  data.close()
  
  cellList = []
  minValList = []
  for iFile, fpath in enumerate(fnames):
    data = output_data.open_netcdf_data(fpath)

    #timeInds = xrange(0,28,4)
    #timeInds = [0]
    nTimes = len(data.dimensions['Time'])
    timeInds = xrange(0,nTimes,1)
    #timeInds = xrange(2,28,4)
    
    #search for cyclone center
    cClose = conn.findOwner_horizNbrs_latLon(pt_ll, cClose, latCell, lonCell, nEdgesOnCell, cellsOnCell)

    cCloseList, pt_ll, minVals = driver_vortexLoc(data, timeInds, pt_ll, cClose, cellsOnCell, nEdgesOnCell, latCell, lonCell, c0IsMin)
    
    cClose = cCloseList[-1]
    cellList.extend(cCloseList)
    minValList.extend(minVals)
    
    data.close()
    
  print cellList
  print minValList
'''
def driver_vortexLoc(data, timeInds, pt_ll, cClose, cellsOnCell, nEdgesOnCell, latCell, lonCell, isMin):
  
  cCloseList = []
  valList = []
  for tInd in timeInds:
    #use theta on 2 pvu
    output_data.print_xtime(data, tInd)
    
    epv_ht, theta_trop = compare.calc_height_theta_2PVU(data, tInd)
    #theta_trop = data.variables['temperature_500hPa'][tInd,:]

    #radius = 1000.e3
    radius = 500.e3
    cClose = findLocalExtremum(pt_ll, cClose, radius, cellsOnCell, nEdgesOnCell, latCell, lonCell,
                                 theta_trop, isMin)
    pt_ll[0] = latCell[cClose]; pt_ll[1] = lonCell[cClose]
    print "{0} {1} {2}".format(pt_ll[0], pt_ll[1], theta_trop[cClose])
    cCloseList.append(cClose)
    valList.append(theta_trop[cClose])
    
    #estimate cyclone region properties
    #tInd=0
    #epv_ht, theta_trop = compare.calc_height_theta_2PVU(data, tInd)
    rShoot = calc_objectRadius_shoot(cClose, cellsOnCell, nEdgesOnCell, latCell, lonCell, theta_trop, isMin)
    cycloneCells = gather_regionGrow_nbrValues(cClose, cellsOnCell, nEdgesOnCell, theta_trop, isMin)
    print "Radius shoot and cells in cyclone region: ", rShoot, cycloneCells
    
  return (cCloseList, pt_ll, valList)

def extrap_objectNewLocation(latlon1, latlon2, dt12, dt2New):
  #given latitude and longitude in radians of 2 previous object locations as numpy arrays,
  #extrapolate to the next time using that "polar velocity".
  #Readjust location for motion around poles

  #For an extreme case against xyz velocity, consider what happens if the object
  #moves a half circle around a latitude

  #the following is dll_dt = (ll2-ll1)/dt12. rNew = rOld + u*dt with less loss of precision
  dll = (latlon2-latlon1)
  newLatLon = latlon2+dll*(dt2New/dt12)

  #correct for cases where object moves across poles just in case
  poleN = np.pi/2.0; poleS = -poleN
  if (newLatLon[0]>poleN):
    latDiff = newLatLon[0]-poleN; #if extrapolated to 93 degrees, should be 87
    newLatLon[0] = poleN-latDiff
  elif (newLatLon[0]<poleS):
    latDiff = poleS-newLatLon[0]; #if extrapolated to -93 degrees, should be -87
    newLatLon[0] = poleS+latDiff
    
  #could bound longitude in same fashion but it should wrap fine
  
  return newLatLon

def findLocalExtremum(pt_ll, cClose, radius, cellsOnCell, nEdgesOnCell, latCell, lonCell,
                     vals, isMin):
  #gather cells within a heuristic search area of point and cell.
  #The seed cell for searching the region needs to be within radius of the point.
  #return the cell with the most greatest or smallest value as specified.

  cellsRegion = conn.gatherCells_radius(pt_ll, radius, cClose, cellsOnCell, nEdgesOnCell, latCell, lonCell)
  cellInd = 0 #junk value not used
  if (isMin):
    cellInd = np.argmin(vals[cellsRegion])
  else:
    cellInd = np.argmax(vals[cellsRegion])
  return cellsRegion[cellInd]

def calc_objectRadius_shoot(c0, cellsOnCell, nEdgesOnCell, latCell, lonCell, vals, c0IsMin):
  #Given initial cell and a discrete surrounding field,
  #shoot lines out to try to characterize the local extrema.
  #Natural choice for search directions is each nbr, and we'll create a single
  #measure out of the various results.

  #cellsOnCell -= 1
  rEarth = 6371.e3
  nNbrs = nEdgesOnCell[c0]
  pt_ll = np.empty(2,dtype='float');
  pt_ll[0]=latCell[c0]; pt_ll[1]=lonCell[c0];
  radiusDirs = np.empty(nNbrs,dtype='float')
  nbrs0 = cellsOnCell[c0,0:nNbrs]
  for i in xrange(nNbrs):
    #walk away from this nbr until cross value. know center is extremum
    iNbr = nbrs0[i]; c0New = iNbr

    flag = True
    while (flag==True):
      nNbrsNbr = nEdgesOnCell[c0New]; nbrsNbr = cellsOnCell[c0New,0:nNbrsNbr]
      (iNbr, d) = conn.find_farthestCellToPoint(pt_ll, nbrsNbr, nNbrsNbr, latCell, lonCell)
      d = d*rEarth
      cFar = nbrsNbr[iNbr]
      val0 = vals[c0New]; valFar = vals[cFar]
      c0New = cFar
      if (c0IsMin): #quit if decrease outside
        flag = valFar<val0
      else: #quit if increase outside
        flag = valFar>val0
        
    radiusDirs[i] = d
  print "Shooting radius in each direction: ", radiusDirs
  rMedian = np.median(radiusDirs, overwrite_input=True) #if overwrite, input array gets partially sorted
  return rMedian

def gather_regionGrow_nbrValues(c0, cellsOnCell, nEdgesOnCell, vals, c0IsMin):
  #Given initial cell and a discrete surrounding field,
  #identify a region by a flood fill/region grow.
  #a cell is added to the flood fill based on a local condition:
  #-if it's value is above/below nbr avg as specified
  #-if more neighbors are above/below

  cellsRegion = [c0]
  oldLen=0; newLen=1;
  while(oldLen<newLen): #quit when stopped adding cells
    newInds = range(oldLen,newLen,1)
    oldLen = newLen

    for i in newInds: #add nbrs of just added
      c0 = cellsRegion[i]
      nNbrs = nEdgesOnCell[c0]
      nbrs = cellsOnCell[c0,0:nNbrs]
      for n in nbrs:
        nNbrNbrs = nEdgesOnCell[n]; nbrNbrs = cellsOnCell[n,0:nNbrNbrs]
        c0Val = vals[n]; valNbrs = np.mean(vals[nbrNbrs])
      if (c0IsMin): #add it if less than avg
        if (c0Val<valNbrs and n not in cellsRegion):
          cellsRegion.append(n)
      else: #add it if greater than avg
        if (c0Val>valNbrs and n not in cellsRegion):
          cellsRegion.append(n)

    newLen = len(cellsRegion)
  return cellsRegion
  
def calc_regionProperties(cells, areaCell, dvEdge, nEdgesOnCell, edgesOnCell, cellsOnEdge):
  #given a connected region of cells,
  
  #net area is sum of cell areas
  netArea = 0.0
  for iCell in cells:
    netArea = netArea + areaCell[iCell]
  
  #perimeter is defined by voronoi vertices shared by <=2 cells
  #if shared by 3 cells, have dual triangle on that corner.
  #have distances in dvEdge
  netPerim = 0.0
  for iCell in cells:
    ne = nEdgesOnCell[iCell]
    for e in xrange(ne):
      iEdge = edgesOnCell[e]
      if((cellsOnEdge[iEdge,0] not in cells) or (cellsOnEdge[iEdge,1] not in cells)):
        #it's a boundary
        netPerim = netPerim+dvEdge[iEdge]


  #hydraulic radius or half of hydraulic diameter depending on terminology
  #r=2*area/perim since 2(pi*r*r)/(2pi*r) = r
  #(a la Reynolds number for fluids through engineering pipes)
  r_h = 2.*netArea/netPerim

  #eccentricity could be longest edge wrt shortest edge
  
  return r_h

def driver_tracers():
  '''
  given .nc (MPAS NETCDF) file and initial location for tracer, we want to 
  advect that tracer with the wind velocity through all time steps.
  '''
  
  #file properties
  ncfname = '/home/nickszap/research/mpas/output.2010-10-23_00:00:00.nc' #input file
  data = output_data.open_netcdf_data(ncfname)
  
  tSteps = len(data.dimensions['Time'])
  vertLevels = len(data.dimensions['nVertLevels']);
  
  r = np.empty(3);
  r[0] = 50.; r[1]=50.; r[2] = 5000.;
  dt = 6.*60.*60.; #can parse xtime in file as actual timestep
  
  #loop through time steps
  vel = np.empty(3);
  for i in range(tSteps):
    (hCell, vCell) = findOwner_coord(r, data, vertLevels);
    vel[0]=data.variables['uReconstructX'][i,hCell,vCell]
    vel[1]=data.variables['uReconstructY'][i,hCell,vCell]
    vel[2]=data.variables['uReconstructZ'][i,hCell,vCell]
    r = integrate(r,u,dt);
    print r
    
  data.close()
  
def integrate(r, u, dt):
  '''
  given a position and velocity as numpy arrays, get new position at next time
  '''
  return(r+u*dt) #if these aren't numpy,netcdf4,..., + will append them!

def findSurfaceCell_bruteForce(r, data):
  '''
  Find the cell on a given spherical level that contains the x,y,z coordinate.
  '''
  #as a quick first try: 
  #1) we're going to assume the mesh is voronoi so the closest cell
  #  is also the cell that contains the point.
  #2) find closest by calculating 'distance' to all cell centers. a global argmin seems silly since we have a mesh!
  #   at the least, we could do the simple walk (especially since we're w/o boundaries)
  
  pass

def findOwner_coord(r, data, vLevels):
  #find the horizontal owner and then the owner in that column.
  #return (horizIndex, vertLevel)
  
  hCell = findOwner_horizNbrs(r, data, 0); #initial guess of cell 0 could be better?
  zInterface = data.variables['zgrid'][:]
  (l,dl) = output_data.calcIndexOfValue(r[2],zInterface[hCell,:], vLevels+1) #this might be buggy outside top and bottom boundaries
  return(hCell,l)

def findOwner_horizNbrs(r, data, cellId):
  '''
  find the cell that contains the xyz coordinate by walking in the direction
  of the closest neighbor (in spherical distance) to the initial guess. 
  Owner is the cell closer to point than any nbr (local min).
  '''
  
  #if we have the memory, we can avoid recomputing distances for nbrs of previous
  #nbrs by storing the whole mesh. say d=Inf until i calc it.
  
  nEdgesOnCell = data.variables['nEdgesOnCell'][:];
  cellsOnCell = data.variables['cellsOnCell'][:];
  latCell = data.variables['latCell'][:];
  lonCell = data.variables['lonCell'][:];
  
  #convert r to spherical coords
  rll = calc_xyz2Spherical(r);
  llNbr = np.empty(2);
  llNbr[0]=latCell[cellId]; llNbr[1]=lonCell[cellId];
  dCell = calc_distSphere(1., rll[1:], llNbr);
  
  flag = 1;
  while (flag==1):
    #keep going towards cell closer to point until no nbr is closer
    flag =0;
    nbrs = getNbrs_cell(cellId, nEdgesOnCell[cellId], cellsOnCell[cellId,:]);
    for i in range(nEdgesOnCell[cellId]):
      nbrId = nbrs[i];
      llNbr[0]=latCell[nbrId]; llNbr[1]=lonCell[nbrId];
      dNbr = calc_distSphere(1., rll[1:], llNbr);
      if (dNbr<dCell):
        dCell = dNbr;
        cellId = nbrId;
        flag=1;
      #running min
    #end of getting min of self and neighbors
  return(cellId);


def getNbrs_cell(cell, nEdgeOnCell, cellsOnCell):
  '''
  return the indices of all neighbors to a given cell
  '''
  #since know that have nEdgesOnCells neighbors since we don't have lateral boundaries,
  #can loop through cellsOnCell to get nbrs on level
  nNbrs = nEdgeOnCell[cell];
  #return (cellsOnCell[0:nNbrs-1]);
  return (cellsOnCell[0:nNbrs]-1); #0:6 give [0,1,...,5] I believe. We also change 1-index to 0
  pass

def calc_xyz2Spherical(r):
  '''
  return the (radius,lat,lon) of an xyz point as the mesher does
  '''
  s = np.empty(3)
  s[0] = np.sqrt(np.dot(r,r))
  s[1] = np.arcsin(r[2]/s[0]) #this is the mesher way
  #s[1] = np.arccos(r[2]/s[0]) #this is the wiki way
  
  #maybe beware divide by 0, but numpy handled some basic cases like atan2(-1,0), atan2(0,0) fine
  s[2] = np.arctan2(r[1],r[0])
  return(s)

def calc_distSphere(r, ll1, ll2):
  '''
  calculate the distance on the sphere between two latLon points on the same radius.
  could use haversine or vincenty formulas for better conditioning/more accuracy apparently.
  '''
  
  '''
  in some function of some part of mpas init...:
  arg1 = sqrt( sin(0.5*(lat2-lat1))**2 +  &
              cos(lat1)*cos(lat2)*sin(0.5*(lon2-lon1))**2 )
 sphere_distance = 2.*radius*asin(arg1)
  '''
  
  dll = ll2-ll1;
  latTerm = np.sin(.5*dll[0]); latTerm = latTerm*latTerm;
  lonTerm = np.sin(.5*dll[1]); lonTerm = lonTerm*lonTerm*cos(ll1[0])*cos(ll2[0]);
  dAngle = np.sqrt(latTerm+lonTerm)
  
  dist = 2.*r*np.arcsin(dAngle)
  
  return(dist)

def test_spherical(data):
  '''
  test whether our lat lon calcs match the lat lon of the mesh
  surprise, they don't!!! apparently lon is sometimes of by -2pi. correct that!!!!
  '''
  
  #cell locations on unit sphere
  xc = data.variables['xCell'][:]
  yc = data.variables['yCell'][:]
  zc = data.variables['zCell'][:]
  
  latc = data.variables['latCell'][:]
  lonc = data.variables['lonCell'][:]
  
  #calc rms difference for surface layer
  dims = xc.shape
  nCells = dims[0]
  rms = 0.
  xyz = np.empty(3)
  dlat = np.empty(nCells); dlon = np.empty(nCells)
  for i in range(nCells):
    xyz[0]=xc[i]; xyz[1]=yc[i]; xyz[2]=zc[i];
    rll = calc_xyz2Spherical(xyz);
    
    dlat[i] = rll[1]-latc[i]; dlon[i] = rll[2]-lonc[i];
    #it doesn't matter if we're off by multiples of 2pi
    eps = float(1e-6); #if we mod -.0001 by 2pi, it becomes 2pi-.0001. we want it to be 0
    dlat[i] = np.mod(dlat[i]+eps, 2.*np.pi); dlon[i]=np.mod(dlon[i]+eps,2.*np.pi);
    
    rms =rms + dlat[i]*dlat[i] + dlon[i]*dlon[i]
    
  rms = rms/nCells; rms = np.sqrt(rms)
  print "RMS of latLon between mesh and calc is: ",rms

  import matplotlib.pyplot as plt
  plt.plot(dlat,'b', dlon, 'g');
  plt.show()
  return (dlat,dlon)

def cfsrFiles():
  s = '''/arctic1/nick/cases/cfsr/output/x4.cfsr.output.2006-08-01_00.00.00.nc
/arctic1/nick/cases/cfsr/output/x4.cfsr.output.2006-08-02_00.00.00.nc
/arctic1/nick/cases/cfsr/output/x4.cfsr.output.2006-08-03_00.00.00.nc
/arctic1/nick/cases/cfsr/output/x4.cfsr.output.2006-08-04_00.00.00.nc
/arctic1/nick/cases/cfsr/output/x4.cfsr.output.2006-08-05_00.00.00.nc
/arctic1/nick/cases/cfsr/output/x4.cfsr.output.2006-08-06_00.00.00.nc
/arctic1/nick/cases/cfsr/output/x4.cfsr.output.2006-08-07_00.00.00.nc
/arctic1/nick/cases/cfsr/output/x4.cfsr.output.2006-08-08_00.00.00.nc
/arctic1/nick/cases/cfsr/output/x4.cfsr.output.2006-08-09_00.00.00.nc
/arctic1/nick/cases/cfsr/output/x4.cfsr.output.2006-08-10_00.00.00.nc
/arctic1/nick/cases/cfsr/output/x4.cfsr.output.2006-08-11_00.00.00.nc
/arctic1/nick/cases/cfsr/output/x4.cfsr.output.2006-08-12_00.00.00.nc
/arctic1/nick/cases/cfsr/output/x4.cfsr.output.2006-08-13_00.00.00.nc
/arctic1/nick/cases/cfsr/output/x4.cfsr.output.2006-08-14_00.00.00.nc
/arctic1/nick/cases/cfsr/output/x4.cfsr.output.2006-08-15_00.00.00.nc
'''.split()
  return s

def cfsrFiles2():
  s = '''/arctic1/nick/cases/cfsr/output/x4.cfsr.output.2006-07-25_00.00.00.nc
/arctic1/nick/cases/cfsr/output/x4.cfsr.output.2006-07-26_00.00.00.nc
/arctic1/nick/cases/cfsr/output/x4.cfsr.output.2006-07-27_00.00.00.nc
/arctic1/nick/cases/cfsr/output/x4.cfsr.output.2006-07-28_00.00.00.nc
/arctic1/nick/cases/cfsr/output/x4.cfsr.output.2006-07-29_00.00.00.nc
/arctic1/nick/cases/cfsr/output/x4.cfsr.output.2006-07-30_00.00.00.nc
/arctic1/nick/cases/cfsr/output/x4.cfsr.output.2006-07-31_00.00.00.nc
/arctic1/nick/cases/cfsr/output/x4.cfsr.output.2006-08-01_00.00.00.nc
/arctic1/nick/cases/cfsr/output/x4.cfsr.output.2006-08-02_00.00.00.nc
/arctic1/nick/cases/cfsr/output/x4.cfsr.output.2006-08-03_00.00.00.nc
/arctic1/nick/cases/cfsr/output/x4.cfsr.output.2006-08-04_00.00.00.nc
/arctic1/nick/cases/cfsr/output/x4.cfsr.output.2006-08-05_00.00.00.nc
/arctic1/nick/cases/cfsr/output/x4.cfsr.output.2006-08-06_00.00.00.nc
/arctic1/nick/cases/cfsr/output/x4.cfsr.output.2006-08-07_00.00.00.nc
/arctic1/nick/cases/cfsr/output/x4.cfsr.output.2006-08-08_00.00.00.nc
'''.split()
  return s

'''
Routine from the module_sphere_utilities.F in mesher to check whether we're working with the same lat/lon:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SUBROUTINE CONVERT_XL
!
! Convert (x, y, z) to a (lat, lon) location on a sphere with 
!    radius sqrt(x^2 + y^2 + z^2).
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine convert_xl(x, y, z, latlon)

   use data_types

   implicit none

   real, intent(in) :: x, y, z
   type (geo_point), intent(out) :: latlon

   real :: dl, clat, pii, rtod
   real :: eps
   parameter (eps=1.e-10)

   pii = 2.*asin(1.0)
   rtod=180./pii
   dl = sqrt(x*x + y*y + z*z)

   latlon%lat = asin(z/dl)

!  check for being close to either pole

   if (abs(x) > eps) then

      if (abs(y) > eps) then

         latlon%lon = atan(abs(y/x))

         if ((x <= 0.) .and. (y >= 0.)) then
            latlon%lon = pii-latlon%lon
         else if ((x <= 0.) .and. (y < 0.)) then
            latlon%lon = latlon%lon+pii
         else if ((x >= 0.) .and. (y <= 0.)) then
            latlon%lon = 2*pii-latlon%lon
         end if

      else ! we're either on longitude 0 or 180

         if (x > 0) then
            latlon%lon = 0.
         else
            latlon%lon = pii
         end if

      end if

   else if (abs(y) > eps) then  

      if (y > 0) then 
         latlon%lon = pii/2.
      else
         latlon%lon = 3.*pii/2.
      end if

   else  ! we are at a pole

      latlon%lon = 0.

   end if

end subroutine convert_xl


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FUNCTION SPHERE_DISTANCE
!
! Given two (latitude, longitude) coordinates on the surface of a sphere,
!    plus the radius of the sphere, compute the great circle distance between
!    the points.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real function sphere_distance(p1, p2, radius)

   use data_types

   implicit none

   type (geo_point), intent(in) :: p1, p2
   real, intent(in) :: radius

   real :: arg1

   arg1 = sqrt( sin(0.5*(p2%lat-p1%lat))**2 +  &
                cos(p1%lat)*cos(p2%lat)*sin(0.5*(p2%lon-p1%lon))**2 )
   sphere_distance = 2.*radius*asin(arg1)

end function sphere_distance
'''

if __name__=='__main__':
  #driver_tracers()
  run_vortexLoc()
  