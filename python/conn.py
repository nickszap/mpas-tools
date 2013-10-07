# -*- coding: utf-8 -*-

# imports
import numpy as np
import netCDF4

import output_data

def findOwner_coord(r, data, vLevels):
  #find the horizontal owner and then the owner in that column.
  #return (horizIndex, vertLevel)
  
  hCell = findOwner_horizNbrs(r, data, 0); #initial guess of cell 0 could be better?
  zInterface = data.variables['zgrid'][:]
  (l,dl) = output_data.calcIndexOfValue(r[2],zInterface[hCell,:], vLevels+1) #this might be buggy outside top and bottom boundaries
  return(hCell,l)

def find_closestCellToPoint(pt_ll, cells, nCells, latCell, lonCell):
  #find the cell closest to the point defined by lat-lon.
  #return index of cell in cells and distance
  minId = 0; dMin = 1.e10 #init as big number
  ll_cell = np.empty(2);
  for i in xrange(nCells):
    cellId = cells[i]
    ll_cell[0]=latCell[cellId]; ll_cell[1]=lonCell[cellId];
    dCell = calc_distSphere(1., pt_ll, ll_cell);
    if (dCell<dMin): #running min
      minId = i; dMin = dCell;
  return (minId, dMin)

def find_farthestCellToPoint(pt_ll, cells, nCells, latCell, lonCell):
  #find the cell closest to the point defined by lat-lon.
  #return index of cell in cells and distance
  minId = 0; dMax = -1.e10 #init as big number
  ll_cell = np.empty(2);
  for i in xrange(nCells):
    cellId = cells[i]
    ll_cell[0]=latCell[cellId]; ll_cell[1]=lonCell[cellId];
    dCell = calc_distSphere(1.0, pt_ll, ll_cell);
    if (dCell>dMax): #running min
      minId = i; dMax = dCell;
  return (minId,dMax)

def find_nbrInDirection(vec, xyz0, nbrs, xCell, yCell, zCell):
  #find the nbr most aligned with the specified vector wrt specified coordinate.
  #Use xyz space under assumption of small curvature since polar is a headache
  #around the pole.
  #Return the nbr

  nNbrs = len(nbrs)
  #vectors from point to nbrs
  pt2Nbr = np.empty((nNbrs,3))
  for i,nbr in enumerate(nbrs):
    pt2Nbr[i,0] = xCell[nbr]-xyz0[0]
    pt2Nbr[i,1] = yCell[nbr]-xyz0[1]
    pt2Nbr[i,2] = zCell[nbr]-xyz0[2]
    #normalize so actual distance from point is not a factor
    mag = np.linalg.norm(pt2Nbr[i,:])
    pt2Nbr[i,:] /= mag
  #

  mag = np.empty(nNbrs)
  for i in xrange(nNbrs):
    mag[i] = np.dot(vec,pt2Nbr[i,:])
  #
  ind = np.argmax(mag)
  return(nbrs[ind])

def gatherCellsOnLine(c0, ll, cellsOnCell, nEdgesOnCell, latCell, lonCell):
  #the "line" is defined by cell0 and latLon of point outside cell.
  #return a list of cells crossed to point.
  
  cOnLine = [c0]; dMin=1.e10
  while (True):
    nNbrs = nEdgesOnCell[c0]
    nbrs = cellsOnCell[c0,0:nEdgesOnCell[c0]]
    (ind,dist) = find_closestCellToPoint(ll, nbrs, nNbrs, latCell, lonCell)
    c0 = nbrs[ind]
    if (dist<dMin):
      cOnLine.append(c0); dMin = dist;
    else:
      break
  return cOnLine

def gatherCells_radius(pt_ll, radius, c0, cellsOnCell, nEdgesOnCell, latCell, lonCell):
  #using seed cell, gather all points with cell centers within radius of the specified point.
  #return a list of those cells

  candidates = [c0]
  closeCells = []; farCells=[]
  ll_cell = np.empty(2);
  while (len(candidates)>0):
    c0 = candidates.pop(0)
    ll_cell[0]=latCell[c0]; ll_cell[1]=lonCell[c0];
    #dCell = calc_distSphere(1., pt_ll, ll_cell);
    dCell = calc_distSphere(6371000., pt_ll, ll_cell);

    if (dCell<radius):
      closeCells.append(c0)
      nbrs = cellsOnCell[c0,0:nEdgesOnCell[c0]]
      for n in nbrs:
        if ((n not in closeCells) and (n not in farCells) and (n not in candidates)):
          candidates.append(n)
    else:
      farCells.append(c0)

  return closeCells

def floodFillRegion(flags, bFlag, c0, cellsOnCell,nEdgesOnCell):
  #given a connected boundary (flags[i]=bFlag) and a cell defining inside the boundary,
  #and flags[i] set not to flags[c0] or bFlag,
  #fill flags[i] for all i connected to inside cell with flags[c0]
  val0 = flags[c0]
  nbrs = cellsOnCell[c0,0:nEdgesOnCell[c0]]
  candidates = nbrs.tolist()
  while(len(candidates)>0):
    c1 = candidates.pop(0)
    if (flags[c1]!=bFlag and flags[c1]!=val0):#nbr of c0, not boundary, not visited
      flags[c1] = val0;
      nbrs = cellsOnCell[c1,0:nEdgesOnCell[c1]]
      for i in nbrs:
        if (flags[i]!=bFlag and flags[i]!=val0):
          candidates.append(i)
  #
  return flags

def gatherVerticesInRegion(nCells, verticesOnCell, nEdgesOnCell, nVerticesTotal, countThreshold):
  #given connected cells that define a region, have local verticesOnCell and nEdgesOnCell.
  #a vertex is in a connected region if it's shared by countThreshold cells.
  #the return is used as an array index for the proper vertices

  vCount = np.zeros(nVerticesTotal,dtype=int)
  for cell in xrange(nCells):
    nVerts = nEdgesOnCell[cell]
    verts = verticesOnCell[cell,0:nVerts]#-1
    for v in verts:
      vCount[v]+=1
  vIn = np.array(xrange(nVerticesTotal))[vCount>(countThreshold-1)] #vcount >= threshold
  return vIn

def testFloodFill():
  #flood fill a region of map.  
  #we'll viz by overwriting a section in the file
  ncfname = '/home/nickszap/research/mpas/output.2010-10-23_00:00:00.nc'
  data = netCDF4.Dataset(ncfname,'a')
  nEdgesOnCell = data.variables['nEdgesOnCell'][:];
  cellsOnCell = data.variables['cellsOnCell'][:]-1;
  latCell = data.variables['latCell'][:];
  lonCell = data.variables['lonCell'][:];
  
  #we'll overwrite ke at level 33
  t=0; l=33;
  ke = data.variables['ke'][t,:,l] #can access like data.variables['ke'][0][2][[4,5]] too
  
  #make a bounding box ----------------------
  ll1 = np.array([49.2505,-123.1119])*np.pi/180. #vancouver
  ll2 = np.array([32.7153,-117.1564])*np.pi/180. #san diego
  ll3 = np.array([33.7489,-84.3881])*np.pi/180. #atlanta
  ll4 = np.array([45.5081,-73.5550])*np.pi/180. #montreal
  llc = np.array([43.6136,-116.2025])*np.pi/180. #boise
  
  #gather cells on that box -------
  cellsOnBox = []
  #find owner of corner then 1->2->3...
  cOnLine = gatherCellsOnLine(0, ll1, cellsOnCell, nEdgesOnCell, latCell, lonCell) #initial guess of cell 0
  own1 = cOnLine[-1] #closest cell to ll1
  cellsOnBox.append(own1)
  
  cOnLine = gatherCellsOnLine(own1, ll2, cellsOnCell, nEdgesOnCell, latCell, lonCell)
  own2 = cOnLine[-1]
  cellsOnBox+= cOnLine
  cOnLine = gatherCellsOnLine(own2, ll3, cellsOnCell, nEdgesOnCell, latCell, lonCell)
  own3 = cOnLine[-1]
  cellsOnBox+= cOnLine
  cOnLine = gatherCellsOnLine(own3, ll4, cellsOnCell, nEdgesOnCell, latCell, lonCell)
  own4 = cOnLine[-1]
  cellsOnBox+= cOnLine
  cOnLine = gatherCellsOnLine(own4, ll1, cellsOnCell, nEdgesOnCell, latCell, lonCell)
  cellsOnBox+= cOnLine[0:-1] #already included start. could just list.pop() at end of loop
  endCell = cOnLine[-1]
  if (own1 != endCell):
    print "Uhoh. Walked in loop and start={0} is not end={1}\n".format(own1,endCell)
  
  #for say acute angles, we can cut back through a cell.
  #we want unique cells in list and don't need to preserve order
  bCells = list(set(cellsOnBox))
  
  #flood fill interior------------
  valInit=-1; val0=42; valBoundary=2; 
  ke[:] = valInit;  
  #find owner of center
  cOnLine = gatherCellsOnLine(own1, llc, cellsOnCell, nEdgesOnCell, latCell, lonCell)
  ownc = cOnLine[-1]
  #set flags for center and boundary
  ke[bCells]=valBoundary
  ke[ownc]=val0
  
  ke = floodFillRegion(ke, valBoundary, ownc, cellsOnCell,nEdgesOnCell)
  
  #write out result
  data.variables['ke'][t,:,l] = ke[:]
  data.close()
  
def findOwner_horizNbrs_cartesian(r, data, cellId):
  '''
  find the cell that contains the xyz coordinate by walking in the direction
  of the closest neighbor (in spherical distance) to the initial guess. 
  Owner is the cell closer to point than any nbr (local min).
  '''
  
  #if we have the memory, we can avoid recomputing distances for nbrs of previous
  #nbrs by storing the whole mesh. say d=Inf until i calc it.
  
  nEdgesOnCell = data.variables['nEdgesOnCell'][:];
  cellsOnCell = data.variables['cellsOnCell'][:]-1;
  latCell = data.variables['latCell'][:];
  lonCell = data.variables['lonCell'][:];
  
  #convert r to spherical coords
  rll = calc_xyz2Spherical(r);
  
  cellId = findOwner_horizNbrs_latLon(rll[1:], cellId, latCell, lonCell, nEdgesOnCell, cellsOnCell)
  return cellId

def findOwner_horizNbrs_latLon(pt_ll, cellId, latCell, lonCell, nEdgesOnCell, cellsOnCell):

  llNbr = np.empty(2);
  llNbr[0]=latCell[cellId]; llNbr[1]=lonCell[cellId];
  dCell = calc_distSphere(1., pt_ll, llNbr);

  flag = 1;
  while (flag==1):
    #keep going towards cell closer to point until no nbr is closer
    flag =0;
    nNbrs = nEdgesOnCell[cellId]
    nbrs = cellsOnCell[cellId,0:nNbrs]
    #nbrs = getNbrs_cell(cellId, nEdgesOnCell[cellId], cellsOnCell[cellId,:]);
    for i in range(nNbrs):
      nbrId = nbrs[i];
      llNbr[0]=latCell[nbrId]; llNbr[1]=lonCell[nbrId];
      dNbr = calc_distSphere(1., pt_ll, llNbr);
      if (dNbr<dCell):
        dCell = dNbr;
        cellId = nbrId;
        flag=1;
      #running min
    #end of getting min of self and neighbors
  return(cellId);

def getNbrs_cell(cell, nEdgeOnCell, cellsOnCell):
  '''
  return the indices of horizontal neighbors to a given cell
  '''
  #since know that have nEdgesOnCells neighbors since we don't have lateral boundaries,
  #can loop through cellsOnCell to get nbrs on level
  nNbrs = nEdgeOnCell[cell];
  #return (cellsOnCell[0:nNbrs-1]);
  return cellsOnCell[cell,0:nNbrs]; #0:6 give [0,1,...,5] I believe. We also change 1-index to 0-based

def getHorizNbrs(data, hCell):
  cellsOnCell = data.variables['cellsOnCell'][hCell] #this has "junk" to maxEdges dimension
  nNbrs = data.variables['nEdgesOnCell'][hCell]
  return (cellsOnCell[0:nNbrs]-1) #0-indexing

def get_allHorizNbrs(data, nCells):
  cellsOnCell = [data.variables['cellsOnCell'][i][0:data.variables['nEdgesOnCell'][i]]-1 for i in range(nCells)]
  return cellsOnCell

def get_domainHorizNbrs(data, hCells, g2lMap):
  #nbr connectivity within the local domain (w/ halos?)
  nEdgesOnCell = data.variables['nEdgesOnCell'][hCells]
  cellsOnCell = data.variables['cellsOnCell'][hCells,:]-1
  
  #convert to localInd
  cellsOnCell = make_localDomainNbrs(len(hCells), cellsOnCell, nEdgesOnCell, g2lMap)
  
  return cellsOnCell

#meshes can be too large for RAM. we'll use domain decomposition and cycle over domains (serially) to reduce cost  
def partition_max(seed0, cellsOnCell, nEdgesOnCell, nCells, nSeeds):
  #perfect load balancing means each site gets same number of cells.
  #we won't get that, but if use number of cells to site as distance metric,
  #get natural region growing procedure.
  #return (map cell->site, cell sites). site should also be map[c]=c
  
  d2Site = np.ones(nCells,dtype=int)*int(1e10) #big so gets changed on first seed
  closestSite = np.ones(nCells,dtype=int)*-1 #should all be >=0 in returned
  
  seeds = [seed0]; #nSeeds=100;
  newSeed=seed0; #maybe max density cell is better choice?
  while (len(seeds)<=nSeeds): #use <= so last seed fills its nbrs
    #flood fill out with this new seed to tag closest
    d2Site[newSeed]=0; closestSite[newSeed]=newSeed;
    added = [newSeed]
    while(len(added)>0):
      c = added.pop(0); dist = d2Site[c]+1 #go through c to get to nbr cell
      nNbrs = nEdgesOnCell[c]
      nbrs = cellsOnCell[c,0:nNbrs]
      for n in nbrs:
        if (dist<d2Site[n]):
          #print "Cell {0} updating from seed {1} to {2} from dist {3} to dist {4}\n".format(n, closestSite[n], newSeed,d2Site[n],dist)
          closestSite[n]=newSeed; d2Site[n]=dist
          added.append(n)
    
    #make a new seed by taking the farthest distance cell.
    #seems like a logical, global heuristic but is there any guarantee?!!!?
    newSeed = np.argmax(d2Site)
    seeds.append(newSeed); #print "partition seed: "+str(newSeed)+'\n'
  
  return (closestSite,seeds[0:-1]) #since appended extra 1 for last seed to fill nbrs

def test_partition():
  ncfname = '/home/nickszap/research/mpas/output.2010-10-23_00:00:00.nc'
  data = netCDF4.Dataset(ncfname,'a')
  nCells = len(data.dimensions['nCells'])
  nEdgesOnCell = data.variables['nEdgesOnCell'][:];
  cellsOnCell = data.variables['cellsOnCell'][:]-1;
  
  nSeeds = 10; seed0 = 0 #maybe max density cell is better choice?
  closestSite,seeds = partition_max(seed0, cellsOnCell, nEdgesOnCell,nCells, nSeeds)
  print seeds
  
  #we'll overwrite ke at level 33
  t=0; l=33;
  closestSite[seeds]=-10 #make seeds stand out
  data.variables['ke'][t,:,l] = closestSite[:]
  data.close()

def getHalo(seed, closestSite, cellsOnCell, nEdgesOnCell, nCells):
  #we horizontally partition into connected domains. if we want to calculate based on all nbrs,
  #we'll also need the neighbors that are in other domains.
  #return the cellInds in the halo
  
  #loop through cells in domain. if nbr not in partition, it's in the halo.
  halo = []
  for c in xrange(nCells):
    if (closestSite[c]==seed): #cell in domain
      nNbrs = nEdgesOnCell[c]
      nbrs = cellsOnCell[c,0:nNbrs]
      for n in nbrs:
        if (closestSite[n]!=seed):
          halo.append(n)
  return list(set(halo))

def make_global2localMap(cellsDomain, cellsHalo, nCells):
  #for domain, want to have cellId be index into array. this really affects the nbr connectivity
  #return a global to local index map that's only valid for domain and halo cells
  
  ''' #any worthwhile way to save memory? 
  #only need entries up to largest id for array access as g2lmap[globalID]
  nCellsD = np.max(cellsDomain); nCellsH = np.max(cellsHalo)
  nCells = np.max((nCellsD,nCellsH))
  '''
  g2lMap = np.zeros(nCells,dtype=int)
  lInd=0 #local index
  for c in cellsDomain:
    g2lMap[c]=lInd; lInd=lInd+1
  for c in cellsHalo:
    g2lMap[c]=lInd; lInd=lInd+1
  
  return g2lMap
  
def make_localDomainNbrs(nCells, cellsOnCell_local, nEdgesOnCell_local, g2lMap):
  #given global nbr indices and a mapping to local cell inds for this domain,
  #reassign the neighbor map sliced for the proper cells with local indices.
  
  for c0Local in xrange(nCells):
    #c0Local = g2lMap[c0]
    nNbrs = nEdgesOnCell_local[c0Local]
    nbrs = cellsOnCell_local[c0Local,0:nNbrs]
    #for i in xrange(nNbrs):
      #cellsOnCell_local[c0Local,i] = g2lMap[nbrs[i]]
    cellsOnCell_local[c0Local,0:nNbrs] = g2lMap[nbrs][:]
  return cellsOnCell_local

def gatherArcticCells(latCell, nCells, latThresh):
  #return the cells above a given latitude threshold (in radians)
  cells = np.array(xrange(nCells))[latCell>latThresh]
  return cells

def get_arcticHalo(cells, latCell, latThresh, cellsOnCell, nEdgesOnCell):
  #a nbr of cell above latitude that is not above lat
  haloCells = []
  for c in cells:
    nNbrs = nEdgesOnCell[c]
    nbrs = cellsOnCell[c,0:nNbrs]
    for n in nbrs:
      if (latCell[n]<=latThresh):
        haloCells.append(n)
  return list(set(haloCells))

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
  input is lat-lon
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
  lonTerm = np.sin(.5*dll[1]); lonTerm = lonTerm*lonTerm*np.cos(ll1[0])*np.cos(ll2[0]);
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
  print "Nope. Call as helper fcts\n"

