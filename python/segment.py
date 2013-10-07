import numpy as np
import netCDF4
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

import compare
import output_data

#watershed has a few options for implementation:
#-for every cell, walk down steepest gradient to the basin
#-every local min is a voronoi map with distance as the sum of ?
#along path to min/voronoi site. benefit of this is ability to remove voronoi sites for iterative merging.
#could remove based on min size of watershed,...

#active contours to identify regions:
#-minimize: internal energy + external energy + curveCost(eg smoothness)
#can create logically square say nearest nbr resample of area N of 75N and ship over to image processing+

def find_minCells(cellsOnCell, nEdgesOnCell, vals, nCells):
  #return array[nCells] with 1 if cell is min

  isMin = np.zeros(nCells,dtype=int)
  for iCell in xrange(nCells):
    nNbrs = nEdgesOnCell[iCell]
    nbrs = cellsOnCell[iCell,0:nNbrs]
    valNbrs = vals[nbrs]
    val0 = vals[iCell]
    if (False not in (val0<=valNbrs)): #is site if can't descend from it
      isMin[iCell] = 1

  nSites = np.sum(isMin)
  print "Number of basin sites: ", nSites
  return isMin

def watershed(vals, cellIsMin, cellsOnCell, nEdgesOnCell, edgesOnCell, dcEdge, nCells):
  #to make adding/deleting basins simple, follow gradient until reach a site.
  #map every cell to follow local steepest gradient. basins go to self.
  #return map of cell to basin

  cell2Site = np.zeros(nCells,dtype=int)
  for iCell in xrange(nCells):
    if (cellIsMin[iCell]>0):
      cell2Site[iCell]= iCell

  #get local steepest path
  for iCell in xrange(nCells):
    if (cellIsMin[iCell]>0):
      continue
    nNbrs = nEdgesOnCell[iCell]
    nbrs = cellsOnCell[iCell,0:nNbrs]
    valNbrs = vals[nbrs]
    val0 = vals[iCell]

    #correspondence is towards minimum gradient.
    for iNbr in xrange(nNbrs):
      iEdge = edgesOnCell[iCell,iNbr]
      dx = dcEdge[iEdge]
      valNbrs[iNbr] = (valNbrs[iNbr]-val0)/dx

    iNbr = np.argmin(valNbrs)
    cell2Site[iCell] = nbrs[iNbr]

  #follow local steepest path to site
  for iCell in xrange(nCells):
    nextCell = cell2Site[iCell]
    while (not cellIsMin[nextCell]>0):
      nextCell = cell2Site[nextCell]
      #print "Cell {0} going to {1}".format(iCell,nextCell)

    cell2Site[iCell] = nextCell

  return cell2Site

def plotSegment(lat, lon, var, isMarker):  
  #map = Basemap(projection='ortho',lon_0=100,lat_0=60, resolution='l')
  map = Basemap(projection='ortho',lon_0=0,lat_0=90, resolution='l')
  x,y = map(lon, lat)
  
  fig1 = plt.figure(1)
  map.drawcoastlines()
  map.drawmapboundary()
  map.pcolor(x,y,var,tri=True,shading='flat',edgecolors='none',cmap=plt.cm.jet) #cmap=plt.cm.hot_r) #,vmin=100,vmax=1000)
  
  xMarker = x[isMarker>0]
  yMarker = y[isMarker>0]
  map.scatter(xMarker,yMarker,marker="o")
  
  #plt.colorbar()
  plt.show()

#filter field
def smooth_mean(cellsOnCell, nEdgesOnCell, valsIn, nCells):
  #val[c0] is average of neighbor values
  
  vals = np.empty_like (valsIn)
  #np.copyto(vals, valsIn)
  vals[:] = valsIn[:]
  
  for iCell in xrange(nCells):
    nNbrs = nEdgesOnCell[iCell]
    nbrs = cellsOnCell[iCell,0:nNbrs]
    valNbrs = valsIn[nbrs]
    val0 = valsIn[iCell]
    
    vals[iCell] = (nNbrs*np.mean(valNbrs)+val0)/(nNbrs+1)
    
  return vals

#reconstruction

def example_segment():
  #we want to associate features with cyclones and anticyclones.
  #we'll create basins for local minima and maxima separately so each cell
  #can be part of 2 basins. To avoid over-segmentation, smooth field first.
  #Then, we map the cell to the proper type of basin based on additional information (eg vorticity)

  ncfname = '/arctic1/nick/cases/v1.0/x4/longer/x4.kf.output.2006-08-15_00.00.00.nc'
  data = netCDF4.Dataset(ncfname,'r')

  nCells = len(data.dimensions['nCells'])
  nEdgesOnCell = data.variables['nEdgesOnCell'][:];
  cellsOnCell = data.variables['cellsOnCell'][:]-1;
  edgesOnCell = data.variables['edgesOnCell'][:]-1;
  dcEdge = data.variables['dcEdge'][:]
  
  #field
  #tInd = 0
  tInd = 12
  epv_ht, theta_trop = compare.calc_height_theta_2PVU(data, tInd)
  #vort_trop, theta_trop = calc_vorticity_theta_2PVU(data, t0)
  #segment based on min or max
  theta_trop = -theta_trop
  
  
  #smooth
  nSmooth = 10
  for iSmooth in xrange(nSmooth):
    theta_trop = smooth_mean(cellsOnCell, nEdgesOnCell, theta_trop, nCells)
  
  #segment
  cellIsMin = find_minCells(cellsOnCell, nEdgesOnCell, theta_trop, nCells)
  cell2Site = watershed(theta_trop, cellIsMin, cellsOnCell, nEdgesOnCell, edgesOnCell, dcEdge, nCells)

  latCell = data.variables['latCell'][:];
  lonCell = data.variables['lonCell'][:];

  latCell *= 180./np.pi
  lonCell *= 180./np.pi

  plotSegment(latCell, lonCell, cell2Site, cellIsMin)
  
  data.close()

def segment_high_low_watershed(data, tInd):
  #get high and low basin seeds, associate cells to both high and low basins if not extrema.
  #to decide whether "really" part of high or low basin, we have options:
  #-(anti-)cyclonic for (high) low...is local vorticity noisy?
  #-closer theta value to maxima a la color scale grouping...huge min or max value now matters
  #-whether steeper gradient is to high or low
  #-physical distance

  nCells = len(data.dimensions['nCells'])
  nEdgesOnCell = data.variables['nEdgesOnCell'][:];
  cellsOnCell = data.variables['cellsOnCell'][:]-1;
  edgesOnCell = data.variables['edgesOnCell'][:]-1;
  dcEdge = data.variables['dcEdge'][:]

  #field
  #tInd = 0
  #tInd = 12
  vort_trop, theta_trop = calc_vorticity_theta_2PVU(data, tInd)

  #smooth
  nSmooth = 10
  for iSmooth in xrange(nSmooth):
    theta_trop = smooth_mean(cellsOnCell, nEdgesOnCell, theta_trop, nCells)
  '''  
  printCells = [0,79440,114005]
  print theta_trop[printCells]
  '''
  #segment
  #mins
  cellIsMin = find_minCells(cellsOnCell, nEdgesOnCell, theta_trop, nCells)
  cell2SiteMin = watershed(theta_trop, cellIsMin, cellsOnCell, nEdgesOnCell, edgesOnCell, dcEdge, nCells)

  #maxs
  theta_trop = -theta_trop
  cellIsMax = find_minCells(cellsOnCell, nEdgesOnCell, theta_trop, nCells)
  cell2SiteMax = watershed(theta_trop, cellIsMax, cellsOnCell, nEdgesOnCell, edgesOnCell, dcEdge, nCells)
  
  #"voting" procedure for criteria
  cell2Site = np.empty(nCells, dtype=int)
  for iCell in xrange(nCells):
    if (cellIsMin[iCell]>0 or cellIsMax[iCell]>0):
      cell2Site[iCell] = iCell
    else:
      if (vort_trop[iCell]<0): #anticyclonic
        cell2Site[iCell] = cell2SiteMax[iCell]
      else: #0 or cyclonic
        cell2Site[iCell] = cell2SiteMin[iCell]

  return cell2Site

def calc_vorticity_theta_2PVU(data, t0):
  #vertical height of sign(lat)*2PVU surface.
  #return height of cell center in column using interpolation from top

  pvuTrop = 2.0

  nCells = len(data.dimensions['nCells'])
  nLevels = len(data.dimensions['nVertLevels'])
  epv = data.variables['ertel_pv'][t0,:,:]
  latCell = data.variables['latCell'][:]

  theta = data.variables['theta'][t0,:,:]
  vortVertex = data.variables['vorticity'][t0,:,:]
  nVertexOnCell = data.variables['nEdgesOnCell'][:]
  vertexOnCell = data.variables['verticesOnCell'][:,:]-1

  vort_trop = np.empty(nCells, dtype='float')
  theta_trop = np.empty(nCells,dtype='float')
  for iCell in xrange(nCells):
    pvuVal = pvuTrop
    if (latCell[iCell]<0):
      pvuVal = -pvuTrop

    #don't really trust top and bottom levels so don't count those for interpolation
    lev0 = 1
    interpLevs = range(lev0,nLevels-1); nInterpLevs = len(interpLevs)
    (l,dl) = output_data.calcIndexOfValue(pvuVal,epv[iCell,interpLevs], nInterpLevs)
    #print "Cell {0} has l,dl = {1},{2}".format(i, l, dl)
    theta_trop[iCell] = output_data.calcValueOfIndex(l,dl,theta[iCell,interpLevs])
    
    if (l<0): #missing val
      vort_trop[iCell] = output_data.missingVal
    else:
      #average vorticity to cell centers and interpolate to trop
      vortl = form_vertVorticityCell(iCell, lev0+l, vortVertex, vertexOnCell, nVertexOnCell)
      vorth = form_vertVorticityCell(iCell, lev0+l+1, vortVertex, vertexOnCell, nVertexOnCell)

      dvort_dl = vorth-vortl;
      vort_trop[iCell] = vortl+dvort_dl*dl
    
  return (vort_trop, theta_trop)


def form_vertVorticityCell(iCell, iLevel, vortVertex, vertexOnCell, nVertexOnCell):
  #approx way of getting vorticity at cell center from voronoi cell vertices

  nVerts = nVertexOnCell[iCell]
  verts = vertexOnCell[iCell,0:nVerts]
  vortVerts = vortVertex[verts,iLevel]
  vortMean = np.mean(vortVerts)
  return vortMean

'''
real(kind=RKIND) function calc_verticalVorticity_cell(c0, level, nVerticesOnCell, verticesOnCell, cellsOnVertex, &
                                                         kiteAreasOnVertex, areaCell, vVortVertex)
      !area weighted average of vorticity at vertices (really midpts of faces) to cell center for the specified cell
      !
      implicit none

      real(kind=RKIND), intent(in) :: areaCell
      integer, intent(in) :: c0, level, nVerticesOnCell
      integer, dimension(:,:), intent(in) :: verticesOnCell, cellsOnVertex
      real(kind=RKIND), dimension(:,:), intent(in) :: kiteAreasOnVertex, vVortVertex

      real(kind=RKIND) :: vVortCell
      integer :: i, iVertex, cellIndOnVertex

      vVortCell = 0.0
      do i = 1,nVerticesOnCell
         iVertex = verticesOnCell(i,c0)
         cellIndOnVertex = elementIndexInArray(c0, cellsOnVertex(:,iVertex), 3)
         vVortCell = vVortCell + kiteAreasOnVertex(cellIndOnVertex, iVertex)*vVortVertex(level, iVertex)/areaCell
      end do

      calc_verticalVorticity_cell = vVortCell
   end function calc_verticalVorticity_cell
'''

if __name__ == '__main__':
  #example_segment()
  
  ncfname = '/arctic1/nick/cases/v1.0/x4/longer/x4.kf.output.2006-08-15_00.00.00.nc'
  data = netCDF4.Dataset(ncfname,'r')
  
  tInd = 12
  
  cell2Site = segment_high_low_watershed(data, tInd)
  
  nCells = len(data.dimensions['nCells'])
  latCell = data.variables['latCell'][:];
  lonCell = data.variables['lonCell'][:];

  latCell *= 180./np.pi
  lonCell *= 180./np.pi
  
  isSite = cell2Site==range(nCells) #site goes to self
  plotSegment(latCell, lonCell, cell2Site, isSite)