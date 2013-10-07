# -*- coding: utf-8 -*-
'''
Remember that the MPAS .nc file comes from a 1-based indexing system (Fortran).
Python and vtk are 0-based indexing.
Using the mesh connectivities as indices in Python is a gotcha!!!!!!!
'''

# imports
import netCDF4
import numpy as np

import vars_column as vars

def mpasInd2PythonInd(i):
  #mpas is 1-based for, say, the cell or vertex index.
  #be careful if cell or vertex reordering starts happening
  return i-1

def hitCycle(index, nn):
  #return true if have cycled.
  #we'll use this to print strings to file since i'm a little weary of having billions of numbers fit into a string
  (d,m) = divmod(index,nn)
  if (m>0 or index==0):
    return False
  return True

def print_xtime(data, tInd):
  s = ''.join(data.variables['xtime'][tInd][0:13]) #prints as '2006-07-24_12'
  print s

def example_nc2vtk():
  '''
  given .nc (MPAS NETCDF) file, output data into classic vtk format.
  This way, we can read into VisIt
  '''
  
  #file properties
  ncfname = '/home/nickszap/research/mpas/output.2010-10-23_00:00:00.nc' #input file
  vtkfname = 'polyTest.vtk' #output file
  
  data = open_netcdf_data(ncfname)
  
  #open the output vtk file and write header.
  fvtk = write_vtk_header_polydata(vtkfname, ncfname)
  
  #write nodes and cells
  nNodes =  write_vtk_xyzNodes(fvtk, data)
  nCells = write_vtk_polygons(fvtk, data)
  vLevels =  len(data.dimensions['nVertLevels'])
  
  #write some cell data
  fvtk.write('\nCELL_DATA '+str(nCells)+'\n')
  write_vtk_cellLandType(fvtk, data, nCells) 
  time = 0
  vLevel = 0
  write_vtk_cellCenterVelocity(fvtk, data, time,vLevel, nCells)
  write_vtk_pressureHeights(fvtk, data, nCells, time, vLevels, 50000) #500mb = 50,000Pa
  
  #write some node data
  
  #close the files
  data.close()
  fvtk.close()
#end driver

def example_nc2vtk_domain():
  '''
  given .nc (MPAS NETCDF) file, output data into classic vtk format.
  This way, we can read into VisIt
  '''

  import conn

  #file properties
  ncfname = '/home/nickszap/research/mpas/output.2010-10-23_00:00:00.nc' #input file
  vtkfname = 'domainTest5.vtk' #output file

  data = open_netcdf_data(ncfname)

  #partition mesh
  seed0 = 0; nSeeds = 10
  nCellsTotal = len(data.dimensions['nCells']); nVerticesTotal = len(data.dimensions['nVertices']);
  nEdgesOnCell = data.variables['nEdgesOnCell'][:];
  cellsOnCell = data.variables['cellsOnCell'][:]-1;
  cell2Site,seeds = conn.partition_max(seed0, cellsOnCell, nEdgesOnCell, nCellsTotal, nSeeds)

  #let's pick out a specific domain
  domainInd = 5; countThreshold=1
  cells = np.array(xrange(nCellsTotal))[cell2Site==seeds[domainInd]]
  nCells=len(cells)
  vOnCell = data.variables['verticesOnCell'][cells,:]-1

  verts = conn.gatherVerticesInRegion(nCells, vOnCell,
                                      nEdgesOnCell[cells,:], nVerticesTotal, countThreshold)

  #open the output vtk file and write header.
  fvtk = write_vtk_header_polydata(vtkfname, ncfname)

  #write nodes and cells
  nNodes = write_vtk_xyzNodes_domain(fvtk, data, verts)
  #need local indices for polygon vertices
  g2lVertex = conn.make_global2localMap(verts, [], nVerticesTotal)
  vOnCell = conn.make_localDomainNbrs(nCells, vOnCell, nEdgesOnCell[cells,:], g2lVertex)
  write_vtk_polygons_domain(fvtk, nCells, vOnCell, nEdgesOnCell[cells,:])

  #write some cell data
  fvtk.write('\nCELL_DATA '+str(nCells)+'\n')
  write_vtk_cellLandType_domain(fvtk, data, cells)

  #write some node data

  #close the files
  data.close()
  fvtk.close()
#end driver

def write_vtk_polyHorizConn_domain(data, fvtk, cells, nEdgesOnCell,nVerticesTotal):
  #write the polygons in this domain to file in legacy vtk format
  
  import conn
  
  #write nodes and cells
  nCells = len(cells)
  vOnCell = data.variables['verticesOnCell'][cells,:]-1
  verts = conn.gatherVerticesInRegion(nCells, vOnCell,nEdgesOnCell[cells,:], nVerticesTotal, 1)
  nNodes = write_vtk_xyzNodes_domain(fvtk, data, verts)

  #need local indices for polygon vertices
  g2lVertex = conn.make_global2localMap(verts, [], nVerticesTotal)
  vOnCell = conn.make_localDomainNbrs(nCells, vOnCell, nEdgesOnCell[cells,:], g2lVertex)
  write_vtk_polygons_domain(fvtk, nCells, vOnCell, nEdgesOnCell[cells,:])
  #done with header and conn for domain mesh-----------------------

def open_netcdf_data(fpath):
  try:
    data = netCDF4.Dataset(fpath,'r') #could try/catch but this kicks out some error message if no file
  except IOError as e:
    print "I/O error({0}): {1}".format(e.errno, e.strerror)
    pass
  print 'File information of: ', fpath
  #print data, data.variables #i think this gives all the keys to the data in the file
  print 'nTimes, nCells, nVertices, nVertLevels: ', len(data.dimensions['Time']),len(data.dimensions['nCells']), len(data.dimensions['nVertices']), len(data.dimensions['nVertLevels']) #can also continue line with \
  
  #do we want to change connectivites to 0-based indexing here????
  #verticesOnCell, cellsOnCell,...
  
  #print out the available vars (data structure is a dictionary)
  #mydict=data.variables
  #for key, value in sorted(mydict.iteritems(), key=lambda (k,v): (v,k)):
    #print "%s: %s" % (key, value)

  return data
#end

def write_vtk_xyzNodes(f, data):
  '''
  write out the xyz coordinates of all of the vertices.
  return the number of nodes
  '''
  
  #node locations
  xv = data.variables['xVertex'][:]
  yv = data.variables['yVertex'][:]
  zv = data.variables['zVertex'][:]
  
  dims = xv.shape
  nNodes = dims[0]
  
  #string = 'POINTS %s float\n' %nNodes
  string  = 'POINTS '+str(nNodes)+' float\n'
  #import string as s...string.format()
  f.write(string)
  for i in xrange(nNodes):
    #string = '%s %s %s\n'%(float(xv[i]), float(yv[i]), float(zv[i]))
    string = str(float(xv[i]))+' '+str(float(yv[i]))+' '+str(float(zv[i]))+'\n'
    f.write(string)
  #xyz vertices written
  
  return nNodes
#end

def write_vtk_xyzNodes_domain(f, data, verts):
  '''
  write out the xyz coordinates of all of the vertices in the domain.
  return the number of nodes
  '''

  #node locations
  xv = data.variables['xVertex'][verts]
  yv = data.variables['yVertex'][verts]
  zv = data.variables['zVertex'][verts]

  dims = xv.shape
  nNodes = dims[0]

  #string = 'POINTS %s float\n' %nNodes
  s  = 'POINTS '+str(nNodes)+' float\n'
  for i in xrange(nNodes):
    s += str(float(xv[i]))+' '+str(float(yv[i]))+' '+str(float(zv[i]))+'\n'
    if(hitCycle(i, 1000)):
      f.write(s); s=''
  f.write(s)
  #xyz vertices written

  return nNodes
#end

def formStringFromArray(vals, nVals):
  string = ''
  for i in xrange(nVals):
    string += str(vals[i])+' ' #creates a trailing space. problem?
  return string
#end

def test_cell():
  '''
  the base mesh is composed of layers of voronoi cells stacked radially out.
  here, we're seeing what info verticesOnCell contains...'2D' or 3d cells
  '''
  
  ncfname = '/home/nickszap/mpas/output.2010-10-23_00:00:00.nc' 
  data = open_netcdf_data(ncfname)
  
  vOnCell = data.variables['verticesOnCell'][:]
  dims = vOnCell.shape
  nCells=dims[0]
  maxEdges=dims[1]
  
  print "nCells and maxEdges: ",nCells,maxEdges
  for cell in xrange(4090):
    nv = verticesInCell(vOnCell[cell,:], maxEdges)
    print vOnCell[cell,0:nv]
#end  

def count_sizeVTKPolygons(nCells, eOnCell):
  '''
  To write the vtk file, we have to count the total amount of data for the polygons first.
  This is nPoly+nVerticesOnAllPolys
  '''
  count = nCells
  #count = nNodes
  for cell in xrange(nCells):
    nv = eOnCell[cell]
    count = count+nv
  return count
#end

def write_vtk_polygons(f, data):
  '''
  write the polygons to file. how are the cells outside of the surface level defined?
  return the # of cells
  '''
  vOnCell = data.variables['verticesOnCell'][:]-1
  eOnCell = data.variables['nEdgesOnCell'][:] #edges on each cell = num vertices on each cell
  dims = vOnCell.shape
  nCells=dims[0]
  maxEdges=dims[1]
  print "nCells and maxEdges: ",nCells,maxEdges
  f.write('\nPOLYGONS '+str(nCells)+' '+str(count_sizeVTKPolygons(nCells, eOnCell))+'\n')
  
  #vtk is 0 based but .nc has nodes as 1 based
  #vOnCell -= 1
  for cell in xrange(nCells):
    nv = eOnCell[cell]
    s = formStringFromArray(vOnCell[cell,:], nv)
    f.write(str(nv)+' '+s+'\n') #7 is cell type of polygon if change to unstructured
  return nCells
#end

def write_vtk_polygons_domain(f, nCells, verticesOnCell, nEdgesOnCell):
  '''
  write the polygons to file.
  use local vOnCell and nEdgeOnCell
  '''

  sz = count_sizeVTKPolygons(nCells, nEdgesOnCell)
  s = '\nPOLYGONS '+str(nCells)+' '+str(sz)+'\n'
  for cell in xrange(nCells):
    nv = nEdgesOnCell[cell]
    poly = formStringFromArray(verticesOnCell[cell,:], nv)
    s += str(nv)+' '+poly+'\n'
    if(hitCycle(cell, 1000)):
      f.write(s); s=''
  f.write(s)
#end

def write_vtk_header_polydata(vtkName, ncName):
  try:
    f = open(vtkName, 'w')
  except IOError as e:
    print "I/O error({0}): {1}".format(e.errno, e.strerror)
    pass

  #Print Header
  s = "# vtk DataFile Version 3.0\n"
  s+= 'Cell data from' + ncName+'\n'
  s+= "ASCII\n"
  s+= "DATASET POLYDATA\n"
  f.write(s)
  
  return f
#end

def write_levelData_float(info, f,vals,nCells):
  #write variable data on a level where field name is in info (w/o spaces i believe)
  #Interestingly enough, the writes may not show up in the file "immediately" since the
  #file may be buffered!!!
  s = '\nSCALARS '+info+' float 1\nLOOKUP_TABLE default\n'
  f.write(s)
  #
  s = ''
  for i in xrange(nCells):
    '''
    if (np.isnan(vals[i])): #writing nan to ascii .vtk seems not to work well
      vals[i] = missingVal
    '''
    s += str(float(vals[i]))+'\n'
    if(hitCycle(i, 1000)):
      f.write(s); s=''
  f.write(s)
#

def write_vtk_cellLandType(f, data, nCells):
  #string = '\nCELL_DATA '+str(nCells)+'\n'
  #f.write(string)
  f.write('\nSCALARS waterLand int 1\nLOOKUP_TABLE default\n')
  wl = data.variables['landmask'][:]
  for i in xrange(nCells):
    string  = str(wl[i])+'\n'
    f.write(string)

def write_vtk_cellLandType_domain(f, data, cells):
  #string = '\nCELL_DATA '+str(nCells)+'\n'
  #f.write(string)
  #f.write('\nSCALARS waterLand int 1\nLOOKUP_TABLE default\n')
  #wl = data.variables['landmask'][cells]
  f.write('\nSCALARS ter float 1\nLOOKUP_TABLE default\n')
  wl = data.variables['ter'][cells]
  for i in xrange(len(cells)):
    s  = str(wl[i])+'\n'
    f.write(s)

def write_vtk_cellCenterVelocity(f, data, time, vLevel, nCells):
  #write some velocity
  
  #string = '\nCELL_DATA '+str(nCells)+'\n'
  #f.write(string)
  label = 'cellVelocity_t'+str(time)+'_level'+str(vLevel)
  #f.write('\nVECTORS '+label+' float\nLOOKUP_TABLE default\n')
  f.write('\nVECTORS '+label+' float\n')
  
  uc = data.variables['uReconstructX'][time,:,vLevel]
  vc = data.variables['uReconstructY'][time,:,vLevel]
  wc = data.variables['uReconstructZ'][time,:,vLevel]
  
  for i in xrange(nCells):
    string = '%s %s %s\n'%(float(uc[i]), float(vc[i]), float(wc[i]))
    f.write(string)

def write_vtk_pressureHeights_eqnState(f, data, nCells, time, nLevels, pLevel):
  #IC files don't have pressure_p
  
  f.write('\nSCALARS z_p'+str(pLevel)+' float 1\nLOOKUP_TABLE default\n')
  s = ''
  for hcell in xrange(nCells):
    ht = vars.calc_htCellCenter_column(data,hcell,nLevels)
    p = vars.calcPressure_column(data,hcell,time, nLevels)
    
    (l,dl) = calcIndexOfValue(pLevel, p, nLevels)
    z = calcValueOfIndex(l,dl,ht)
    
    s += str(float(z))+'\n'
    if(hitCycle(hcell,1000)):
      f.write(s);
      s=''
  f.write(s)
    
  
def write_vtk_pressureHeights(f, data, nCells, time, vLevels, pLevel):
  zInterface = data.variables['zgrid'][:]
  pressP = data.variables['pressure_p'][time,:,:] #perturbation from base pressure
  pressBase = data.variables['pressure_base'][time,:,:]
  pressure = [None]*vLevels #just calculate column at a time
  zMid = [None]*vLevels #heights of midpoints of cells
  
  f.write('\nSCALARS z_p'+str(pLevel)+' float 1\nLOOKUP_TABLE default\n')
  for cell in xrange(nCells):
    #calc pressure and cell heights at cell centers
    for l in xrange(vLevels):
      pressure[l] = pressP[cell,l]+pressBase[cell,l]
    #
    for l in xrange(vLevels):
      zMid[l] = .5*(zInterface[cell,l]+zInterface[cell,l+1])
    #
    
    (l,dl) = calcIndexOfValue(pLevel, pressure, vLevels)
    z = calcValueOfIndex(l,dl,zMid)
    
    string  = str(float(z))+'\n'
    f.write(string)

missingVal = -99999
def calcIndexOfValue(val0,vals, nVals):
  #do a linear interpolation to trap val0 in a column by going top to surface.
  #return,say, (1,.5) if val is halfway between ind1 and ind2.
  #will break by returning (missingVal,string) instead if val0 isn't in (min(vals),max(vals))
  for l in xrange(nVals-1,0,-1): #[nVals,nVals-1,...,1]
    vH = vals[l]; vL = vals[l-1]; #high and low
    sgn1 = val0-vH; sgn2 = val0-vL;
    if (sgn1*sgn2<=0):
      #sandwiched value. equal in case val0 is a vals[l].
      #get linear interpolation: val0 = vals[l]+dvals/dl * dl
      #Avoid divide by 0 by just assuming value is halfway between...
      dv_dl = vH-vL;
      dl = 0;
      if (abs(dv_dl)<1.e-12):
        dl = .5;
      else:
        dl = (val0-vL)/dv_dl
      return (l-1,dl)
  
  s = "Uhoh. val0 not found in column\n"
  #print s
  return (missingVal,s)

def calcIndexOfValue_fromBottom(val0,vals, nVals):
  #do a linear interpolation to trap val0 in a column by going surface to top.
  #return,say, (1,.5) if val is halfway between ind1 and ind2.
  #will break by returning (missingVal,string) instead if val0 isn't in (min(vals),max(vals))
  for l in xrange(0,nVals-1,1): #[nVals,nVals-1,...,1]
    vH = vals[l+1]; vL = vals[l]; #high and low
    sgn1 = val0-vH; sgn2 = val0-vL;
    if (sgn1*sgn2<=0):
      #sandwiched value. equal in case val0 is a vals[l].
      #get linear interpolation: val0 = vals[l]+dvals/dl * dl
      #Avoid divide by 0 by just assuming value is halfway between...
      dv_dl = vH-vL;
      dl = 0;
      if (abs(dv_dl)<1.e-12):
        dl = .5;
      else:
        dl = (val0-vL)/dv_dl
      return (l,dl)

  s = "Uhoh. val0 not found in column\n"
  #print s
  return (missingVal,s)

def calcValueOfIndex(l,dl,vals):  
  #do a linear interpolation to find the value corresponding to, say, ind=1.5
  #val0 = vals[l]+dvals/dl * dl. return NaN if index not in column
  
  if (l==missingVal):
    #return float('nan') #looks like paraview doesn't support ascii nan's
    return l
  
  dv_dl = vals[l+1]-vals[l];
  val0 = vals[l]+dv_dl*dl;
  return val0
  
if __name__ == '__main__':
  driver_nc2vtk()

def write_vtk_var_timeLevelCells(f, var, data, level, time, nCells):
  #write single time,level data for all cells
  
  f.write('\nSCALARS '+var+' float 1\nLOOKUP_TABLE default\n')
  vals = data.variables[var][time,:,level]
  for cell in xrange(nCells):
    string  = str(float(vals[cell]))+'\n'
    f.write(string)
  #
  
def write_vtk_staticGeoFields(f,data,nCells):
  #write some surface fields of the cells in our scvt mesh that don't change over time
  
  #could use a dict if want titles of vars different from keys
  vars = ['ter','landmask','snoalb']
  for var in vars:
    f.write('\nSCALARS '+var+' float 1\nLOOKUP_TABLE default\n')
    vals = data.variables[var][:]
    for cell in xrange(nCells):
      string  = str(float(vals[cell]))+'\n'
      f.write(string)
#end

def write_vtk_timeGeoFields(f,data,nCells, nTimes):
  #write some surface fields of the cells in our scvt mesh that do change over time
  vars = ['vegfra','xice','sfc_albbck', 'skintemp', 'sst', 'tmn'] #use seaice if want xice binary
  for t in range(nTimes):
    for var in vars:
      f.write('\nSCALARS '+var+str(t)+' float 1\nLOOKUP_TABLE default\n')
      vals = data.variables[var][t,:]
      for cell in xrange(nCells):
        string  = str(float(vals[cell]))+'\n'
        f.write(string)
#end
  
def write_vtk_cellLandType(f, data, nCells):
  #string = '\nCELL_DATA '+str(nCells)+'\n'
  #f.write(string)
  f.write('\nSCALARS waterLand int 1\nLOOKUP_TABLE default\n')
  wl = data.variables['landmask'][:]
  for i in xrange(nCells):
    string  = str(wl[i])+'\n'
    f.write(string)

#------------------ 3d field of prisms ----------------------------
def example_prism_single(f):
  #write out a vtk file with a single to make sure format is proper for viz
  #Print Header
  f.write("# vtk DataFile Version 3.0\n")
  foutInfo = 'Cell data from test\n'
  f.write(foutInfo)
  f.write("ASCII\n")
  f.write("DATASET UNSTRUCTURED_GRID\n")
  
  numPts = 6
  pxyz = np.array([0.,0.,0., 1.,0.,0., 0.,1.,0.,
                  0.,0.,1., 1.,0.,1., 0.,1.,1.])
  s  = 'POINTS '+str(numPts)+' float\n'
  f.write(s)
  for p in xrange(numPts):
    s = str(float(pxyz[3*p+0]))+' '+str(float(pxyz[3*p+1]))+' '+str(float(pxyz[3*p+2]))+'\n'
    f.write(s)
  
  nPrisms = 1
  s = '\nCELLS '+str(nPrisms)+' '+str(nPrisms*(6+1))+'\n' #6 pts per prism +1 int for #conn, ie '6 vertex0 v1...v5' for prism
  f.write(s)
  for p in xrange(nPrisms):
    s = '6'
    for i in xrange(6): #
      s+=' '+str(i)
    s+='\n'
    f.write(s)
  #
  s = '\nCELL_TYPES '+str(nPrisms)+'\n'
  for i in range(nPrisms):
    s+= '13\n'
  f.write(s)

def example_prism_column():
  #write out a column of the mesh as prisms.
  #the column is a triangle radially out identified by primal (eg hexagonal) vertex
  import vars_column
  
  ncfname = '/home/nickszap/research/mpas/output.2010-10-23_00:00:00.nc' #input file
  vtkfname = 'prismTest.vtk' #output file
  data = open_netcdf_data(ncfname)
  nCells = len(data.dimensions['nCells']); nLevels = len(data.dimensions['nVertLevels'])
  nVertices = len(data.dimensions['nVertices'])
  
  f = write_vtk_header_unstructured(vtkfname, ncfname)
  
  #hack for this case ---------------
  nVertices = 1; vertexId = 10
  nCells = 3; cells = data.variables['cellsOnVertex'][vertexId]-1 #0-indexing
  xyzc = vars_column.calc_cellCenterColumn(data,cells,nLevels) #index as xyzc[celli,level,dirs]
  #end hack ----------------
  
  numPts = nCells*nLevels
  s  = 'POINTS '+str(numPts)+' float\n'
  f.write(s)
  for c in xrange(nCells):
    for l in xrange(nLevels):
      s = str(float(xyzc[c,l,0]))+' '+str(float(xyzc[c,l,1]))+' '+str(float(xyzc[c,l,2]))+'\n'
      f.write(s)
  #
  nPrisms = nVertices*(nLevels-1) #1 triangle per vertex with prisms centered over interfaces
  s = '\nCELLS '+str(nPrisms)+' '+str(nPrisms*(6+1))+'\n' #6 pts per prism +1 int for #conn, ie '6 vertex0 v1...v5' for prism
  f.write(s)
  for v in xrange(nVertices):
    for l in xrange(nLevels-1):
      s = '6'
      for i in xrange(3): #visit 3 cells on vertex at base. unchecked winding
        c = cells[i] #would be cells[v,i] for multiple vertices
        ind = cellToTriangleInd(i,l, nLevels) #would be cellToTriangleInd(c,l, nLevels)
        s+= ' '+str(ind)
      for i in xrange(3): #visit 3 cells on vertex at top. unchecked winding
        c = cells[i] #would be cells[v,i] for multiple vertices
        ind = cellToTriangleInd(i,l+1, nLevels) #would be cellToTriangleInd(c,l, nLevels)
        s+= ' '+str(ind)
      s+='\n'
      f.write(s)
  #
  s = '\nCELL_TYPES '+str(nPrisms)+'\n'
  for i in range(nPrisms):
    s+= '13\n'
  f.write(s)
  
  data.close()
  f.close()
  
def example_prismsRegion():
  import vars_column
  import conn
  ncfname = '/home/nickszap/research/mpas/output.2010-10-23_00:00:00.nc' #input file
  vtkfname = 'prismTestRegion.vtk' #output file
  data = open_netcdf_data(ncfname)
  nCellsTotal = len(data.dimensions['nCells']); nLevels = len(data.dimensions['nVertLevels'])
  nVerticesTotal = len(data.dimensions['nVertices'])

  f = write_vtk_header_unstructured(vtkfname, ncfname)

  #hack for this case ---------------
  c0 = 10 #do cell and nbrs
  nEdgesOnCell = data.variables['nEdgesOnCell'][:]
  nNbrs = nEdgesOnCell[c0]
  nbrs = data.variables['cellsOnCell'][c0,0:nNbrs]-1
  cells = [c0]+nbrs.tolist()

  verticesOnCell = data.variables['verticesOnCell'][cells,:]-1
  vertices = conn.gatherVerticesInRegion(len(cells), verticesOnCell, nEdgesOnCell[cells], nVerticesTotal, 3)
  nVertices = len(vertices);
  nCells = len(cells);
  xyzc = vars_column.calc_cellCenterColumn(data,cells,nLevels) #index as xyzc[celli,level,dirs]
  #end hack ----------------

  numPts = nCells*nLevels
  s  = 'POINTS '+str(numPts)+' float\n'
  f.write(s)
  for c in xrange(nCells):
    for l in xrange(nLevels):
      s = str(float(xyzc[c,l,0]))+' '+str(float(xyzc[c,l,1]))+' '+str(float(xyzc[c,l,2]))+'\n'
      f.write(s)
  #
  #for prisms need global cell index to local map.
  #only need map for cells in region so can possibly save a bit on memory
  g2lCell = np.zeros(np.max(cells),dtype=int) #global to local map for cells on horizontal
  for i,c in enumerate(cells):
    g2lCell[c]=i
  cellsOnVertex = data.variables['cellsOnVertex'][vertices,:]-1

  nPrisms = nVertices*(nLevels-1) #1 triangle per vertex with prisms centered over interfaces
  s = '\nCELLS '+str(nPrisms)+' '+str(nPrisms*(6+1))+'\n' #6 pts per prism +1 int for #conn, ie '6 vertex0 v1...v5' for prism
  f.write(s)
  for v in xrange(nVertices):
    for l in xrange(nLevels-1):
      s = '6'
      for i in xrange(3): #visit 3 cells on vertex at base. unchecked winding
        c = cellsOnVertex[v,i] #would be cells[v,i] for multiple vertices
        ind = g2lCell[c]
        ind = cellToTriangleInd(ind,l, nLevels) #would be cellToTriangleInd(c,l, nLevels)
        s+= ' '+str(ind)
      for i in xrange(3): #what happens if 4 equidistant cell centers. still 3???
        c = cellsOnVertex[v,i] #would be cells[v,i] for multiple vertices
        ind = g2lCell[c]
        ind = cellToTriangleInd(ind,l+1, nLevels) #would be cellToTriangleInd(c,l, nLevels)
        s+= ' '+str(ind)
      s+='\n'
      f.write(s)
  #
  s = '\nCELL_TYPES '+str(nPrisms)+'\n'
  for i in range(nPrisms):
    s+= '13\n'
  f.write(s)

  data.close()
  f.close()

def example_plot3d_epv_arctic():
  import vars_column
  import conn
  #ncfname = '/home/nickszap/research/mpas/output.2010-10-23_00:00:00.nc' #input file
  #ncfname = '/arctic1/nick/cases/cfsr/output/x4.cfsr.output.2006-08-01_00.00.00.nc'
  #ncfname = '/arctic1/nick/cases/v1.0/x4/longer/x4.kf.output.2006-08-15_00.00.00.nc'
  ncfname = '/arctic1/nick/cases/650k/x1.t.output.2006-08-15_00.00.00.nc'
  #vtkfname = '/data01/epv-prismTestArctic.vtk' #output file
  data = open_netcdf_data(ncfname)
  nCellsTotal = len(data.dimensions['nCells']); nLevels = len(data.dimensions['nVertLevels'])
  nVerticesTotal = len(data.dimensions['nVertices']); nTimes = len(data.dimensions['Time'])

  nEdgesOnCell = data.variables['nEdgesOnCell'][:]

  #hack for this case ---------------
  latCell = data.variables['latCell'][:]
  #latThresh = 85.*np.pi/180.
  latThresh = 75.*np.pi/180.
  cells = conn.gatherArcticCells(latCell, nCellsTotal, latThresh)
  #should get halo as well
  #conn.get_arcticHalo(cells, latCell, latThresh, cellsOnCell, nEdgesOnCell)

  verticesOnCell = data.variables['verticesOnCell'][cells,:]-1
  vertices = conn.gatherVerticesInRegion(len(cells), verticesOnCell, nEdgesOnCell[cells], nVerticesTotal, 3)
  nVertices = len(vertices);
  nCells = len(cells);
  print "Number of cells in region: ", nCells
  xyzc = vars_column.calc_cellCenterColumn(data,cells,nLevels) #index as xyzc[celli,level,dirs]
  #end hack ----------------
  
  for timeInd in xrange(nTimes):
    vtkfname = '/data01/epv-prismTestArctic_650k-'+str(timeInd)+'.vtk' #output file
    f = write_vtk_header_unstructured(vtkfname, ncfname)
  
    #Here we'll manipulate xyz for more clear viz
    #ince we're in shallow water land, the vertical is quite scrunched wrt horizontal.
    #We could do this anywhere over earth by taking a "central" cell in the region's
    #tangent plane as the reference.
    #We'll stretch it a bit to help with viz, since over arctic magnify z
    zFactor = 100.
    numPts = nCells*nLevels
    s  = 'POINTS '+str(numPts)+' float\n'
    f.write(s)
    zgrid = data.variables['zgrid'][cells,:]
    for c in xrange(nCells):
      s = ''
      for l in xrange(nLevels):
        #s += str(float(xyzc[c,l,0]))+' '+str(float(xyzc[c,l,1]))+' '+str(float(xyzc[c,l,2]))+'\n'
        s += str(float(xyzc[c,l,0]))+' '+str(float(xyzc[c,l,1]))+' '+str(float(zFactor*zgrid[c,l]))+'\n'
      f.write(s)
    #
    #for prisms need global cell index to local map.
    #only need map for cells in region so can possibly save a bit on memory
    g2lCell = np.zeros(np.max(cells)+1,dtype=int) #global to local map for cells on horizontal
    for i,c in enumerate(cells):
      g2lCell[c]=i
    cellsOnVertex = data.variables['cellsOnVertex'][vertices,:]-1

    nPrisms = nVertices*(nLevels-1) #1 triangle per vertex with prisms centered over interfaces
    s = '\nCELLS '+str(nPrisms)+' '+str(nPrisms*(6+1))+'\n' #6 pts per prism +1 int for #conn, ie '6 vertex0 v1...v5' for prism
    f.write(s)
    for v in xrange(nVertices):
      for l in xrange(nLevels-1):
        s = '6'
        for i in xrange(3): #visit 3 cells on vertex at base. unchecked winding
          c = cellsOnVertex[v,i] #would be cells[v,i] for multiple vertices
          ind = g2lCell[c]
          ind = cellToTriangleInd(ind,l, nLevels) #would be cellToTriangleInd(c,l, nLevels)
          s+= ' '+str(ind)
        for i in xrange(3): #what happens if 4 equidistant cell centers. still 3???
          c = cellsOnVertex[v,i] #would be cells[v,i] for multiple vertices
          ind = g2lCell[c]
          ind = cellToTriangleInd(ind,l+1, nLevels) #would be cellToTriangleInd(c,l, nLevels)
          s+= ' '+str(ind)
        s+='\n'
        f.write(s)
    #
    s = '\nCELL_TYPES '+str(nPrisms)+'\n'
    for i in xrange(nPrisms):
      s+= '13\n'
    f.write(s)

    f.write('\nPOINT_DATA '+str(numPts)+'\n')

    s = '\nSCALARS epv float 1\nLOOKUP_TABLE default\n'
    f.write(s)
    #timeInd=0
    ertel_pv = data.variables['ertel_pv'][timeInd,cells,:]
    for c in xrange(nCells):
      s = ''
      for l in xrange(nLevels):
        s += str(float(ertel_pv[c,l]))+'\n'
      f.write(s)

    f.close()
    
  data.close()

def example_smallRegion_3dplane(data, cells): #write this!!!!
  #write a 3d volume of a region. stretch vertical for easier viz
  #since ~30km atmosphere height is super thin wrt horizontal spacing
  #considering it's chopped up with nVertLevels
  #1 idea for use is showing evolution of >2pvu atmosphere.
  
  #plane vertical is average of vertical normals in region.
  #xyzp = xyz-xyz.*planeVertical removes component in direction of planeVertical.
  #xyzp += planeVertical*zgrid*factor adds in a stretched vertical.
  #a number of options for factor...say:
  # - make lowest level unit aspect ratio wrt cell width (sqrt(cellArea)?)
  pass

#-------------------section for writing prisms ---------------------------------
def write_vtk_header_unstructured(vtkName, ncName):
  try:
    f = open(vtkName, 'w')
  except IOError as e:
    print "I/O error({0}): {1}".format(e.errno, e.strerror)
    pass

  #Print Header
  f.write("# vtk DataFile Version 3.0\n")
  foutInfo = 'Cell data from' + ncName+'\n'
  f.write(foutInfo)
  f.write("ASCII\n")
  f.write("DATASET UNSTRUCTURED_GRID\n")
  
  return f
  
def write_vtk_headerConn_prisms(vtkName, ncName, data, nCells, nLevels):
  f = write_vtk_header_unstructured(vtkName, ncName)
  
  write_vtk_cellCenters(f, data, nCells, nLevels)
  write_vtk_prisms(f,data, data.variables['xCell'][:],data.variables['yCell'][:],data.variables['zCell'][:],nCells,nLevels)
  
  return f

def cellToTriangleInd(hCell,l, nLevels):
  #index cell center [hCell,level] flatly as hCell*nLevels+level
  #this has to match the way they're written as points.
  return hCell*nLevels+l

def write_vtk_cellCenters(f, data, nCells, nLevels):
  '''
  The cell centers form the vertices of the triangles.
  We index cell center [hCell,level] flatly as hCell*nLevels+level
  return the number of cells (triangle vertices)
  '''
  import sys; sys.path.append("/home/nickszap/Dropbox/pythonScripts/mem")
  import varsM
  
  xyz, ht = varsM.get_allCellCenters(data,nCells,nLevels)
  
  s  = 'POINTS '+str(nCells*nLevels)+' float\n'
  for hCell in xrange(nCells):
    for l in xrange(nLevels):
      s += str(float(xyz[hCell,l,0]))+' '+str(float(xyz[hCell,l,1]))+' '+str(float(xyz[hCell,l,2]))+'\n'
      if(hitCycle(hCell, 100)):
        f.write(s); s=''
  f.write(s)
  
  return (nCells*nLevels,xyz)
#end

def formTriangleString(inds, l, nLevels):
  s = str(cellToTriangleInd(inds[0],l,nLevels))+' '+str(cellToTriangleInd(inds[1],l,nLevels))+' '+str(cellToTriangleInd(inds[2],l,nLevels))
  return s

def write_vtk_prisms(f,data, x,y,z,nCells, nLevels):
  #vtk has a finite set of volume cells. we're "lucky" that the dual is prisms.
  #x,y,z are of the surface
  #each corner of the prism is a cell center
  
  import plotPython
  tri = plotPython.recoverTriangles(data, nCells) #surface triangulation
  nTri = len(tri)
  tri = plotPython.windTriangles(tri, nTri, x,y,z)
  
  nPrisms = nTri*(nLevels-1)
  s = '\nCELLS '+str(nPrisms)+' '+str(nPrisms*(6+1))+'\n' #6 pts per prism +1 int for #conn, ie '6 vertex0 v1...v5' for prism
  for t in xrange(nTri):
    nodes = tri[t,:] #surface nodes
    for l in xrange(nLevels-1):
      s += '6 ' #number of points
      s+= formTriangleString(nodes, l, nLevels)
      s+= ' '+formTriangleString(nodes, l+1, nLevels)+'\n'
      if(hitCycle(t, 1000)):
        f.write(s); s=''
  f.write(s)
  
  wedgeID = 13
  s = '\nCELL_TYPES '+str(nPrisms)+'\n'
  for i in range(nPrisms):
    s+= '13\n'
    if(hitCycle(i, 1000)):
      f.write(s); s=''
  f.write(s)

  
def modifyNetcdfField_example():
  fname = '/home/nickszap/research/mpas/init.nc'
  data = netCDF4.Dataset(fname,'a') #r+ same as a for append?
  temp = data.variables['cf1'][:]
  temp[0] = -999
  data.variables['cf1'][:] = temp[:] #data.variables['cf1'][0]=-999 also works
  data.close()
  
#old (and likely untested) functions ----------------------------------------------------
def calcTheta_m(theta, qv):
  #modified moist potential temperature given potential temperature and water vapor mixing ratio
  Rv = 461.6; 
  Rd = 287.04;
  return theta*(1.+qv*Rv/Rd)
def calcPressure_thermo(p0, dzeta_dz, theta_m):
  #Return the pressure in Pa at a given state.
  #from skamarock A Multi-scale Nonhydrostatic Atmospheric Model:
  #Using Centroidal Voronoi Tesselations and C-Grid Staggering
  #p = p0 (Rd dZeta/dz Thetam / p0)^(cp/cv)
  
  #Hacking through model fields, I get the following:
  #dZeta/dz is a vertical grid field stored at cell faces. To get at cell center,
  #average the lower and upper faces.
  #Is p0 standard atmosphere (101325Pa), 100 000 Pa, surface pressure, or the reference-state dry pressure????!!!!
  
  #!!!!!!!!!!!!!!!Check these with cavallo!!!!!!!!!
  from math import pow
  Cp = 1004.5;
  Cv = 717.5;
  Rd = 287.04;
  base = Rd*dzeta_dz*theta_m/p0
  power = Cp/Cv
  print '\nFor pressure,base^power: ', base, power
  p = p0*pow(base, power)
  return p
def pascalLevel2Height(pLevel, theta_m, dzeta_dz, zInterface, vLevels):
  #find the height corresponding to the pressure level in a column by scanning down from the top of the domain.
  #we'll do linear interpolation between the bounding cells,
  #but fancier ways exist in numpy.interp and scipy.interp1d with root finding (scipy.optimize.newton or scipy.optimize.brentq)
  #return (index of cell with p>pLevel, z height of pLevel)
  
  #calc pressures of cell centers in column
  p0 = 101325.0 #Who knows if this is the right value to use???
  
  p = [None]*vLevels
  for i in range(vLevels): #go from surface to top cell (p~0) stepping by 1
    p[i] = calcPressure(p0, .5*(dzeta_dz[i+1]+dzeta_dz[i]), theta_m[i]) #i believe the model averages faces to get cell values!!!
  
  #make sure pLevel falls in pressure range of column
  pMax = max(p); pMin = min(p);
  if ((pLevel<pMin) or (pLevel>pMax)):
    print "Uhoh. Pressure level of %s isn't bracketed by column range (%s,%s)" %(pLevel, pMin, pMax)
    return pMin*(pLevel<pMin)+pMax*(pLevel>pMax) #doesn't return (ind,z) like rest of function!!!
  
  dpOld = p[0]-pLevel
  for i in range(1,vLevels,1):
    #find where pLevel falls in column
    dp = p[i]-pLevel
    
    if (dp*dpOld>0):
      #haven't changed sign on dp
      dpOld = dp
    else :
      #we're on the other side of the pressure level. linear interpolation.
      zim1 = .5*(zInterface[i-1]+zInterface[i]) #height of lower cell
      zi = .5*(zInterface[i]+zInterface[i+1]) #height of higher cell
      dp_dz = (p[i]-p[i-1])/(zi-zim1)
      dp = pLevel-p[i-1]
      dz = dp/dp_dz
      return (i,zim1+dz)
    
    #won't return anything if don't hit else
def write_vtk_pressureHeightsOld(f, data, nCells, time, vLevels, pLevel):
  zInterface = data.variables['zgrid'][:]
  dzeta_dz = data.variables['zz'][:] #at cell centers! top value is junk
  
  #theta_m = data.variables['theta_m'][time,:,:] #Does theta_m not exist in the output file???????
  qv = data.variables['qv'][time,:,:]
  theta = data.variables['theta'][time,:,:]
  theta_m = [None]*vLevels #just calculate column at a time
  
  f.write('\nSCALARS height_'+str(pLevel)+' float 1\nLOOKUP_TABLE default\n')
  for cell in range(nCells):
    #calc theta_m at cell centers
    for l in range(vLevels):
      theta_m[l] = calcTheta_m(theta[cell,l], qv[cell,l])
    #calculated theta_m for column
    
    (ind, z) = pascalLevel2Height(pLevel, theta_m[:], dzeta_dz[cell,:], zInterface[cell,:], vLevels)
    string  = str(float(z))+'\n'
    f.write(string)
#end of old an untested --------------------------------------
