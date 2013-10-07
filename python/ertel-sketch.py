#------------- This below is to convert into fortran --------------------------
#epv = np.dot(3dAbsoluteVorticity,gradTheta)/density   *(10^6) to scale to proper units

#This is rewritten to be clear so it may not even run in Python.
#A reasonable place for this to live would be in the diagnostics near pv_cell
#in core_atmosphere/dynamics/mpas_atm_time_integration.F
#Note that python is 0-based indexing (all of the cell,edge,... indices get decremented)

#We'll use the local xyz system for each cell since it's well aligned with physical gradients.
#x,y are the tangent plane and z is the radially out.

#Full curl is (dw/dy-dv/dz)i + (du/dz-dw/dx)j + (dv/dx-du/dy)k
#k component is solved for already through circulation to vertices. We'll use
#the kite_areas weighting as in pv_cell to get it at the cell centers.

#The other components of vorticity are:
#(dw/dy-dv/dz)i + (du/dz-dw/dx)j where velocities are in the local coordinate system.
#We'll also need dTheta/dX, dTheta/dY, dTheta/dZ.
#For all these, it's finite volume gradients in the
#horizontal and finite differencing in the vertical.


#imports
import numpy as np
import netCDF4

def formErtelPV(gradxu, gradtheta, density, unitX, unitY, unitZ):
  #Return epv in pvu's for a given cell
  omeg_e = (2.*np.pi) / (24*3600); #maybe keeping denom as integers reduces roundoff
  eVort = 2.*omeg_e #earth vorticity in [0,0,1] direction
  eVortDir = np.array([0,0,1])
  eVortComponents = np.array([eVort*np.dot(eVortDir,unitX), eVort*np.dot(eVortDir,unitY),eVort*np.dot(eVortDir,unitZ)])

  gradxu += eVortComponents #elementwise add for arrays
  epv = np.dot(gradxu,gradtheta)/density
  epv *= 1.e6 #in potential vorticity units, PVU is just a scaling of the SI units
  return epv

def coordinateSystem_cell(data, c0):
  #Local xyz system for cell is well aligned with physical gradients
  #Return 3x3 with x as [0,:], z as [2,:]...
  xyz = np.empty((3,3), dtype='float')
  xyz[0:2,:] = data.variables['cellTangentPlane'][c0,:,:] #cell, direction, 3
  xyz[2,:] = data.variables['localVerticalUnitVectors'][c0,:] #cell, 3
  #already checked that cross(x,y) = z
  return xyz

def form_faceVelocity(u, normalFace):
  #for c-grid staggering, have normal velocities at faces.
  return u*normalFace

def fluxSign(c0, cellsOnEdge):
  #normal should point out of cell for finite volume gradients.
  #connectivity info looks like normals point cellOnEdge[0]->cellsOnEdge[1], right?
  if (c0==cellsOnEdge[0]):
    return 1
  elif (c0==cellsOnEdge[1]):
    return -1
  else:
    print "\nUhoh. Cell {0} not on edge\n".format(c0)

#def calc_horizDeriv_fv(valEdges, nEdges, dvEdge, dhEdge, normalEdge, unitDeriv, volumeCell):
def calc_horizDeriv_fv(val0, valNbrs, nEdges, dvEdge, dhEdge, normalEdge, unitDeriv, volumeCell):
  #Component version of divergence theorem: SSS dF/dxi dV = SS F.ni dA
  #Assume dvEdges don't change between levels since curvature of earth is small.
  #Don't need top and bottom faces since no component of area in direction of deriv.
  #The result is dVal/dUnitDeriv

  #Inputs are cell center values at center cell and nbrs ordered by edge, length
  #and height of faces, and normal directions for face and derivative
  #normalEdges should point out of cell

  sum = 0.0;
  for e in xrange(nEdges):
    vale = .5*(val0+valNbrs[e])
    areaFace = dvEdge[e]*dhEdge[e] #treat as a rectangle
    areaFace *= np.dot(normalEdge[e,:],unitDeriv) #area of face in derivative direction
    sum += vale*areaFace
  sum /= volumeCell

  return sum

def calc_horizVorticity_wComponents(w0, wNbrs, unitX, unitY, nEdges, dvEdge, dhEdge, normalEdge, volumeCell):
  #w should be averaged from the faces to the cell centers.
  #should volumeCell = areaCell*hCell or areaCell*(mean(dhEdges))? dhEdge = z height of face on edge
  #normalEdges point out of cell

  #w
  dw_dx = calc_horizDeriv_fv(w0, wNbrs, nEdges, dvEdge, dhEdge, normalEdge, unitX, volumeCell)
  dw_dy = calc_horizDeriv_fv(w0, wNbrs, nEdges, dvEdge, dhEdge, normalEdge, unitY, volumeCell)

def driver_verticalDeriv_velocity(c0, data, unitVelocity, time0, level0,nLevels):
  #for vertical derivatives, we'll just use centered differences
  #and 1-sided for top and bottom cells to get estimates of derivatives at midpts
  #of faces. Then we'll average those by face area to estimate cell center value.

  #unitVelocity is the direction we'll specify, ie east/north/...

  nNbrs = data.variables['nEdgesOnCell'][c0]
  nbrs = data.variables['cellsOnCell'][c0,0:nNbrs]-1
  edges = data.variables['edgesOnCell'][c0,0:nNbrs]-1 #since python is 0-indexed
  lenEdges = data.variables['dvEdge'][edges]
  edgeNormals = data.variables['edgeNormalVectors'][edges,:]

  hNbrs = zgrid[nbrs,level0+1]-zgrid[nbrs,level0]

  u0 = data.variables['u'][time0,edges,level0]
  um=np.empty(nNbrs); up = um
  if (level0>0):
    um = data.variables['u'][time0,edges,level0-1]
  if (level0<nLevels-1):
    up = data.variables['u'][time0,edges,level0+1]

  #project into the direction of the coordinate system we're looking at
  for e in xrange(nNbrs):
    factor = np.dot(unitVelocity,edgeNormals[e,:])
    u0[e]*=factor; um[e]*=factor; up[e]*=factor

  #vertical finite difference based on zgrid.
  #heights of the edges are solved for in pv_cell diagnostic in mpas code
  dVel_dz = np.empty(nNbrs)
  if (level0==0):
    #1 sided to above
    #height is the distance between midpts of faces on that edge

    h0 = .5*(data.variables['zgrid'][c0,level0+1]+data.variables['zgrid'][c0,level0]) #level
    h0p = .5*(data.variables['zgrid'][c0,level0+2]+data.variables['zgrid'][c0,level0+1]) #level+1 (ie plus)
    for e in xrange(nNbrs):
      nbr = nrbs[e]
      hNbr = .5*(data.variables['zgrid'][hNbr,level0+1]+data.variables['zgrid'][hNbr,level0])
      hNbrp = .5*(data.variables['zgrid'][hNbr,level0+2]+data.variables['zgrid'][hNbr,level0+1])
      he = .5*(h0p+hNbrp)-.5*(h0+hNbr); #midpt+ - midpt-
      dVel_dz[e] = (up[e]-u0[e])/he
  elif(level==nLevels-1):
    #1 sided to below
    same procedure as above but to cell below

  else:
    #centered difference: (u_+1/2 - u_-1/2)/h where u_+1/2 = .5*(u_+1+u_0)
    #So, center value doesn't matter
    for e in xrange(nNbrs):
      get h0m, h0p, hNbrm, hNbrp #cell centers above and below this level
      hem = .5(h0m+hNbrm)
      hep = .5*(h0p+hNbrp)
      he = hep-hem
      dVel_dz[e] = (up[e]-um[e])/he
      
  #average face derivatives to cell center.
  #the question is to weigh by face area, inverse distance, mass (eg area*density on face)...
  sum = 0.0; sumArea = 0.0;
  for e in xrange(nNbrs):
    faceArea = lenEdge[e]*h_edge[e] #h_edge in pv_cell in fortran code
    sumArea+= faceArea
    sum += dVel_dz[e]*faceArea
  sum /= sumArea
  return sum

def driver_horizDeriv_cell_scalar(c0, data, var, unitDeriv, time0, level0):
  #call with var='theta' to get horizontal theta derivatives

  val0 = data.variables[var][time0, c0, level0] #for w, avg levels to cell center
  nNbrs = data.variables['nEdgesOnCell'][c0]
  nbrs = data.variables['cellsOnCell'][c0,0:nNbrs]-1
  valNbrs = data.variables[var][time0,nbrs,level0]

  h0 = zgrid[c0,level0+1]-zgrid[c0,level0]
  hNbrs = zgrid[nbrs,level0+1]-zgrid[nbrs,level0]
  #make hNbrs in height of edges
  hNbrs = .5*(hNbrs+h0)

  edges = data.variables['edgesOnCell'][c0,0:nNbrs]-1
  lenEdges = data.variables['dvEdge'][edges]
  edgeNormals = data.variables['edgeNormalVectors'][edges,:]
  cellsOnEdge = data.variables['cellsOnEdge'][edges,:]-1
  for e in xrange(nNbrs):
    edgeNormals[e] *= fluxSign(c0, cellsOnEdge[e,:])

  areaCell = data.variables['areaCell'][c0]
  volumeCell = h0*areaCell

  dval_ddir = calc_horizDeriv_fv(val0, valNbrs, nNbrs, lenEdges, hNbrs, edgeNormals, unitDeriv, volumeCell)
  return dval_ddir


#Use these to get dTheta/dz
def calc_vertDeriv_center(val0, valp, valm, htCell):
  #3 point finite difference where we first average values to the faces.
  #alternatives are least squares, fit parabola,...
  #valp is the value in the positive direction from the central value
  fp = .5*(val0+valp); fm = .5*(val0+valm)
  df = fp-fm;
  df_dx = df/htCell
  return df_dx

def calc_vertDeriv_one(valp, valm, dz):
  #2 point finite difference (for boundaries, may not have external nbr)
  df = valp-valm
  df_dx = df/dz
  return df_dx




#Old way here, don't use ----------------------------------------------
def calc_horizVorticity(vel0, velp, velm, velNbrs, unitX, unitY, unitZ, nEdges, dvEdge, dhEdge, normalEdge, volumeCell,
                        level, nLevels, dhCell):
  #velocities come in as uReconstructXYZ+w.nxyz
  #return i,j components of curl: (dw/dy-dv/dz)i + (du/dz-dw/dx)j + (dv/dx-du/dy)k
  #since k comes from averaged circulations on triangular dual

  #for inputs, dhCell should be (1)if middle, ht of cell (2)if bottom or top, ht between cells
  #volumeCell = areaCell*dzwCell
  #normalEdges point out of cell

  #w
  val0 = np.dot(vel0,unitZ)
  valNbrs = np.empty(nEdges,dtype='float')
  for e in xrange(nEdges):
    valNbrs[e] = np.dot(velNbrs,unitZ)
  dw_dx = calc_horizDeriv_fv(val0, valNbrs, nEdges, dvEdge, dhEdge, normalEdge, unitX, volumeCell)
  dw_dy = calc_horizDeriv_fv(val0, valNbrs, nEdges, dvEdge, dhEdge, normalEdge, unitY, volumeCell)

  #v
  dv_dz=0.0
  val0 = np.dot(vel0,unitY)
  if (level==0):
    valp = np.dot(velp,unitY)
    dv_dz = calc_vertDeriv_one(valp, val0, dhCell)
  elif (level==nLevels-1):
    valm = np.dot(velm,unitY)
    dv_dz = calc_vertDeriv_one(val0, valm, dhCell)
  else:
    valp = np.dot(velp,unitY)
    valm = np.dot(velm,unitY)
    dv_dz = calc_vertDeriv_center(val0, valp, valm, dhCell)

  #u
  du_dz=0.0
  val0 = np.dot(vel0,unitX)
  if (level==0):
    valp = np.dot(velp,unitX)
    du_dz = calc_vertDeriv_one(valp, val0, dhCell)
  elif (level==nLevels-1):
    valm = np.dot(velm,unitX)
    du_dz = calc_vertDeriv_one(val0, valm, dhCell)
  else:
    valp = np.dot(velp,unitX)
    valm = np.dot(velm,unitX)
    du_dz = calc_vertDeriv_center(val0, valp, valm, dhCell)

  #(dw/dy-dv/dz)i + (du/dz-dw/dx)j
  curli = dw_dy-dv_dz
  curlj = du_dz-dw_dx
  return (curli,curlj)