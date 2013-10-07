# Define where my python is located
#!/usr/local/epd/bin/python

#Don't have enough memory to load full meshes of variables. We'll exploit
#horizontal domain decomposition since nVertLevels is "small"

#use to check available RAM
#import memUtils
#memUtils.calcNDecomp(nVals)

#If we don't want to go through domain decomp, virtual memory would be no slower than file IO, right?

#imports
import numpy as np
import netCDF4
import time
import gc

import sys; sys.path.append("/home/nickszap/Dropbox/pythonScripts");
import conn
import output_data

#air,thermo,earth,... constants copied from one of Cavallo's codes
Cp = 1004.5;
Cv = 717.5;
gamma = Cp/Cv;
Rd = 287.04;
Rv = 461.6;
RvRd = Rv / Rd;
g = 9.81;
L = 2.50e6;
Talt = 288.1500;
Tfrez = 273.1500;
To = 300;
Po = 101325.;
Pr = 100000.;
lapsesta = 6.5 / 1000;
kappa = Rd / Cp;
epsil = Rd/Rv;
R_earth = 6371200.;
omeg_e = (2.*np.pi) / (24*3600); #maybe keeping denom as integers reduces roundoff
eo = 6.11;

def calc_theta_m(theta, qv):
  #modified moist potential temperature given potential temperature and water vapor mixing ratio
  return theta*(1.+qv*RvRd)
  
def calc_temperature(theta, p):
  p0 = 1.e5 #this is in the MPAS code
  fac = (p/p0)**kappa
  return theta*fac

def calc_coriolis(lat):
  #coriolis parameter for lat in radians
  f = 2.*omeg_e*np.sin(lat)
  return f

def calc_pressure_mpasEqnState(zz, theta_m, rho):
  #calculate the pressure given the equation of state
  p0 = 1.e5;
  #num = Rd*data.variables['zz'][hcell,level]*data.variables['theta_m'][time,hcell,level]
  num = Rd*zz*theta_m*rho
  fac = np.power(num/p0,gamma)
  return p0*fac
  
def calc_pressure_perturb(pressure_base, pressure_p):
  p = pressure_base+pressure_p
  return p

def calc_moistDensity(data, hcell, timeInd, level):
  #rho_moist = rho_dry(1+qv+qc+qr+...)
  #or should we compute a density consistent with ideal gas law from pressure and temperature?
  pass

'''
There are a number of approaches to calculating SLP (p at 0 meters). Standard atmosphere based
on a curve fit of a global climatological distribution was decribed here on 2/7/2013:
http://en.wikipedia.org/wiki/Atmospheric_pressure#Altitude_atmospheric_pressure_variation

However, putting the validity of the hydrostatic, adiabatic lapse rate,... assumptions aside,
our complication may be potential temperature being the model variable
'''

def calc_pressureHypsometric(T0, p0, dz):
  #dz is positive radially out
  scale = -g*dz/(Rd*T0)
  p = p0*np.exp(scale)
  return p

def calc_slp_surface(pSurf, zSurf, tSurf):
  #using surface or first cell values, estimate p at height=0m

  if (0):
    tol = 5. #tolerance to sea level in m
    if (abs(zSurf)<tol):
      #this "is" 0 meters, just to save some time and not lose accuracy
      return pSurf

  #move temperature to surface
  tRef = calc_temperatureFromLapse(tSurf, -.5*zSurf) #half way between surface and sea level
  slp = calc_pressureHypsometric(tRef, pSurf, -zSurf)
  return slp
  
def calc_temperatureFromLapse(t0, dz):
  #return the temperature at the changed height in meters using standard atmospheric lapse rate.
  #dz is positive up and t decreases
  return t0-lapsesta*dz

def calc_slp_ben(tSurf, hSurf, pSurf):
  
  t0 = calc_temperatureFromLapse(tSurf,-.5*hSurf) #"reference temperature" halfway between ground and 0m
  p = calc_pressureHypsometric(t0, pSurf, -hSurf); #dz is positive up
  return p
  
def calc_slp_hydro(ph, rho, dz):
  #Hydrostatic with fixed density
  #dz is positive up
  
  dp_dz = -rho*g;
  p0 = ph+dp_dz*dz #where dz=ht-0m= ht
  return p0

def form_curl(u_xyz, v_xyz, w_xyz):
  #form the vector grad x u given all of the spatial derivs
  
  vort = np.empty(3)
  vort[0] = w_xyz[1]-v_xyz[2]
  vort[1] = u_xyz[2]-w_xyz[0]
  vort[2] = v_xyz[0]-u_xyz[1]
  return vort

def form_curl_column(u_xyz, v_xyz, w_xyz, nLevels):
  #form the vector grad x u given all of the spatial derivs

  vort = np.empty(nLevels,3)
  for l in xrange(nLevels):
    vort[l,0] = w_xyz[l,1]-v_xyz[l,2]
    vort[l,1] = u_xyz[l,2]-w_xyz[l,0]
    vort[l,2] = v_xyz[l,0]-u_xyz[l,1]
  return vort

def calc_vorticitySCVT():
  #using the normal and tangential velocities on the mesh
  pass

def calc_absoluteVorticity(vort):
  #return the absolute vorticity at a given cell center
  vort[2] = vort[2]+2.*omeg_e
  return vort

def calc_ertelPV(vort, theta_xyz, rho):
  #Ertel PV is (wLocal+wEarth).grad(tracer)/rho
  #variables are dimensioned [level(, xyz)]
  #Opposite sign in the southern hemisphere
  
  #local to absolute vorticity
  vort[2] += 2.*omeg_e #absolute with earth

  #form epv
  epv = np.dot(vort,theta_xyz)/rho
  epv = epv*1e6 #to PVU
  return epv

def calc_ertelPV_column(vort, theta_xyz, rho, nLevels):
  #Ertel PV is (wLocal+wEarth).grad(tracer)/rho
  #variables are dimensioned [level(, xyz)]
  #Opposite sign in the southern hemisphere
  
  #local to absolute vorticity
  for l in xrange(nLevels):
    vort[l,2] += 2.*omeg_e #absolute with earth

  #form epv
  epv = np.empty(nLevels)
  for l in xrange(nLevels):
    epv[l] = np.dot(vort[l,:],theta_xyz[l,:])/rho[l]
  epv = epv*1e6 #to PVU
  return epv
  
def calc_thetaOnPVU_column(pvuVal, epv, theta, nLevels):
  #Ertel PV is (wLocal+wEarth).grad(tracer)/rho
  #pvuVal in [1.5,2.0] for dynamic tropopause
  (l,dl) = calcIndexOfValue(pvuVal,epv, nLevels) #index on 2PVU surface
  val = calcValueOfIndex(l,dl,theta)
  return val

#def fill_lstsqMatrix(A, xyzNbrs, xyz0, nHNbrs, nbrInds, level, nLevels):
def fill_lstsqMatrix(A, xyzNbrs, xyz0, nHNbrs, level, nLevels):
  #form the dxyz matrix for fitting a plane.
  
  for i in xrange(nHNbrs):
    for dir in xrange(3):
      #A[i,dir] = xyzNbrs[nbrInds[i],level,dir]-xyz0[level,dir]
      A[i,dir] = xyzNbrs[i,level,dir]-xyz0[level,dir]

  i=nHNbrs
  if (level>0): #have cell beneath
    for dir in xrange(3):
        A[i,dir] = xyz0[level-1,dir]-xyz0[level,dir]
    i=i+1
  if (level<nLevels-1): #have cell above
    for dir in xrange(3):
        A[i,dir] = xyz0[level+1,dir]-xyz0[level,dir]
    i=i+1
  return A

#def fill_lstsqVector(fillVec, level, vals0, valsNbrs, nHNbrs, nbrInds,nLevels):
def fill_lstsqVector(fillVec, level, vals0, valsNbrs, nHNbrs, nLevels):
  #form the dValue vector for fitting the lstsq plane from the column of information.
  #horizontal neighbors, bottom, then top
  row=0; #dirs = [0,1,2];
  for i in xrange(nHNbrs):
    #dval
    #fillVec[row] = valsNbrs[nbrInds[i],level]-vals0[level]
    fillVec[i] = valsNbrs[i,level]-vals0[level]
  row=nHNbrs
  if (level>0): #have cell beneath
    fillVec[row] = vals0[level-1]-vals0[level]
    row=row+1
  if (level<nLevels-1): #have cell above
    fillVec[row] = vals0[level+1]-vals0[level]
    row=row+1
  return fillVec

def calc_varsGradients_leastSquare(level, nVars, xyz0, xyzNbrs,vals0,valsNbrs,nHNbrs, nLevels):
  '''
  Pass in the columns of own and nbrs data.
  vals is the field values (eg theta) indexed as valsNbrs[variableInd, nbr, level]
  
  Calculating all of the derivs together should be more efficient since we do svd 1x
  
  Index return as [(dvar/dx, dvar/dy, dvar/dz), varIndex]
  '''
  
  #get layer neighbors. know vertical by level in nLevels (l=0 is bottom, l=nlevels-1 is top).
  #think about using vertically diagonal nbrs...
  nNbrs = nHNbrs+2 #for top and bottom
  if (level==0): #split if in case have only 1 level
    nNbrs = nNbrs-1;
  if (level==nLevels-1):
    nNbrs = nNbrs-1;
  
  #fit a plane: y-y0 = cx *(x-x0) + cy *(y-y0) +cz *(z-z0) to likely non-linear data
  #(dF/dxi) dxi = F-F0
  
  #Setup matrix system with horizontal then vertical 
  A = np.empty((nNbrs, 3), np.dtype('float64')) #should we specify dtype of float? what precision, np.dtype('float64')?
  y = np.empty((nNbrs,nVars),np.dtype('float64'))
  yTemp = np.empty(nNbrs, np.dtype('float64'))

  fill_lstsqMatrix(A, xyzNbrs, xyz0, nHNbrs, level, nLevels)
  for v in xrange(nVars):
    fill_lstsqVector(yTemp, level, vals0[v][:], valsNbrs[v][:,:], nHNbrs,nLevels)
    for i in xrange(nNbrs):
      y[i,v] = yTemp[i]

  #solve the least squares system for derivs
  #cx,cy,cz = np.linalg.lstsq(A, y)[0] #lstsq also gives residues, rank, singular values and can raise exception
  return np.linalg.lstsq(A, y)[0]

def calc_varsGradients_leastSquare_column(ht0, htNbrs, nVars, xyz0, xyzNbrs,vals0,valsNbrs,nHNbrs, nLevels):
  '''
  Pass in the columns of own and nbrs data.
  vals is the field values (eg theta) indexed as valsNbrs[variableInd, nbr, level].
  xyz and xyz0 get altered in what follows.
  
  Fit a plane: y-y0 = cx *(x-x0) + cy *(y-y0) +cz *(z-z0) to likely non-linear data
  ie dxi (dF/dxi) = F-F0 -> Ax=b
  solved in a least squares sense.
  
  There are two "optimizations here":
  -Calculating the derivs of ux,uy,uz,theta lets us do svd 1x per cell
  -We don't *need* to do an SVD for each cell center since (for l in [1,nLevels-1]):
    For level l, hl = altitude/RadiusEarth; xl = x0(1+hl) since sphere
    xla-xlb = x0a(1+hla)-x0b(1+hlb) = x0a-x0b+x0a*hla-x0b*hlb
    ~ (x0a-x0b)(1+hl?)...so, df/dxl=(1+hl?)*df/dx0
    So, we'll need to break the column into bottom, middle, and top SVDs valid for given mesh
    since the # of nbrs changes.
    For the middle section, we'll scale the F-F0 to the distance of the reference level
  
  Index return as [level,varIndex,(dvar/dx, dvar/dy, dvar/dz)]
  '''
  gradients = np.empty((nLevels,nVars,3), np.dtype('float64'))
  
  #for bottom and top, use single cell approach
  l=0;
  soln = calc_varsGradients_leastSquare(l, nVars, xyz0, xyzNbrs,vals0,valsNbrs,nHNbrs, nLevels) #form_curl(coeffs[:,0],coeffs[:,1],coeffs[:,2])
  for v in xrange(nVars):
    for i in xrange(3):
      gradients[l,v,i]=soln[i,v]
  l = nLevels-1;
  soln = calc_varsGradients_leastSquare(l, nVars, xyz0, xyzNbrs,vals0,valsNbrs,nHNbrs, nLevels)
  for v in xrange(nVars):
    for i in xrange(3):
      gradients[l,v,i]=soln[i,v]
      
  #for the cells in between, we can use the ht approx for a single SVD
  #get layer neighbors. know vertical by level in nLevels (l=0 is bottom, l=nlevels-1 is top).
  nNbrs = nHNbrs+2 #for top and bottom
  
  #Setup matrix system with horizontal then vertical
  nCellsL = nLevels-2 #cells in the "middle" of level
  A = np.empty((nNbrs, 3), np.dtype('float64')) #should we specify dtype of float? what precision, np.dtype('float64')?
  y = np.empty((nNbrs,nCellsL*nVars),np.dtype('float64')) #accessed as y[nbr, cell*nVars+v
  yTemp = np.empty(nNbrs, np.dtype('float64'))
  
  lRef = nLevels/2; #lRef=1
  fill_lstsqMatrix(A, xyzNbrs, xyz0, nHNbrs, lRef, nLevels)
  for l in xrange(1,nLevels-1):
    for v in xrange(nVars):
      fill_lstsqVector(yTemp, l, vals0[v][:], valsNbrs[v][:,:], nHNbrs,nLevels)
      for i in xrange(nNbrs):
        y[i,(l-1)*nVars+v] = yTemp[i]
  #
  #Rescale the dF's to the 'reference' locations
  for l in xrange(1,nLevels-1):
    h0 = ht0[lRef]; hl = ht0[l];
    scale = (R_earth+h0)/(R_earth+hl)
    for i in xrange(nNbrs):
      for v in xrange(nVars):
        y[i,(l-1)*nVars+v] *= scale
        
  ''' # 0
  for l in xrange(1,nLevels-1):
    #horizontal
    for i in xrange(nHNbrs):
      h0 = R_earth+.5*(htNbrs[i,lRef]+ht0[lRef]);
      hl = R_earth+.5*(htNbrs[i,l]+ht0[l]);
      scale = h0/hl
      for v in xrange(nVars):
        y[i,nVars*(l-1)+v] *= scale
    #vertical
    i = nHNbrs #under
    h0 = R_earth+.5*(ht0[lRef]+ht0[lRef-1]);
    hl = R_earth+.5*(ht0[l]+ht0[l-1]);
    scale = h0/hl
    for v in xrange(nVars):
      y[i,nVars*(l-1)+v] *= scale
    
    i = nHNbrs+1 #above
    h0 = R_earth+.5*(ht0[lRef]+ht0[lRef+1]);
    hl = R_earth+.5*(ht0[l]+ht0[l+1]);
    scale = h0/hl
    for v in xrange(nVars):
      y[i,nVars*(l-1)+v] *= scale
  '''
  
  #solve the least squares system for derivs
  #cx,cy,cz = np.linalg.lstsq(A, y)[0] #lstsq also gives residues, rank, singular values and can raise exception
  soln =  np.linalg.lstsq(A, y)[0]
  for l in xrange(1,nLevels-1):
    for v in xrange(nVars):
      for i in xrange(3):
        gradients[l,v,i]=soln[i,(l-1)*nVars+v]
  #gradients proper scale since dF rescaled
  
  return gradients

def driverErtel(state, hCell, hNbrs, nLevels):
  #put variable values in data structures for gradient calcs.
  #pass gradients to calc ertel pv.
  #Return the epv in the column at cell centers

  #hNbrs = state.cellsOnCells[hCell,:]
  nHNbrs = len(hNbrs)

  nVars = 4 #u's and theta
  vals0 = [None]*nVars
  vals0[0] = state.ux[hCell,:]
  vals0[1] = state.uy[hCell,:]
  vals0[2] = state.uz[hCell,:]
  vals0[3] = state.theta[hCell,:]
  xyz0 = state.xyz[hCell,:,:]

  xyzNbrs = state.xyz[hNbrs,:,:]
  valsNbrs = [None]*nVars
  valsNbrs[0] = state.ux[hNbrs,:]
  valsNbrs[1] = state.uy[hNbrs,:]
  valsNbrs[2] = state.uz[hNbrs,:]
  valsNbrs[3] = state.theta[hNbrs,:]
  
  epv = np.empty(nLevels,dtype=float)
  for l in xrange(nLevels):
    coeffs = calc_varsGradients_leastSquare(l, nVars, xyz0, xyzNbrs,vals0,valsNbrs,nHNbrs, nLevels)
    vort = form_curl(coeffs[:,0],coeffs[:,1],coeffs[:,2])
    #print coeffs; print vort
    epv[l] = calc_ertelPV(vort, coeffs[:,3], state.rho[hCell,l])
  return epv

def driverErtel_column(state, hCell, hNbrs, nLevels):
  #put variable values in data structures for gradient calcs.
  #pass gradients to calc ertel pv.
  #Return the epv in the column at cell centers

  #hNbrs = state.cellsOnCells[hCell,:]
  nHNbrs = len(hNbrs)

  nVars = 4 #u's and theta
  vals0 = [None]*nVars
  vals0[0] = state.ux[hCell,:]
  vals0[1] = state.uy[hCell,:]
  vals0[2] = state.uz[hCell,:]
  vals0[3] = state.theta[hCell,:]
  xyz0 = state.xyz[hCell,:,:] #.copy() #these get altered!
  ht0 = state.ht[hCell,:]

  valsNbrs = [None]*nVars
  valsNbrs[0] = state.ux[hNbrs,:]
  valsNbrs[1] = state.uy[hNbrs,:]
  valsNbrs[2] = state.uz[hNbrs,:]
  valsNbrs[3] = state.theta[hNbrs,:]
  xyzNbrs = state.xyz[hNbrs,:,:] #.copy()
  htNbrs = state.ht[hNbrs,:]
  
  epv = np.empty(nLevels,dtype=float)
  grads = calc_varsGradients_leastSquare_column(ht0, htNbrs, nVars, xyz0, xyzNbrs,vals0,valsNbrs,nHNbrs, nLevels)
  for l in xrange(nLevels):
    #gradients[l,v,i]
    vort = form_curl(grads[l,0,:], grads[l,1,:], grads[l,2,:])
    #print coeffs; print vort
    epv[l] = calc_ertelPV(vort, grads[l,3,:], state.rho[hCell,l])
  return epv

def get_allCellCenters(data,nCells,nLevels):
  #xyz coordinates of cell centers indexed as [hCell, level,xyz]
  xs = data.variables['xCell'][:] #surface layer
  ys = data.variables['yCell'][:]
  zs = data.variables['zCell'][:]
  zgrid = data.variables['zgrid'][:,:]

  xyz = np.empty((nCells,nLevels,3),dtype=float)
  ht = np.empty((nCells,nLevels), dtype=float)
  for i in xrange(nCells):
    for l in xrange(nLevels):
      ht[i,l] = .5*(zgrid[i,l]+zgrid[i,l+1])
      xyzs = [xs[i],ys[i],zs[i]];
      d = np.linalg.norm(xyzs);
      if (d<6e3):
        print "\nUhoh. Coordinate is probably on unit sphere. Scale by radius?\n"
        return
      fac = 1+ht[i,l]/d
      for dir in range(3):
        xyz[i,l,dir] = xyzs[dir]*fac
  #
  return (xyz,ht)

def get_domainCellCenters(data,hCells,nLevels):
  #xyz coordinates of cell centers indexed as [hCellIndex, level,xyz]
  xs = data.variables['xCell'][hCells] #surface layer
  ys = data.variables['yCell'][hCells]
  zs = data.variables['zCell'][hCells]
  zgrid = data.variables['zgrid'][hCells,:]
  
  nCells = len(hCells)
  xyz = np.empty((nCells,nLevels,3),dtype=float)
  ht = np.empty((nCells,nLevels), dtype=float)
  for i in xrange(nCells):
    xyzs = [xs[i],ys[i],zs[i]];
    d = np.linalg.norm(xyzs);
    if (d<6e3):
      print "\nUhoh. Coordinate is probably on unit sphere. Scale by radius?\n"
      return
    for l in xrange(nLevels):
      ht[i,l] = .5*(zgrid[i,l]+zgrid[i,l+1])
      fac = 1+ht[i,l]/d
      for dir in range(3):
        xyz[i,l,dir] = xyzs[dir]*fac
  #
  return (xyz,ht)

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
  print s
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

class dummyObject(object):
  #can add attributes by doing: a = dummyObject(); a.somefield = 1;
  pass

def loadFields(data, timeInd, hCells, nLevels): #, g2lMap):
  #Since we can't store whole fields in memory, get a subset to reduce repeated reads from file.
  #Return "state" that acts like loaded data field. Can add attributes to returned variable too
  
  import time #timing this function
  tStart = time.time()
  
  #Cell coordinates
  xyz, ht = get_domainCellCenters(data,hCells,nLevels)
  
  #velocity
  #The uReconstruct(X,Y,Z) are the horizontal velocity normal to the face reconstructed
  #to the cell center. For full velocity field, include w.
  ux = data.variables['uReconstructX'][timeInd,hCells,:]
  uy = data.variables['uReconstructY'][timeInd,hCells,:]
  uz = data.variables['uReconstructZ'][timeInd,hCells,:]

  n = data.variables['localVerticalUnitVectors'][hCells,:]
  wFace = data.variables['w'][timeInd,hCells,:]
  nCells = len(hCells)
  w = np.empty((nCells,nLevels))
  for l in xrange(nLevels):
    w[:,l] = .5*(wFace[:,l]+wFace[:,l+1])
  #ux[:,:] += n[:,0]*w[:,:] #multiplies element wise by row
  #uy[:,:] += n[:,1]*w[:,:]
  #uz[:,:] += n[:,2]*w[:,:]
  for l in xrange(nLevels):
    ux[:,l] += n[:,0]*w[:,l]
    uy[:,l] += n[:,1]*w[:,l]
    uz[:,l] += n[:,2]*w[:,l]

  #eqn of state vars
  #density
  rho = data.variables['rho'][timeInd,hCells,:]
  
  #temperature
  theta = data.variables['theta'][timeInd,hCells,:]
  
  #integration outputs have pressure. If from init, calc pressure
  #pSfc = data.variables['surface_pressure'][time,:]
  p = []
  try:
    p = data.variables['pressure_base'][timeInd,hCells,:]+data.variables['pressure_p'][timeInd,hCells,:]
    #break
  except:
    #init files don't have pressure field
    zz = data.variables['zz'][hCells,:]
    theta_m = data.variables['theta_m'][timeInd,hCells,:]
    p = calc_pressure_mpasEqnState(zz, theta_m, rho)
  
  #wrap up output so don't have extensive return tuple
  state = dummyObject()
  state.xyz=xyz; state.ht=ht;
  state.ux=ux; state.uy=uy; state.uz=uz;
  state.rho=rho; state.theta=theta;
  state.p=p; #state.pSfc = pSfc;

  #state.cellsOnCells = conn.get_domainHorizNbrs(data, hCells, g2lMap)
  
  tEnd = time.time()
  print "Loaded fields in {0} seconds\n".format(tEnd-tStart)
  
  return state

def driver_domains(nSeeds):
  #
  ncfname = '/arctic1/mduda/60km/x1.163842.output.2006-07-09_12.00.00.nc'
  #ncfname = '/home/nickszap/research/mpas/output.2010-10-23_00:00:00.nc'
  data = netCDF4.Dataset(ncfname,'r')
  
  nCellsTotal = len(data.dimensions['nCells'])
  nVerticesTotal = len(data.dimensions['nVertices']);
  nLevels = len(data.dimensions['nVertLevels'])
  nEdgesOnCell = data.variables['nEdgesOnCell'][:];
  cellsOnCell = data.variables['cellsOnCell'][:]-1;
  
  seed0 = 0
  #seed0 = np.argmax(data.variables['meshDensity'][:]) #seems like a decent heuristic
  cell2Site,seeds = conn.partition_max(seed0, cellsOnCell, nEdgesOnCell,nCellsTotal, nSeeds)
  
  for domainInd in xrange(nSeeds):
    #output domain mesh ------------------------
    cells = np.array(xrange(nCellsTotal))[cell2Site==seeds[domainInd]]
    nCells = len(cells)
    
    #open the output vtk file and write header.
    vtkfname = 'test'+str(domainInd)+'.vtk'
    fvtk = output_data.write_vtk_header_polydata(vtkfname, ncfname)
    
    #write nodes and cells
    output_data.write_vtk_polyHorizConn_domain(data, fvtk, cells, nEdgesOnCell,nVerticesTotal)
    
    #cell values and connectivity for this domain
    haloCells = conn.getHalo(seeds[domainInd], cell2Site, cellsOnCell, nEdgesOnCell, nCellsTotal)
    g2lCell = conn.make_global2localMap(cells, haloCells, nCellsTotal)
    c2c = conn.make_localDomainNbrs(nCells, cellsOnCell[cells,:], nEdgesOnCell[cells], g2lCell)
    neededCells = cells.tolist(); neededCells.extend(haloCells) #has to be domain then halo (as in g2l map)
    
    #I'm having memory errors.
    #gc.collect()
    
    #load data for domain and halo -----------------------------
    timeInd = 0
    print "Loading data for domain {0} with {1} cells\n".format(domainInd, len(neededCells))
    #print neededCells
    state = loadFields(data, timeInd, neededCells, nLevels)  
    
    #compute derived variables -------------------------------
    #theta on dynamic tropopause
    pv = np.empty((nCells,nLevels), dtype=float)
    #make_localDomainNbrs(nCells, cellsOnCell_local, nEdgesOnCell_local, g2lMap)
    for hCell in xrange(nCells):
      hNbrs = c2c[hCell,0:nEdgesOnCell[cells[hCell]]]
      pvColumn = driverErtel(state, hCell, hNbrs, nLevels)
      for l in range(nLevels):
        pv[hCell,l] = pvColumn[l]
    #
    pvuVal = 2.; pv = np.abs(pv) #questionable hack for southern hemisphere
    thetaVal = np.empty(nCells)
    for hCell in range(nCells):
      (l,dl) = output_data.calcIndexOfValue(pvuVal,pv[hCell,:], nLevels)
      thetaVal[hCell] = output_data.calcValueOfIndex(l,dl,state.theta[hCell,:])

    #write some cell data ----------------
    fvtk.write('\nCELL_DATA '+str(nCells)+'\n')
    output_data.write_levelData_float('theta_2pvu', fvtk, thetaVal, nCells)
    
    fvtk.close()
  data.close()
  
def driver_arctic():
  #plot epv on a polar cap
  
  ncfname = '/arctic1/nick/cases/163842/testDuda/x1.163842.output.2006-07-15_00.00.00.nc'
  #ncfname = '/arctic1/mduda/60km/x1.163842.output.2006-07-08_00.00.00.nc'
  #ncfname = '/home/nickszap/research/mpas/output.2010-10-23_00:00:00.nc'
  data = netCDF4.Dataset(ncfname,'r')

  nCellsTotal = len(data.dimensions['nCells'])
  nVerticesTotal = len(data.dimensions['nVertices']);
  nLevels = len(data.dimensions['nVertLevels'])
  nEdgesOnCell = data.variables['nEdgesOnCell'][:];
  cellsOnCell = data.variables['cellsOnCell'][:]-1;
  
  latThresh = 45.*np.pi/180.
  #latThresh = 70.*np.pi/180.
  latCell = data.variables['latCell'][:]
  cells = conn.gatherArcticCells(latCell, nCellsTotal, latThresh)
  nCells = len(cells)

  #open the output vtk file and write header.
  vtkfname = 'test.arctic.pv_approx.vtk'
  fvtk = output_data.write_vtk_header_polydata(vtkfname, ncfname)

  #write nodes and cells
  output_data.write_vtk_polyHorizConn_domain(data, fvtk, cells, nEdgesOnCell,nVerticesTotal)

  #cell values and connectivity for this domain
  haloCells = conn.get_arcticHalo(cells, latCell, latThresh, cellsOnCell, nEdgesOnCell)
  g2lCell = conn.make_global2localMap(cells, haloCells, nCellsTotal)
  c2c = conn.make_localDomainNbrs(nCells, cellsOnCell[cells,:], nEdgesOnCell[cells], g2lCell)
  neededCells = cells.tolist(); neededCells.extend(haloCells) #has to be domain then halo (as in g2l map)

  #I'm having memory errors.
  #gc.collect()

  #load data for domain and halo -----------------------------
  timeInd = 0
  print "Loading data for domain {0} with {1} cells\n".format('arctic', len(neededCells))
  #print neededCells
  state = loadFields(data, timeInd, neededCells, nLevels)  

  #compute derived variables -------------------------------
  #theta on dynamic tropopause
  pv = np.empty((nCells,nLevels), dtype=float)
  #make_localDomainNbrs(nCells, cellsOnCell_local, nEdgesOnCell_local, g2lMap)
  for hCell in xrange(nCells):
    hNbrs = c2c[hCell,0:nEdgesOnCell[cells[hCell]]]
    pvColumn = driverErtel(state, hCell, hNbrs, nLevels)
    #pvColumn = driverErtel_column(state, hCell, hNbrs, nLevels)
    for l in range(nLevels):
      pv[hCell,l] = pvColumn[l]
  #
  pvuVal = 2.; #pv = np.abs(pv) #don't need questionable hack for southern hemisphere
  thetaVal = np.empty(nCells)
  for hCell in xrange(nCells):
    (l,dl) = output_data.calcIndexOfValue(pvuVal,pv[hCell,:], nLevels)
    thetaVal[hCell] = output_data.calcValueOfIndex(l,dl,state.theta[hCell,:])

  #write some cell data ----------------
  fvtk.write('\nCELL_DATA '+str(nCells)+'\n')
  output_data.write_levelData_float('theta_2pvu', fvtk, thetaVal, nCells)

  fvtk.close()
  data.close()
  
def test_deriv():
  #we'll use analytic functions to test derivs
  
  ncfname = '/arctic1/mduda/60km/x1.163842.output.2006-07-08_00.00.00.nc'
  #ncfname = '/home/nickszap/research/mpas/output.2010-10-23_00:00:00.nc'
  data = netCDF4.Dataset(ncfname,'r')

  nCellsTotal = len(data.dimensions['nCells'])
  nVerticesTotal = len(data.dimensions['nVertices']);
  nLevels = len(data.dimensions['nVertLevels'])
  nEdgesOnCell = data.variables['nEdgesOnCell'][:];
  cellsOnCell = data.variables['cellsOnCell'][:]-1;
  
  c0 = 0
  cells = np.array([c0])
  nCells = len(cells)
  
  haloCells = cellsOnCell[c0,0:nEdgesOnCell[c0]]
  #cell values and connectivity for this domain
  g2lCell = conn.make_global2localMap(cells, haloCells, nCellsTotal)
  c2c = conn.make_localDomainNbrs(nCells, cellsOnCell[cells,:], nEdgesOnCell[cells], g2lCell)
  neededCells = cells.tolist(); neededCells.extend(haloCells) #has to be domain then halo (as in g2l map)

  #I'm having memory errors.
  #gc.collect()

  #load data for domain and halo -----------------------------
  timeInd = 0
  print "Loading data for domain {0} with {1} cells\n".format('arctic', len(neededCells))
  #print neededCells
  state = loadFields(data, timeInd, neededCells, nLevels)
  
  #edit state with analytic functions
  state.ux[:,:] = state.xyz[:,:,0]; print state.ux
  state.uy[:,:] = state.xyz[:,:,1]
  state.uz[:,:] = state.xyz[:,:,2]
  
  out = [None]*nCells
  for hCell in xrange(nCells):
    grads = []
    gInd = cells[hCell]
    hNbrs = c2c[hCell,0:nEdgesOnCell[gInd]]
    nHNbrs = len(hNbrs)
    
    nVars = 4 #u's and theta
    vals0 = [None]*nVars
    vals0[0] = state.ux[hCell,:]
    vals0[1] = state.uy[hCell,:]
    vals0[2] = state.uz[hCell,:]
    vals0[3] = state.theta[hCell,:]
    xyz0 = state.xyz[hCell,:,:] #.copy() #these get altered!
    #ht0 = state.ht[hCell,:]

    valsNbrs = [None]*nVars
    valsNbrs[0] = state.ux[hNbrs,:]
    valsNbrs[1] = state.uy[hNbrs,:]
    valsNbrs[2] = state.uz[hNbrs,:]
    valsNbrs[3] = state.theta[hNbrs,:]
    xyzNbrs = state.xyz[hNbrs,:,:] #.copy()
    #htNbrs = state.ht[hNbrs,:]
    
    for l in xrange(nLevels):
      grads.append(calc_varsGradients_leastSquare(l, nVars, xyz0, xyzNbrs,vals0,valsNbrs,nHNbrs, nLevels))
    #grads = calc_varsGradients_leastSquare_column(ht0, htNbrs, nVars, xyz0, xyzNbrs,vals0,valsNbrs,nHNbrs, nLevels)
    
    out[hCell] = grads
    print grads
    
  data.close()
  
  return out

