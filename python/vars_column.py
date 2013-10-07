#module for computing meteorological variables (especially from mpas fields)
#In "A Multi-Resolution Approach to Global Ocean Modeling", Tables A.5 and A.6 
#have variables at the center of the layer in the vertical

#imports
import numpy as np

import output_data

#print "Did you add w to uReconstruct(x,y,z) to make it full velocity?\n"

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
Po = 101325;
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
  fac = p/p0
  return theta*(fac**(Rd/Cp))

def calc_coriolis(lat):
  #coriolis parameter for lat in radians
  f = 2.*omeg_e*np.sin(lat)
  return f

def calcPressure_mpasEqnState(data,hcell,time,level):
  #calculate the pressure given the equation of state
  p0 = 1.e5;
  #num = Rd*data.variables['zz'][hcell,level]*data.variables['theta_m'][time,hcell,level]
  num = Rd*data.variables['zz'][hcell,level]*data.variables['theta_m'][time,hcell,level]*data.variables['rho'][time,hcell,level]
  fac = np.power(num/p0,gamma)
  return p0*fac

def calcPressure_column(data,hcell,time, nLevels):
  #return the pressure at cell centers. the IC files don't have pressure_p so do it other way
  p0 = 1.e5;
  #num = Rd*data.variables['zz'][hcell,0:nLevels]*data.variables['theta_m'][time,hcell,:] #zz is dimensioned levels+1 but top is junk
  num = Rd*data.variables['zz'][hcell,0:nLevels]*data.variables['theta_m'][time,hcell,:]*data.variables['rho'][time,hcell,:]
  fac = np.power(num/p0,gamma)
  return p0*fac
  
def calc_pressure(data, hcell, time, level):
  p = data.variables['pressure_base'][time,hcell,level]+data.variables['pressure_p'][time,hcell,level]
  return p

def calc_moistDensity(data, hcell, time, level):
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
  
def calc_slp(data, hcell, nLevels, time):
  #wrapper for calculating sea level (h=0m) pressure.
  #Given (unjustified) assumptions on lapse rate, hydrostatic,
  #the issue is how to find the pressure at below model levels. This boils down to
  #what reference level to move down from.
  
  #Following that described in slp_ben.m from Ryan Torn,
  #T is moved from 2km above ground to the ground using the standard lapse rate.
  #If surface cell is higher than 2km (eg Tibet), use that cell.
  #Advantage of height is to alleviate diurnal effects (on temperature at surface, right???)
  
  print "This might be buggy!!!"
  
  refHt = 2000.; #2km 
  
  #get t and p at sfc
  theta = data.variables['theta'][time,hcell,:]
  #p = data.variables['pressure_base'][time,hcell,:]+data.variables['pressure_p'][time,hcell,:]
  p = calcPressure_column(data,hcell,time,nLevels)
  t = np.empty(nLevels)
  for l in range(nLevels):
    t[l] = calc_temperature(theta[l], p[l])
  #
  hg = .5*(data.variables['zgrid'][hcell,0]+data.variables['zgrid'][hcell,1])
  
  #these values get changed below.
  #here we're just setting them to my heart's content so they're in the proper scope
  ph = 70000.; th = 280.; h = refHt;
  
  if (hg>refHt):
    #have to use cell as reference to move to sfc (could NaNs if slp not meaningful)
    h = hg
    #ph = p[0]
    th = t[0]
  else:
    #interpolate to height
    #zgrid is at faces of cells instead of midpts of layers.
    #need to adjust index into centered arrays since cell_l0=layer_l0+.5
    h = refHt
    
    (l,dl) = output_data.calcIndexOfValue(refHt,data.variables['zgrid'][hcell,:], nLevels+1)
    if (dl>=.5):
      #same level, just shift
      dl = dl-.5
    else:
      #have to shift into neighbor
      l = l-1
      dl = dl+.5
      
    #account for indices out of bounds. this shouldn't happen w/ height at center.
    #adjusting accounts for: vals[l+1]-vals[l] in calcValueOfIndex
    if (l>nLevels-2):
      #extrapolation up
      dl = dl+(l-(nLevels-1))
      l = nLevels-2;
    elif (l<0):
      #extrapolation down
      dl = dl+(l-0)
      l=0;
    
    #interpolate to height
    #ph = calcValueOfIndex(l,dl,p)
    th = output_data.calcValueOfIndex(l,dl,t)
  #
  
  #move temperature to sfc and calc slp
  tSurf = calc_temperatureFromLapse(th, hg-h)
  pSurf = p[0] #could also data.variables['surface_pressure'][time,ncells]
  tRef = calc_temperatureFromLapse(tSurf, -hg/2.) #half way between surface and sea level
  slp = calc_pressureHypsometric(tRef, pSurf, -hg)
  return slp

def calc_slp_surface(data,hcell,time):
  #using surface or first cell values
  
  #heights of faces of cell 0
  z01 = data.variables['zgrid'][hcell,0:2]
  
  #surface_pressure is at bottom face, right????
  pSurf = data.variables['surface_pressure'][time,hcell]
  if (abs(z01[0])<1):
    #we're at z=0
    return pSurf
  
  #we're either above or below sea level
  #get temperature at cell0
  p0 = calcPressure_mpasEqnState(data,hcell,time,0)
  theta0 = data.variables['theta'][time,hcell,0]
  t0 = calc_temperature(theta0, p0)
  #move temperature to surface
  z0 = .5*(z01[0]+z01[1])
  tSurf = calc_temperatureFromLapse(t0, z01[0]-z0)
  tRef = calc_temperatureFromLapse(tSurf, -.5*z01[0]) #half way between surface and sea level
  slp = calc_pressureHypsometric(tRef, pSurf, -z01[0])
  return slp

def calc_slp_arps(data, hcell, time):
  #pressure at 0 meters using the procedure in ARPS as shared to me by a one Lee Carlaw.
  #This is described in:
  #An alternative sea level pressure reduction and a statistical comparison of geostrophic wind estimates with observed #surface winds: benjamin and miller, 1990
  
  #Using hydrostatic and hypsometric equations, we need a T0
  #The paper way is to get T0 from 700mb. ARPS uses the ground level if Psurface<700
  #2km instead of pressure is used in slp_ben.m created April 2005 Ryan Torn, U. Washington
  pass
  
def calc_temperatureFromLapse(t0, dz):
  #return the temperature at the changed height in meters using standard atmospheric lapse rate.
  #dz is positive up and t decreases
  return t0-lapsesta*dz

def calc_slp_ben(tSurf, hSurf, pSurf):
  
  t0 = calc_temperatureFromLapse(tSurf,-.5*hSurf) #"reference temperature" halfway between ground and 0m
  p = calc_pressureHypsometric(t0, pSurf, -hSurf); #dz is positive up
  return p
  
def calc_slp_hydro(data, hcell, time):
  #pressure at 0 meters. Use hydrostatic with density fixed as at ground
  
  #dp/dz = -rho g
  ht = calc_htCellCenter(data, hcell, 0)
  rhoh = data.variables['rho'][time,hcell,0] #assume constant rho down? dry density???
  ph = calc_pressure(data, hcell, time, 0)
  dp_dz = -rhoh*g;
  p0 = ph-dp_dz*ht #where dz=ht-0m= ht
  return p0

def form_curl(u_xyz, v_xyz, w_xyz):
  #form the vector grad x u given all of the spatial derivs
  
  vort = np.empty(3)
  vort[0] = w_xyz[1]-v_xyz[2]
  vort[1] = u_xyz[2]-w_xyz[0]
  vort[2] = v_xyz[0]-u_xyz[1]
  return vort

def calc_vorticitySCVT():
  #using the normal and tangential velocities on the mesh
  pass

def calc_absoluteVorticity(data, time, hCell, level, nLevels):
  #return the absolute vorticity at a given cell center
  vort = calc_gradxu_leastSquares(data, time, hCell, level, nLevels)
  vort[2] = vort[2]+2.*omeg_e
  return vort
  
def calc_gradxu_leastSquares(data, time, hCell, level, nLevels):
  #Return local vorticity at a cell center in xyz space
  
  #we could do a number of diagnostics to see error (like comparing derivs to divergence)!!!
  print "Did you add w to uReconstruct(x,y,z) to make it full velocity?\n"
  u_xyz = np.empty((3,3),dtype=float) #du/dx,du/dy,du/dz \n dv/dx,...
  keys = ['uReconstructX','uReconstructY','uReconstructZ']
  i=0
  for key in keys:
    gradu = calc_varGradient_leastSquare(data, key, hCell, level, nLevels, time)
    u_xyz[i,:] = gradu[:]
    i = i+1
    
  vort = form_curl(u_xyz[0,:], u_xyz[1,:], u_xyz[2,:])
  return vort

def calc_ertelPV_column(data,time,hCell,nLevels):
  #Ertel PV is (wLocal+wEarth).grad(tracer)/rho
  #Do we care about the sign?????? Is it opposite in the southern hemisphere?
  
  #local and absolute vorticity
  vort = calc_lstsq_vort_column(data, time, hCell, nLevels) #vort[l,dirs]
  for l in xrange(nLevels):
    vort[l,2] += 2.*omeg_e #absolute with earth
  
  #gradient of tracer
  theta_xyz = calc_lstsqDeriv_column(data, time, hCell, 'theta', nLevels) #coeffs[l,dirs]
  
  #form epv
  rho = data.variables['rho'][time,hCell,:]
  epv = np.empty(nLevels)
  for l in xrange(nLevels):
    epv[l] = np.dot(vort[l,:],theta_xyz[l,:])/rho[l]
  epv = epv*1e6 #to PVU
  return epv
  
def calc_thetaOnEPV_column(data,time,hCell,nLevels):
  #Ertel PV is (wLocal+wEarth).grad(tracer)/rho
  
  #local and absolute vorticity
  vort = calc_lstsq_vort_column(data, time, hCell, nLevels) #vort[l,dirs]
  for l in xrange(nLevels):
    vort[l,2] += 2.*omeg_e #absolute with earth
  
  #gradient of tracer
  theta_xyz = calc_lstsqDeriv_column(data, time, hCell, 'theta', nLevels) #coeffs[l,dirs]
  
  #form epv
  rho = data.variables['rho'][time,hCell,:]
  epv = np.empty(nLevels)
  for l in xrange(nLevels):
    epv[l] = np.dot(vort[l,:],theta_xyz[l,:])/rho[l]  
  #looks like sign depends on hemisphere, but actually going by north/south is baffling about the equator
  epv = abs(epv*1e6)
  '''
  #going by hemisphere is tricky since I'm not sure what to assign and sign(0)=0
  sgn = np.sign(data.variables['zCell'][hCell])
  if (sgn==0.0):#python treats 0==0.0 as True as well
    sgn=1 #I'm not really sure what to do at equator!!!
  epv *= sgn*1e6
  '''
  (l,dl) = output_data.calcIndexOfValue(2.0,epv, nLevels) #index on 2PVU surface
  #print "For hcell {0}, dynamic tropopause at verticalIndex {1},{2}\n".format(hCell,l,dl)
  theta = data.variables['theta'][time,hCell,:]
  val = output_data.calcValueOfIndex(l,dl,theta)
  return val

'''
def calc_vertErtelPV_column(data, time, hCell, nLevels):
  #just "radial" component of 1/rho (vort)*grad(s)
  
  #MPAS pv is the vertical "potential-ed vorticity", ie divided by height
  #Ertel PV is the Dynamics PV with the passive tracer (eg potential temperature)
  #Return epv in PVUs (ie SI * 10^6)

  #MPAS computations are found in core_nhyd_atmos/mpas_atm_time_integration.F
  #Since the reconstruction of circulation is a big component of the grid (normal vels),
  #we'll use the MPAS computation of vorticity even though we don't get the horiz components!!!
  
  #vertical absolute vorticity by converting the MPAS "potential-ed vorticity"
  #following computations done in atm_compute_solve_diagnostics(dt, s, diag, grid) of
  #mpas_atm_time_integration.F
  vpv = np.empty(nLevels);
  mpv = data.variables['pv_cell'][time, hCell,:]
  h = data.variables['rho_zz'][time,hCell,:] #why do they use this????
  vpv = mpv*h
  
  #radial grid for finite differences
  tmp = data.variables['zgrid'][hCell,:]
  z = np.empty(nLevels)
  for l in range(nLevels):
    z[l] = .5*(tmp[l+1]+tmp[l]);
  
  #(radial) gradient of tracer
  s = data.variables['theta'][time,hCell,:]
  grads = np.empty(nLevels)
  #centered difference except for top and bottom
  l = 0; grads[l] = (s[l+1]-s[l])/(z[l+1]-z[l]);
  l = nLevels-1; grads[l] = (s[l-1]-s[l])/(z[l-1]-z[l])
  for l in range(nLevels):
    grads[l] = (s[l+1]-s[l-1])/(z[l+1]-z[l-1])
  
  rho = data.variables['rho'][time,hCell,:]
  vpv = vpv*grads/rho #elementwise if vectors, right?
    
  return vpv;
'''

def calc_varGradient_leastSquare_data(data, key, hCell, level, nLevels, time):
  '''
  Calculate the gradient of a variable defined by 'key'.
  Remember that the mesh conn may be 1-indexed (as it is in the file)
  Return (dvar/dx, dvar/dy, dvar/dz)
  Should we normalize xyz by, say, radius??? Work in km?
  '''
  
  #get layer neighbors. know vertical by level in nLevels (l=0 is bottom, l=nlevels-1 is top).
  #think about using vertically diagonal nbrs...
  nbrs = getHorizNbrs(data, hCell);
  nNbrs = len(nbrs)+2 #for top and bottom
  if (level==0): #split if in case have only 1 level
    nNbrs = nNbrs-1;
  if (level==nLevels-1):
    nNbrs = nNbrs-1;
  
  #fit a plane: y-y0 = cx *(x-x0) + cy *(y-y0) +cz *(z-z0) to likely non-linear data
  #(dF/dxi) dxi = F-F0
  
  #Setup matrix system
  A = np.empty((nNbrs, 3), np.dtype('float64')) #should we specify dtype of float? what precision, np.dtype('float64')?
  y = np.empty(nNbrs,np.dtype('float64'))
  
  xyz0 = calc_cellCenterCoord(data,hCell,level)
  var0 = data.variables[key][time,hCell, level]
  #horizontal
  i=0
  for cell in nbrs:
    var = data.variables[key][time,cell, level]
    y[i] = var-var0
    xyz = calc_cellCenterCoord(data,cell,level)
    #create dx,dy,dz
    for dir in xrange(3):
      A[i,dir] = xyz[dir]-xyz0[dir]
    #
    i = i+1;
  '''
  var = data.variables[key][time,nbrs,level]
  y[:] = var-var0
  '''
  #vertical
  if (level>0): #have someone beneath
    var = data.variables[key][time,hCell, level-1]
    y[i] = var-var0
    xyz = calc_cellCenterCoord(data,hCell,level-1)
    #create dx,dy,dz
    for dir in xrange(3):
      A[i,dir] = xyz[dir]-xyz0[dir]
    #
    i = i+1;
  if (level<nLevels-1): #have someone above
    var = data.variables[key][time,hCell, level+1]
    y[i] = var-var0
    xyz = calc_cellCenterCoord(data,hCell,level+1)
    #create dx,dy,dz
    for dir in xrange(3):
      A[i,dir] = xyz[dir]-xyz0[dir]
    #
    i = i+1;
  
  #solve the least squares system for derivs
  #cx,cy,cz = np.linalg.lstsq(A, y)[0] #lstsq also gives residues, rank, singular values and can raise exception
  return np.linalg.lstsq(A, y)[0] 

def lstsq_qr():
    print "Not written!\n"
    #solving Ax=b with QR-Factorization
    Q,R = linalg.qr(A) # qr decomposition of A
    Qb = dot(Q.T,b) # computing Q^T*b (project b onto the range of A)
    x_qr = linalg.solve(R,Qb) # solving R*x = Q^T*b

def fill_lstsqMatrix():
  #form the dxyz matrix for fitting a plane
  pass

def fill_lstsqVector(fillVec, level, vals0, valsNbrs, nHNbrs, nLevels):
  #form the dValue vector for fitting the lstsq plane from the column of information
  #horizontal neighbors, bottom, then top
  row=0; #dirs = [0,1,2];
  for i in xrange(nHNbrs):
    #dval
    fillVec[row] = valsNbrs[i,level]-vals0[level]
  print "not written yet"
  
def calc_varGradient_leastSquare(level,xyz0, xyzNbrs,vals0,valsNbrs,nHNbrs, nLevels):
  '''
  Pass in the columns of own and nbrs data.
  Return (dvar/dx, dvar/dy, dvar/dz)
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
  y = np.empty(nNbrs,np.dtype('float64'))
  
  #horizontal
  row=0; #dirs = [0,1,2];
  for i in xrange(nHNbrs):
    #dval
    y[row] = valsNbrs[i,level]-vals0[level]
    #dx
    for dirs in xrange(3):
      A[row,dirs] = xyzNbrs[i,level,dirs]-xyz0[0,level,dirs];
    
    row = row+1;
    
  #vertical
  if (level>0): #have someone beneath
    y[row] = vals0[level-1]-vals0[level]
    for dirs in xrange(3):
      A[row,dirs] = xyz0[0,level-1,dirs]-xyz0[0,level,dirs];
    
    row = row+1
  if (level<nLevels-1): #have someone above
    y[row] = vals0[level+1]-vals0[level]
    for dirs in xrange(3):
      A[row,dirs] = xyz0[0,level+1,dirs]-xyz0[0,level,dirs];
    
    row = row+1
    
  #solve the least squares system for derivs
  #cx,cy,cz = np.linalg.lstsq(A, y)[0] #lstsq also gives residues, rank, singular values and can raise exception
  return np.linalg.lstsq(A, y)[0]

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
  
  #horizontal
  row=0; #dirs = [0,1,2];
  for i in xrange(nHNbrs):
    #dval
    for v in xrange(nVars):
      y[row,v] = valsNbrs[v,i,level]-vals0[v,level]
    #dx
    for dirs in xrange(3):
      A[row,dirs] = xyzNbrs[i,level,dirs]-xyz0[0,level,dirs];
    
    row = row+1;
    
  #vertical
  if (level>0): #have someone beneath
    for v in xrange(nVars):
      y[row,v] = vals0[v,level-1]-vals0[v,level]
    for dirs in xrange(3):
      A[row,dirs] = xyz0[0,level-1,dirs]-xyz0[0,level,dirs];
    
    row = row+1
  if (level<nLevels-1): #have someone above
    for v in xrange(nVars):
      y[row,v] = vals0[v,level+1]-vals0[v,level]
    for dirs in xrange(3):
      A[row,dirs] = xyz0[0,level+1,dirs]-xyz0[0,level,dirs];
    
    row = row+1
    
  #solve the least squares system for derivs
  #cx,cy,cz = np.linalg.lstsq(A, y)[0] #lstsq also gives residues, rank, singular values and can raise exception
  return np.linalg.lstsq(A, y)[0]

def calc_lstsqDeriv_column(data, time, hCell, key, nLevels):
  #fit a plane: y-y0 = cx *(x-x0) + cy *(y-y0) +cz *(z-z0) to likely non-linear data
  #(dF/dxi) dxi = F-F0
  #To do the least squares, we need: coords and vals for self and nbrs.
  #we go by column to reduce the number of file reads (time) and assume the horizontal mesh
  #might get big enough that we can't store it all in memory.
  #Index return as [level,cx or cy or cz]
  
  #neighbors
  hNbrs = getHorizNbrs(data, hCell).tolist();
  nHNbrs = len(hNbrs)
  
  #coordinates will eventually be relative to center cell
  xyzNbrs = calc_cellCenterColumn(data,hNbrs,nLevels)
  xyz0 = calc_cellCenterColumn(data,[hCell],nLevels)
  
  #values will also eventually relative to center cell 
  valsNbrs = data.variables[key][time,hNbrs,:]
  vals0 = data.variables[key][time,hCell,:]
  
  #lstsq solution
  #consider doing linalg.svd and then computing all of the lstsq gradients for that cell!!!
  coeffs = np.empty((nLevels,3)); #dirs = [0,1,2];
  for l in xrange(nLevels):
    cxyz = calc_varGradient_leastSquare(l, xyz0,xyzNbrs,vals0,valsNbrs, nHNbrs, nLevels);
    for dirs in xrange(3):
      coeffs[l,dirs] = cxyz[dirs]
  
  return coeffs

def calc_lstsqDeriv_vars_column(data, time, hCell, keys, nKeys, nLevels):
  #fit a plane: y-y0 = cx *(x-x0) + cy *(y-y0) +cz *(z-z0) to likely non-linear data
  #(dF/dxi) dxi = F-F0
  #To do the least squares, we need: coords and vals for self and nbrs.
  #we go by column to reduce the number of file reads (time) and assume the horizontal mesh
  #might get big enough that we can't store it all in memory.
  #Index return as [level,cx or cy or cz]
  
  #neighbors
  hNbrs = getHorizNbrs(data, hCell).tolist();
  nHNbrs = len(hNbrs)
  
  #coordinates will eventually be relative to center cell
  xyzNbrs = calc_cellCenterColumn(data,hNbrs,nLevels)
  xyz0 = calc_cellCenterColumn(data,[hCell],nLevels)
  
  #values will also eventually relative to center cell 
  valsNbrs = np.empty((nKeys,nHNbrs,nLevels));
  vals0 = np.empty((nKeys,nLevels))
  for k in xrange(nKeys):
    valsNbrs[k,:,:] = data.variables[keys[k]][time,hNbrs,:]
    vals0[k,:] = data.variables[keys[k]][time,hCell,:]
  
  print("continue here!!!!!!!\n")
  
  #lstsq solution
  #consider doing linalg.svd and then computing all of the lstsq gradients for that cell!!!
  coeffs = np.empty((nLevels,3)); #dirs = [0,1,2];
  for l in xrange(nLevels):
    cxyz = calc_varGradient_leastSquare(l, xyz0,xyzNbrs,vals0,valsNbrs, nHNbrs, nLevels);
    for dirs in xrange(3):
      coeffs[l,dirs] = cxyz[dirs]
  
  return coeffs

def calc_lstsq_vort_column(data, time, hCell, nLevels):
  #curl of velocity field
  #index returned as [level,vortx vorty or vortz]
  print "Did you add w to uReconstruct(x,y,z) to make it full velocity?\n"
  u_xyz = calc_lstsqDeriv_column(data, time, hCell, 'uReconstructX', nLevels);
  v_xyz = calc_lstsqDeriv_column(data, time, hCell, 'uReconstructY', nLevels);
  w_xyz = calc_lstsqDeriv_column(data, time, hCell, 'uReconstructZ', nLevels);
  
  vort = np.empty((nLevels,3)); #dirs = [0,1,2]
  for l in xrange(nLevels):
    curl = form_curl(u_xyz[l,:],v_xyz[l,:],w_xyz[l,:])
    for dirs in xrange(3):
      vort[l,dirs] = curl[dirs]
  return vort
  
def getHorizNbrs(data, hCell):
  cellsOnCell = data.variables['cellsOnCell'][hCell] #this has "junk" to maxEdges dimension
  nNbrs = data.variables['nEdgesOnCell'][hCell]
  return cellsOnCell[0:nNbrs]-1 #0-indexing

def calc_htCellCenter(data, hcell, l):
  #radial height of cell center from 0 layer
  ht = .5*(data.variables['zgrid'][hcell,l]+data.variables['zgrid'][hcell,l+1])
  return ht

def calc_htCellCenter_column(data,hCell,nLevels):
  #radial height from sea level
  zInterface = data.variables['zgrid'][hCell,:]
  z = np.empty(nLevels)
  for l in xrange(nLevels):
    z[l] = .5*(zInterface[l]+zInterface[l+1])
  return z
    
def calc_cellCenterCoord(data, hcell, l):
  #return the xyz coordinate of the midpt of the cell.
  #Although the user's guide says that the coords are given on the unit sphere.
  #it seems that they're actually on the spherical earth from the output files.
  
  xyzS = np.array([data.variables['xCell'][hcell], data.variables['yCell'][hcell], data.variables['zCell'][hcell]]);
  d = np.linalg.norm(xyzS) #2-norm is sqrt sum of squares
  if (d<6e3):
    print "\nUhoh. Coordinate is probably on unit sphere. Do we just need to scale by radius?\n"
    return
    
  #project radially out, ie length*dir
  #ht = .5*(data.variables['zgrid'][hcell,l+1]+data.variables['zgrid'][hcell,l])
  ht = calc_htCellCenter(data,hcell,l) #will ht of cell center ever be different than avg of faces?
  xyz = xyzS*(1.+ht/d) #xyzS/d is unit vector.
  return xyz
  
def calc_cellCenterColumn(data,hCells,nLevels):
  #return the xyz coordinates of the nLevels cells in the hCells column
  #Return xyzCellCenters where index returned mat as xyzc[nbr#: 0 1 ..., vertical level, x y or z]
  
  xyzS = np.matrix([data.variables['xCell'][hCells], data.variables['yCell'][hCells], data.variables['zCell'][hCells]]); #surface face xyz
  d = np.linalg.norm(xyzS[:,0]) #2-norm is sqrt sum of squares
  if (d<6e3):
    print "\nUhoh. Coordinate is probably on unit sphere. Do we just need to scale by radius?\n"
    return
  
  htFaces = data.variables['zgrid'][hCells,:]
  nCells = len(hCells); xyzc = np.empty((nCells,nLevels,3),dtype=float) #xyz of cell centers
  for i in xrange(nCells):
    for l in xrange(nLevels):
      htc = .5*(htFaces[i,l]+htFaces[i,l+1])
      cellxyz = xyzS[:,i]*(1.+htc/d); #scale radially out
      #inds = [0,1,2]; xyzc[i,l,inds] = cellxyz[inds]; #something about matrix isn't working here
      for dirs in xrange(3):
        xyzc[i,l,dirs] = cellxyz[dirs,0]
  #
  return xyzc
      
def leastSquaresExample():
  #From http://docs.scipy.org/doc/numpy/reference/generated/numpy-linalg-lstsq-1.py
  
  # Fit a line, ``y = mx + c``, through some noisy data-points:

  x = np.array([0, 1, 2, 3])
  y = np.array([-1, 0.2, 0.9, 2.1])

  # By examining the coefficients, we see that the line should have a
  # gradient of roughly 1 and cut the y-axis at, more or less, -1.

  # We can rewrite the line equation as ``y = Ap``, where ``A = [[x 1]]``
  # and ``p = [[m], [c]]``.  Now use `lstsq` to solve for `p`:

  A = np.vstack([x, np.ones(len(x))]).T
  A
  # array([[ 0.,  1.],
  # [ 1.,  1.],
  # [ 2.,  1.],
  # [ 3.,  1.]])

  m, c = np.linalg.lstsq(A, y)[0]
  print m, c
  # 1.0 -0.95

  # Plot the data along with the fitted line:

  import matplotlib.pyplot as plt
  plt.plot(x, y, 'o', label='Original data', markersize=10)
  plt.plot(x, m*x + c, 'r', label='Fitted line')
  plt.legend()
  plt.show()
  