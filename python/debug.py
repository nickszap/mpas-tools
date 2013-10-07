#python modules
import numpy as np

#my modules
import output_data
import vars

class dummyObject(object):
  #can add attributes by doing: a = dummyObject(); a.somefield = 1;
  pass

class fasterObject(object):
  __slots__ = ('variables')

def evalFct_threaded(ana, data, key, time, hCell, l):
  #playing around with python's multithreading
  
  import multiprocessing as mproc
  
  return ana

def testLSTSQ(data,nCellsRMS,nLevelsRMS):
  #use analytic functions to check accuracy of our gradients.
  #least squares function accesses the data.variables[key] structure so that we don't
  #need to load all nCells x nLevels (x time) values.
 
  #This means that we would have to write over some time varying cell field with our analytic fcts.
  #However, we open the files as write to read only.
  #So, we'll rig up a data structure that includes a dictionary instead.
  
  #calculate analytic fct over domain
  ana = fasterObject()
  key = 'key'; time = 0;
  nCells = len(data.dimensions['nCells']); nLevels = len(data.dimensions['nVertLevels'])
  ana.variables = {key:np.empty((1,nCells,nLevels)), 
                  'xCell':data.variables['xCell'][:],'yCell':data.variables['yCell'][:],'zCell':data.variables['zCell'][:],
                  'zgrid':data.variables['zgrid'][:], 
                  'cellsOnCell':data.variables['cellsOnCell'][:], 'nEdgesOnCell':data.variables['nEdgesOnCell'][:]}
  
  for hCell in range(nCells):
    print str(hCell)+'\n'
    top = min(nLevels,nLevelsRMS+2)
    for l in range(top):
      xyz = vars.calc_cellCenterCoord(data, hCell, l)
      (val,ddxyz) = testLSTQ_f(xyz)
      ana.variables[key][time,hCell,l] = val
  
  #compare least squares approx of deriv to exact
  rmsxyz = np.zeros(3)
  for hCell in range(nCellsRMS):
    for l in range(nLevelsRMS):
      xyz = vars.calc_cellCenterCoord(data, hCell, l)
      (val,ddxyz_e) = testLSTQ_f(xyz) #exact
      ddxyz_a = vars.calc_varGradient_leastSquare(ana, key, hCell, l, nLevels, time) #approximate
      
      #sum squares
      ddxyz = ddxyz_a-ddxyz_e
      
      rmsxyz = rmsxyz+ddxyz*ddxyz
      print "ddxyz: {0} for hcell,level: ({1},{2}) \n".format(ddxyz,hCell,l)
  #
  #mean
  rmsxyz = rmsxyz/(nCellsRMS*nLevelsRMS);
  #root
  rmsxyz = np.sqrt(rmsxyz)
  
  s = 'RMS of lstsq deriv over 3d field in xyz: {0}\n'.format(rmsxyz)
  print s
  
def testLSTQ_f(xyz):
  #analytic functions to test least squares.
  #return (val, (dval/dx, dval/dy, dval/dz))
  
  fVersion = 4
  
  #plane gets worst accuracy of ~ 1e-11
  if (fVersion==1):
    #ci*xi+offset
    c = np.array([3.,-1.,20.])
    val = np.dot(xyz,c) + 22.4;
    return (val, c)
  #
  
  #squares gets errors over 50. we're not really "about 
  if (fVersion==2):
    #ci*xi^2+offset
    c = np.array([3.,-1.,20.])
    fac = 1.e3;
    val = np.dot(xyz*xyz/(fac*fac),c) + 22.4;
    return (val, 2*c*xyz/(fac*fac))
  #
  
  #xyz gets errors of 1.6e-5. 
  if (fVersion==3):
    fac = 1.e6
    xyz = xyz/fac;
    val = xyz[0]*xyz[1]*xyz[2];
    fac = fac*fac*fac; deriv = (xyz[1]*xyz[2]/fac,xyz[0]*xyz[2]/fac,xyz[1]*xyz[0]/fac)
    return (val,deriv)
  #
  
  #sinusoids: rms errors of (dval/dxyz): 0.00063592  0.00068587  0.000641 (refLen=300km)
  #8.51791496e-05   8.49776813e-05   6.93804044e-05 (refLen=1200km)
  if (fVersion==4):
    refLen = 60.e3*20.; scale = 2.*np.pi/refLen
    val = np.sin(scale*xyz[0])+np.cos(scale*xyz[1])+3.*np.sin(scale*xyz[2])
    deriv = (scale*np.cos(scale*xyz[0]),-scale*np.cos(scale*xyz[1]), scale*3.*np.cos(scale*xyz[2]))
    return (val,deriv)
  

if __name__=='__main__':
  #test out the least squares
  fpath = '/arctic1/nick/cases/40962/output.40962.2006-06-01_00.00.00.nc'
  data = output_data.open_netcdf_data(fpath)
  nLevels =  len(data.dimensions['nVertLevels'])
  nCells = len(data.dimensions['nCells'])
  
  nCellsRMS = 20; nLevelsRMS = 8
  testLSTSQ(data,nCellsRMS,nLevelsRMS)
  