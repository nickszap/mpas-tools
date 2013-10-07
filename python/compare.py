
import netCDF4
import numpy as np
import datetime as dt
from time import strftime
'''
import matplotlib #comment out matplotlib above!
matplotlib.use('Agg')
import matplotlib.pyplot as plt2
'''
import matplotlib.pyplot as plt

import os

import cfsr
import output_data

def loc_initFile(t0):
  #file names are like: 'vert_sfc.2006-07-10_00.nc'
  #return the full path to the file of the requested datetime time
  '''
  fdir = '/arctic1/nick/cases/cfsr/output/'
  #timeString = cfsr.form_cfsrTimeString(t0)
  timeString = formMPASTime(t0)
  #fname = 'x4.vert_sfc.'+timeString+'.nc'
  #fname = 'vert_sfc.'+timeString+'.nc'
  fname = 'x4.cfsr.output.'+timeString+'.nc'
  '''
  '''
  fdir = '/arctic1/nick/cases/cfsr/output/'
  timeString = formMPASTime(t0)
  fname = 'x4.cfsr.output.'+timeString+'.nc'
  '''
  fdir = '/arctic1/nick/cases/cfsr/2007/'
  timeString = formMPASTime(t0)
  fname = 'x1.cfsr.output.'+timeString+'.nc'

  return fdir+fname

def formMPASTime(t0):
  tTuple = dt.datetime.timetuple(t0);
  s = strftime('%Y-%m-%d_%H.%M.%S', tTuple)
  return s

def mpasOutputName(t0, iScheme):
  #file names are like output.163842.2006-07-15_00.00.00.nc
  '''
  fdir = '/arctic1/nick/cases/163842/testDuda/'
  tTuple = dt.datetime.timetuple(t0);
  s = strftime('%Y-%m-%d_%H.%M.%S', tTuple)
  fname = 'x1.163842.output.'+s+'.nc'
  '''

  s = formMPASTime(t0)

  fdir = '/'
  fname = "unknown option"

  '''
  if (iScheme==0):
    fdir = '/arctic1/nick/cases/v1.0/x4/2week/kf/'
    fname = 'x4.kf.output.'+s+'.nc'
  elif(iScheme ==1):
    fdir = '/arctic1/nick/cases/v1.0/x4/2week/tiedtke/'
    fname = 'x4.tiedtke.output.'+s+'.nc'
  '''
  '''
  if (iScheme==0):
    fdir = '/arctic1/nick/cases/v1.0/'
    fname = 'x1.163842.gwdo.output.'+s+'.nc'
  elif(iScheme ==1):
    fdir = '/arctic1/nick/cases/v1.0/'
    fname = 'x1.tiedtke.output.'+s+'.nc'
  '''
  '''
  if (iScheme==0):
    fdir = '/arctic1/nick/cases/v1.0/x4/august/kf/v1.1/'
    fname = 'x4.kf.output.'+s+'.nc'
  elif(iScheme ==1):
    fdir = '/arctic1/nick/cases/v1.0/x4/august/tiedtke/'
    fname = 'x4.t.output.'+s+'.nc'
  '''
  if (iScheme==0):
    fdir = '/arctic1/nick/cases/v1.0/2007/'
    fname = 'x1.tiedtke.output.'+s+'.nc'
  
  return fdir+fname

def differenceFileName(t0):
  #file names are like: 'diff.2006-07-10_00.nc'
  s = cfsr.form_cfsrTimeString(t0)
  fname = 'diff.'+s+'.nc'
  return fname

def example_extractFileFields(ncIn, ncOut):
  #extracting just the fields of interest to a new file will reduce the 
  #network transfer time
  cmd = 'ncks -v theta '+str(ncIn)+' '+str(ncOut)
  os.system(cmd)

def calc_errorOutput_refInit_variable(key, mpasData, t0, deltat, timeInds): #, lev):
  #return a list for each statistic of length ntimes

  s = "\nStatistics for MPAS-Init for variable {0}\n".format(key)
  s += "timeInd mean stdev max(diff) min(diff)\n"
  print s; s = ''

  m = []; dev = []; maxd = []; mind = []

  for i in timeInds:
  #for i in xrange(4,5,1):
    t = t0+i*deltat
    initFile = loc_initFile(t)
    #print "Reference file is {0}\n".format(initFile)
    initData = netCDF4.Dataset(initFile,'r')

    varInit = initData.variables[key][0,:] #[0,:,lev] #[0,:]
    varMPAS = mpasData.variables[key][i,:] #[i,:,lev] #[i,:]
    varMPAS -= varInit
    #varMPAS.flatten()
    m.append(np.mean(varMPAS)); dev.append(np.std(varMPAS))
    maxd.append(np.max(varMPAS)); mind.append(np.min(varMPAS))
    #s += "t{0} {1}  {2}  {3}  {4}\n".format(i,m,dev,maxd, mind)

    initData.close()

  return (m,dev,maxd,mind)

def example_stat_files():

  outFileTimes = [dt.datetime(2007,8,01,00), dt.datetime(2007, 8, 8, 00)]

  #schemes = ['x4_kf', 'x4_t'];
  #schemes = ['x1_kf', 'x1_t'];
  schemes = ['x1_t']#, 'x1_t'];
  nSchemes = len(schemes)

  #keys_pressure = ['temperature']
  #pressureLevs = ['500hPa']
  keys_pressure = ['temperature', 'height', 'uzonal', 'umeridional', 'vorticity']
  pressureLevs = ['200hPa', '500hPa', '850hPa']
  #keys_levels = ['theta']
  #modelLevs = [5]
  #keys_levels = ['uReconstructZonal', 'uReconstructMeridional', 'rho', 'rh', 'theta']
  #modelLevs = [3,10,20,30]
  for kp in keys_pressure:
  #for key in keys_levels:
    for pl in pressureLevs:
    #for ml in modelLevs:
      key = kp+'_'+pl
      mv = [[] for i in xrange(nSchemes)]; devv = [[] for i in xrange(nSchemes)];
      maxdv = [[] for i in xrange(nSchemes)]; mindv = [[] for i in xrange(nSchemes)]
      for iScheme  in xrange(nSchemes):
        for fTime in outFileTimes:
          t0 = fTime
          mpasFile = mpasOutputName(t0, iScheme)
          #mpasFile = '/arctic1/nick/cases/v1.0/x4/2week/tiedtke/x4.tiedtke.output.2006-07-24_12.00.00.nc'
          print "Output file is {0}".format(mpasFile)
          mpasData = netCDF4.Dataset(mpasFile,'r')

          #fields: config_frames_per_outfile tells us when we need to get new file
          #config_output_interval time separation between values
          #nTimeInds = 28
          nTimeInds = len(mpasData.dimensions['Time'])
          #h6 = dt.timedelta(hours=6) #can also get deltaTOutput from netcdf attributes
          deltat = dt.timedelta(hours=6) #timestep of model output
          #timeInds = range(2,nTimeInds,4) #index of mpas file to time of existing reference file
          timeInds = range(0,nTimeInds,4)
          
          initFile1 = loc_initFile(t0+timeInds[0]*deltat)
          initFile2 = loc_initFile(t0+timeInds[1]*deltat)
          print "Reference file for first 2 times would be:\n{0}\n{1}".format(initFile1, initFile2)

          (m,dev,maxd,mind) = calc_errorOutput_refInit_variable(key, mpasData, t0, deltat, timeInds) #, ml)
          mv[iScheme].extend(m); devv[iScheme].extend(dev);
          maxdv[iScheme].extend(maxd); mindv[iScheme].extend(mind)

          mpasData.close()
          
      #print mv; print devv; print maxdv; print mindv;
      print(mv[0][:]);
      if (nSchemes>1):
        print(mv[1][:]);

      tfname = formMPASTime(outFileTimes[0])[0:13]
      titlePlot = 'x1_mean'+key; fname = titlePlot+'_'+tfname
      plot_save_lists(mv, nSchemes, schemes, ' ', fname)
      titlePlot = 'x1_stdev'+key; fname = titlePlot+'_'+tfname
      plot_save_lists(devv, nSchemes, schemes, ' ', fname)
      titlePlot = 'x1_maxd'+key; fname = titlePlot+'_'+tfname
      plot_save_lists(maxdv, nSchemes, schemes, ' ', fname)
      titlePlot = 'x1_mind'+key; fname = titlePlot+'_'+tfname
      plot_save_lists(mindv, nSchemes, schemes, ' ', fname)

def plot_save_lists(vals, nVars, labels, title, fname):
  #import matplotlib #comment out matplotlib above!
  #matplotlib.use('Agg')
  #import matplotlib.pyplot as plt2
  
  fig1 = plt2.figure()
  plt2.hold(True)
  for v in xrange(nVars):
    plt2.plot(vals[v][:], label=labels[v])
  
  #save instead of displaying
  plt2.title(title)
  #plt2.savefig(figName+'.png')
  plt2.autoscale(enable=True, axis='both', tight=True)
  plt2.legend()
  plt2.savefig(fname+'.png')
  plt2.close()

def example_difference():
  #For this, we have mpas outputs and CFSR ICs every 6 hours.
  #The files are on the same horizontal and vertical mesh
  '''
  if (True):
    reload(matplotlib)
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt2
  else:
    plt2 = plt
  '''
  import matplotlib #comment out matplotlib above!
  matplotlib.use('Agg')
  import matplotlib.pyplot as plt2
  
  t0 = dt.datetime(2006,7,15,0)
  mpasFile = mpasOutputName(t0)
  mpasData = netCDF4.Dataset(mpasFile,'r')
  
  #fields: config_frames_per_outfile tells us when we need to get new file
  #config_output_interval time separation between values
  nTimeInds = 28
  h6 = dt.timedelta(hours=6)
  
  #store difference fields in a copy of output file
  diffFile = differenceFileName(t0)
  if (False):
    #we'll store the difference fields in a copy of the output file for time animation
    cmd = 'cp '+mpasFile+' '+diffFile
    os.system(cmd)
  #diffData = netCDF4.Dataset(diffFile,'a')
  
  fieldKeys = ['uReconstructX','uReconstructY','uReconstructZ','uReconstructZonal', 'uReconstructMeridional','w','theta','qv','qc','qr']
  for i in xrange(nTimeInds):
  #for i in xrange(4,5,1):
    t = t0+i*h6
    initFile = loc_initFile(t)
    initData = netCDF4.Dataset(initFile,'r')
    
    for key in fieldKeys:
      varInit = initData.variables[key][0,:,:]
      varMPAS = mpasData.variables[key][i,:,:]
      varMPAS -= varInit
      s = "Statistics for MPAS-Init for variable {0} at timeInd {1}\n".format(key,i)
      m = np.mean(varMPAS); dev = np.std(varMPAS); maxd = np.max(np.abs(varMPAS))
      s += "{0},\t{1},\t{2}\n".format(m,dev,maxd)
      print s; s = ''
      #diffData.variables[key][i,:,:] = varMPAS[:,:]

      #histograms
      varMPAS.flatten() #else each column is separate dataset
      nBins = 100
      plt2.hist(varMPAS, bins=nBins, normed=False)
      #save instead of displaying
      figName = key+'.t'+str(i)
      plt2.title(figName)
      #plt2.savefig(figName+'.png')
      plt2.autoscale(enable=True, axis='x', tight=True)
      plt2.savefig(figName+'.png')
      plt2.clf()
    
    initData.close()
    
  mpasData.close(); #diffData.close();

def example_difference_vtk():
  #For this, we have mpas outputs and CFSR ICs every 6 hours.
  #The files are on the same horizontal and vertical mesh
  
  t0 = dt.datetime(2006,7,15,0)
  mpasFile = mpasOutputName(t0)
  mpasData = netCDF4.Dataset(mpasFile,'r')
  
  #fields: config_frames_per_outfile tells us when we need to get new file
  #config_output_interval time separation between values
  nTimeInds = 28
  h6 = dt.timedelta(hours=6)
  
  s = cfsr.form_cfsrTimeString(t0)
  diffFile = 'diff.'+s+'.6hr.vtk'
  
  fvtk = output_data.write_vtk_header_polydata(diffFile, s)
  
  #write nodes and cells
  nNodes =  output_data.write_vtk_xyzNodes(fvtk, mpasData)
  nCells = output_data.write_vtk_polygons(fvtk, mpasData)
  nLevels =  len(mpasData.dimensions['nVertLevels'])
  
  fieldKeys = ['uReconstructX','uReconstructY','uReconstructZ','uReconstructZonal', 'uReconstructMeridional','w','theta']#,'qv','qc','qr']
  for i in xrange(nTimeInds):
  #for i in xrange(4,5,1):
    t = t0+i*h6
    initFile = loc_initFile(t)
    initData = netCDF4.Dataset(initFile,'r')
    
    for key in fieldKeys:
      varInit = initData.variables[key][0,:,:]
      varMPAS = mpasData.variables[key][i,:,:]
      varMPAS -= varInit
      s = "Statistics for MPAS-Init for variable {0} at timInd {1}\n".format(key,i)
      m = np.mean(varMPAS); dev = np.std(varMPAS); maxd = np.max(np.abs(varMPAS))
      s += "{0},\t{1},\t{2}\n".format(m,dev,maxd)
      print s; s = ''
      #diffData.variables[key][i,:,:] = varMPAS[:,:]

      #histograms
      varMPAS.flatten() #else each column is separate dataset
      nBins = 100
      plt2.hist(varMPAS, bins=nBins, normed=False)
      #save instead of displaying
      figName = key+'.t'+str(i)
      plt2.title(figName)
      #plt2.savefig(figName+'.png')
      plt2.autoscale(enable=True, axis='x', tight=True)
      plt2.savefig(figName+'.png')
      plt2.clf()
    
    initData.close()
    
  mpasData.close(); #diffData.close();
  
def example_runDifferences():
  #import netCDF4
  #import datetime as dt
  #from time import strftime
  #import numpy as np
  
  path1 = '/arctic1/nick/cases/163842/testDuda/'
  path2 = '/arctic1/mduda/60km/'
  
  t0 = dt.datetime(2006,7,8,12)
  deltat = dt.timedelta(hours=12)
  tFinal = dt.datetime(2006,7,8,12); #tFinal=t0
  t=t0
  while (t<=tFinal):
    #construct filename of output: eg x1.163842.output.2006-07-13_12.00.00.nc
    tTuple = dt.datetime.timetuple(t);
    s = strftime('%Y-%m-%d_%H.%M.%S', tTuple)
    fname = 'x1.163842.output.'+s+'.nc'
    #fname = x1.163842.static.nc
    
    #load netcdf files
    print "Differencing files for t="+str(t)+'\n'
    data1 = netCDF4.Dataset(path1+fname,'r')
    data2 = netCDF4.Dataset(path2+fname,'r')

    fieldKeys = ['uReconstructX','uReconstructY','uReconstructZ','uReconstructZonal', 'uReconstructMeridional','w','theta','qv','qc','qr']
    #for v in data1.variables:
    for v in fieldKeys:
      try:
        diff = data1.variables[v][:]-data2.variables[v][:]
        maxDiff = np.max(np.abs(diff))
        print "Maxdiff for {0} is {1}".format(v,maxDiff) #print gives a newline
      except Exception:
        print "Couldn't difference variable {0}".format(v)
        continue
      #if (maxDiff>1.e-10):
        #print "Variable {0} has maxDiff={1} between runs\n".format(v,maxDiff)
    
    data1.close(); data2.close()
    t += deltat

def parseStatisticsSection(l):
  #given a line written as 'mean, stdev, max', return these
  s = l.strip().split(',')
  vals = []
  for val in s:
    vals.append(float(val))
  return vals

def parseStatisticsLog(fname,key):
  #A section in the log will generally look like:
  #Statistics for MPAS-Init for variable uReconstructMeridional at timInd 1
  #0.101107448641,	1.57482784472,	24.2203702459
  
  f = open(fname,'r')
  allLines = f.readlines()
  nLines = len(allLines)
  meand = []; stdevd = []; maxd = [];
  searchText = 'variable '+key
  for l in xrange(nLines):
    if (searchText in allLines[l]):
      vals = parseStatisticsSection(allLines[l+1])
      meand.append(vals[0]); stdevd.append(vals[1]); maxd.append(vals[2])
  
  f.close()
  return (meand,stdevd,maxd)

def plotStatistics_example(var):
  import matplotlib.pyplot as plt
  
  stat = parseStatisticsLog('./ComparePlots/stat.txt',var)
  plt.subplot(311) #311 is 3 rows, 1 column, first figure
  plt.plot(stat[0]);
  plt.xlabel('Time index (6 hours)')
  #plt.ylabel('Kelvin')
  plt.ylabel('Variable unit')
  plt.title('Mean of {0} difference'.format(var))
  
  plt.subplot(312) #311 is 3 rows, 1 column, first figure
  plt.plot(stat[1]);
  plt.xlabel('Time index (6 hours)')
  #plt.ylabel('Kelvin')
  plt.title('Standard deviation of {0} difference'.format(var))
  
  plt.subplot(313) #311 is 3 rows, 1 column, first figure
  plt.plot(stat[2]);
  plt.xlabel('Time index (6 hours)')
  #plt.ylabel('Kelvin')
  plt.title('Maximum of absolute {0} difference'.format(var))
  
  plt.tight_layout(pad=1.08)#, h_pad=None, w_pad=None, rect=None) #avoiding overlap of titles and labels
  plt.show()


def calc_height_theta_2PVU(data, t0):
  #vertical height of sign(lat)*2PVU surface.
  #return height of cell center in column using interpolation from top

  pvuTrop = 2.0

  nCells = len(data.dimensions['nCells'])
  nLevels = len(data.dimensions['nVertLevels'])
  epv = data.variables['ertel_pv'][t0,:,:]
  latCell = data.variables['latCell'][:]
  zgrid = data.variables['zgrid'][:,:]
  theta = data.variables['theta'][t0,:,:]

  #hColumn = np.empty(nLevels,dtype='float')
  epv_ht = np.empty(nCells, dtype='float')
  theta_trop = np.empty(nCells,dtype='float')
  for i in xrange(nCells):
    #calc ht of cell centers
    #for k in xrange(nLevels):
      #hColumn[k] = .5*(zgrid[i,k]+zgrid[i,k+1])
    hColumn = .5*(zgrid[i,0:nLevels]+zgrid[i,1:nLevels+1])
    pvuVal = pvuTrop
    if (latCell[i]<0):
      pvuVal = -pvuTrop

    #don't really trust top and bottom levels so don't count those for interpolation
    #interpLevs = range(0,nLevels); nInterpLevs = len(interpLevs)
    interpLevs = range(1,nLevels-1); nInterpLevs = len(interpLevs)
    (l,dl) = output_data.calcIndexOfValue(pvuVal,epv[i,interpLevs], nInterpLevs)
    #(l,dl) = output_data.calcIndexOfValue_fromBottom(pvuVal,epv[i,interpLevs], nInterpLevs)
    #(l,dl) = output_data.calcIndexOfValue(pvuTrop,np.abs(epv[i,interpLevs]), nInterpLevs)
    #print "Cell {0} has l,dl = {1},{2}".format(i, l, dl)
    epv_ht[i] = output_data.calcValueOfIndex(l,dl,hColumn[interpLevs])
    theta_trop[i] = output_data.calcValueOfIndex(l,dl,theta[i,interpLevs])

  return (epv_ht, theta_trop)

def example_2pvu():
  #fpath = '/arctic1/nick/cases/v1.0/x4/august/kf/v1.1/x4.kf.output.2006-08-01_00.00.00.nc'
  fpath = '/arctic1/nick/cases/v1.0/x4/august/tiedtke/v1.1/x4.t.output.2006-08-08_00.00.00.nc'
  #fpath = '/arctic1/nick/cases/v1.0/x4/august/kf/v1.1/x4.kf.output.2006-08-08_00.00.00.nc'
  #fnames = searchFiles()
  #for iFile, fpath in enumerate(fnames):
  data = output_data.open_netcdf_data(fpath)

  for timeInd in xrange(0,28,4):
  #for timeInd in [0]:
    #open the output vtk file and write header, nodes, and cells
    vtkfname = 'x4_t_2006-08-01_1day.'+str(28+timeInd)+'.vtk'
    #vtkfname = 'x4_cfsr_2006-07-25_1day.'+str(iFile)+'.vtk'
    fvtk = output_data.write_vtk_header_polydata(vtkfname, fpath)
    #fvtk = open(vtkfname,'w'); nCells = 163842
    nNodes =  output_data.write_vtk_xyzNodes(fvtk, data)
    nCells = output_data.write_vtk_polygons(fvtk, data)

    #write some cell dataa
    fvtk.write('\nCELL_DATA '+str(nCells)+'\n')

    #calc values
    #timeInd = 0
    epv_ht, theta_trop = calc_height_theta_2PVU(data, timeInd)

    output_data.write_levelData_float('ht_2pvu', fvtk, epv_ht,nCells)
    output_data.write_levelData_float('theta_2pvu', fvtk, theta_trop,nCells)

    fvtk.close()
  data.close()

def searchFiles():
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
/arctic1/nick/cases/cfsr/output/x4.cfsr.output.2006-08-09_00.00.00.nc
/arctic1/nick/cases/cfsr/output/x4.cfsr.output.2006-08-10_00.00.00.nc
/arctic1/nick/cases/cfsr/output/x4.cfsr.output.2006-08-11_00.00.00.nc
/arctic1/nick/cases/cfsr/output/x4.cfsr.output.2006-08-12_00.00.00.nc
/arctic1/nick/cases/cfsr/output/x4.cfsr.output.2006-08-13_00.00.00.nc
/arctic1/nick/cases/cfsr/output/x4.cfsr.output.2006-08-14_00.00.00.nc
/arctic1/nick/cases/cfsr/output/x4.cfsr.output.2006-08-15_00.00.00.nc
'''.split()
  return s

if __name__=='__main__':
  #example_difference()
  example_2pvu()
  #example_stat_files()
