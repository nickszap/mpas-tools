#module for creating mpas fields of cfsr data over time

#python libraries
import os
import datetime as dt

#my libraries
import sys; sys.path.append("/raid1/nick/pythonScripts");
import cfsr
#import output_data
import namelist
import arctic

def ungribbed2mpas(t0):
  #given intermediate (ie ungribbed pgb+flx) cfsr and all files in this directory,
  #generate mpas vertical and surface initial conditions
  #integrate for 0 steps to generate fields like pv
  
  np = 12 #number of procs to use
  nNodes = arctic.calc_numNodes(np)
  timeString = cfsr.form_mpasTimeString(t0)
  
  #geography interpolation already done
    
  #vertical and surface
  fname_namelist = "namelist.init"
  try:
    f = open(fname_namelist, 'w') #i don't understand python's scoping. How does f exist outside this block?
  except IOError as e:
    print "I/O error({0}): {1}".format(e.errno, e.strerror)
    pass
  #namelist.writeNamelist_vertSfc_cfsrCase(f, timeString,nNodes,arctic.nCorePerNode)
  namelist.writeNamelist_vertSfc_cfsrCase(f, timeString,nNodes)
  f.close()
  
  os.system('rm namelist.input')
  os.system('ln -s namelist.init namelist.input')
  
  exe = './init_nhyd_atmos_model.exe'
  #arctic.main(exe, np) #no need to mpdboot since go through hydra
  cmd = 'time mpiexec -f ~/machines -n '+str(np)+' '+exe
  os.system(cmd)
  
  '''
  #we don't integrate any more since 'pv' variable is not Ertel's
  #integration
  fname_namelist = "namelist.integrate"
  try:
    f = open(fname_namelist, 'w')
  except IOError as e:
    print "I/O error({0}): {1}".format(e.errno, e.strerror)
    pass
  namelist.writeNamelist_integrate(f, timeString, timeString, nNodes,arctic.nCorePerNode)
  f.close()
  
  os.system('rm namelist.input')
  os.system('ln -s namelist.integrate namelist.input')
  
  exe = './nhyd_atmos_model.exe'
  arctic.main(exe, np)
  '''
  
def example():
  #for every 6 hours, get the atmospheric field from the cfsr 6 hr forecast
  
  t0 = dt.datetime(2006,6,1,0)
  #tf = dt.datetime(2006,9,30,0)
  tf = dt.datetime(2006,6,1,3)
  h6 = dt.timedelta(hours=6)
  
  i=0; t=t0;
  while (t<=tf):
    '''
    #download pgb and flx
    cfsr.cfsr_from_nomads(t.year, t.month, t.day, t.hour)
    
    #make intermediate files through ungrib
    fname_cfsrIntermediate = cfsr.runUngrib_CFSR(t)
    '''
    #generate initial conditions and run for full field
    ungribbed2mpas(t)
    
    #increment day
    i = i+1;
    t = t0+i*h6;
    
  #done w/ all time steps
#

if __name__=='__main__':
  example()