#download june 2006 cfsr data:
#wget -r --no-parent -e robots=off  --reject "index.html*" http://nomads.ncdc.noaa.gov/modeldata/cmd_pgbh/2006/200606/

#format of MPAS time is '2006-07-06_06:00:00':
#strftime('%Y-%m-%d_%H:%M:%S', datetime.datetime.timetuple(datetime.datetime(y,m,d,h)))

import os #instead of wget, python has its own tools too (eg urllib2 with urlretrieve)
import datetime as dt #the 06 forecast is valid six hours from the url date. we need to work in calendar time
from time import strftime, sleep
import netCDF4
import numpy as np

def example_cfsrNetcdf(data):
  #f = "/arctic1/nick/cases/cfsr/2006-anl/pgbhnl.gdas.2006081500.grb2.trop.nc"
  #data = netCDF4.Dataset(f,'r')
  
  lon1d = data.variables['lon_0'][:]
  lat1d = data.variables['lat_0'][:]
  
  nLon = len(lon1d)
  nLat = len(lat1d)


def cfsr_from_nomads(y, m, d, h):
  #given the calendar info input as integers, download the pressure level and surface data for the indicated time.
  #hour should be one of: 00,06,12,18 since we use the 6 hour forecasts since analysis doesn't have flux data

  #For the right time, we think f means forecast valid at url time and h means hour6 of run started at url time.
  #so, use the time we want in flxf06 url and time-6 in pghbh06 url
  
  #form the URLs for pressure and surface data ---------------
  #go through quite a bit of trouble to get the time six hours before
  tf = dt.datetime(y,m,d,h) #forecast time (ie the time of the initial condition)
  h6 = dt.timedelta(hours=6)
  th = tf-h6 #since want hour6 of forecast started at urlTime
  
  #make strings of time info. just str(tf.day) gives 1 instead of 01 for first day
  tf = dt.datetime.timetuple(tf); th = dt.datetime.timetuple(th)
  yf = strftime('%Y',tf)
  ymf = strftime('%Y%m',tf)
  ymdf = strftime('%Y%m%d',tf)
  ymdhf = strftime('%Y%m%d%H',tf)
  
  yh = strftime('%Y',th)
  ymh = strftime('%Y%m',th)
  ymdh = strftime('%Y%m%d',th)
  ymdhh = strftime('%Y%m%d%H',th)
  
  #build the urls as strings
  url_base = 'http://nomads.ncdc.noaa.gov/modeldata/cmd_flxf/'
  subdir = yf+'/'+ymf+'/'+ymdf+'/'
  fnamef = 'flxf06.gdas.'+ymdhf+'.grb2'
  urlf = url_base+subdir+fnamef
  print(urlf)
  
  url_base = 'http://nomads.ncdc.noaa.gov/modeldata/cmd_pgbh/'
  subdir = yh+'/'+ymh+'/'+ymdh+'/'
  fnameh = 'pgbh06.gdas.'+ymdhh+'.grb2'
  urlh = url_base+subdir+fnameh
  print(urlh)
  
  #download the files -----------------------------------
  #os.chdir(datadir) #maybe want it in a different directory
  #print('Pausing to avoid loading nomads server\n'); sleep(15);
  os.system('wget '+urlf)
  sleep(5)
  os.system('wget '+urlh)
  
#end

def time2String(tf):
  #convert the forecast time to a string at this time
  tf = dt.datetime.timetuple(tf)
  yf = strftime('%Y',tf)
  ymf = strftime('%Y%m',tf)
  ymdf = strftime('%Y%m%d',tf)
  ymdhf = strftime('%Y%m%d%H',tf)
  return (yf,ymf,ymdf,ymdhf)

def time2String_6(th):
  #conver the hour6 started six hours ago to a string 6 hours ago
  h6 = dt.timedelta(hours=6)
  th = th-h6
  return time2String(th)
  
def formName_flx6(t0):
  #make the name of the file as it comes from nomads
  a = time2String(t0)
  ymdh = a[3];
  fname = 'flxf06.gdas.'+ymdh+'.grb2'
  return fname

def formName_pgb6(t0):
  #make the name of the file as it comes from nomads
  a = time2String_6(t0)
  ymdh = a[3];
  fname = 'pgbh06.gdas.'+ymdh+'.grb2'
  return fname

def form_cfsrTimeString(t0):
  #format of CFSR time is '2006-07-06_06'
  tTuple = dt.datetime.timetuple(t0);
  s = strftime('%Y-%m-%d_%H', tTuple)
  return s

def form_mpasTimeString(t0):
  #format of mpas time is 2006-07-06_06:00:00'
  tTuple = dt.datetime.timetuple(t0);
  s = strftime('%Y-%m-%d_%H:%M:%S', tTuple)
  return s

def download_july06():
  #get all vertical and surface data from july 2006
  t0 = dt.datetime(2006,7,1,0)
  tf = dt.datetime(2006,7,21,0)
  h6 = dt.timedelta(hours=6)
  i=0; t=t0;
  while (t<=tf):
    t = t0+i*h6
    cfsr_from_nomads(t.year, t.month, t.day, t.hour)
    i = i+1

if __name__=='__main__':
  download_july06()

def makeIntermediate_CFSR(t0, cfsDataPath):
  #Given downloaded CFSR pressure and surface flux data for the given hour,
  #create intermediate WPS files (FLX and CFSR=FLX+PGB) using ungrib
  
  #cfsDataPath = '/raid1/nick/cfs/july'
  os.chdir(cfsDataPath)
  runUngrib_CFSR(t0)
  
def runUngrib_CFSR(t0):
  #For initialization, need FLX and CFSR=FLX+PGB intermediate files from WPS.
  #First make the intermediate PGB and FLX files. Then concatenate the FLX onto the PGB for CFSR.
  #Return the file name of intermediate pgb+flx
  
  #locations of files. everything should already be linked into this directory
  wpsDir = '/raid1/nick/wps/WPS/'
  ungribExe = wpsDir+'ungrib.exe'
  VtablePGB = 'Vtable.CFSR_press_pgbh06'
  VtableFLX = 'Vtable.CFSR_sfc_flxf06'
  
  runtimeString = form_mpasTimeString(t0);
  fnameTime = form_cfsrTimeString(t0);
  fnamePGB = formName_pgb6(t0); prefixPGB='PGB:'
  fnameFLX = formName_flx6(t0); prefixFLX='FLX:'
  
  #do flx one first
  os.system('rm Vtable') #symbolic link
  os.system('ln -s '+VtableFLX+' Vtable')
  os.system('rm GRIBFILE.*')
  os.system('ln -s '+fnameFLX+' GRIBFILE.AAA')
  writeNamelist_wps(runtimeString, runtimeString)
  os.system(ungribExe)
  os.system('mv FILE:'+fnameTime+' '+prefixFLX+fnameTime)
  
  #next is pgb
  os.system('rm Vtable') #symbolic link
  os.system('ln -s '+VtablePGB+' Vtable')
  os.system('rm GRIBFILE.*')
  os.system('ln -s '+fnamePGB+' GRIBFILE.AAA')
  writeNamelist_wps(runtimeString, runtimeString)
  os.system(ungribExe)
  os.system('mv FILE:'+fnameTime+' '+prefixPGB+fnameTime)
  
  #CFSR=pgb+flx
  fnameCFSR = 'CFSR:'+fnameTime
  cmd = 'cat '+prefixPGB+fnameTime+' '+prefixFLX+fnameTime+' '+'> '+fnameCFSR
  os.system(cmd)
  
  #as of now, can delete the nomads files and PGB intermediate
  os.system('rm '+prefixPGB+fnameTime)
  os.system('rm '+fnameFLX); os.system('rm '+fnamePGB);
  
  return(fnameCFSR)
  
def writeNamelist_wps(t0String, tfString):
  #time format is '2006-07-06_06:00:00'
  f = open('namelist.wps','w') #try, except if want to do properly...
  
  s = '''
&share
 start_date = %s,
 end_date = %s,
 interval_seconds = 21600
/

&ungrib
 out_format = 'WPS',
 prefix = 'FILE',
/
''' %(t0String,tfString) #it looks like the prefix gets ignored
  f.write(s)
  
  f.close()
#end

#------------------------------------------------
#we can request variables by hand from: http://rda.ucar.edu/datasets/ds093.0/#description

'''
old stuff:
#download the file
  url_base = 'http://nomads.ncdc.noaa.gov/modeldata/cmd_pgbh/' #pressure levels
  year = '2006'; month = '06'; day = '01'; hour = '00' #(00,06,12,18)
  ym = year+month; ymd = ym+day; ymdh = ymd+hour;
  fname = 'pgbh06.gdas.'+ymdh+'.grb2'
  url = url_base+year+'/'+ym+'/'+ymd+'/'+fname
  #http://nomads.ncdc.noaa.gov/modeldata/cmd_pgbh/2006/200606/20060601/pgbh06.gdas.2006060100.grb2 for pressure level data

  #http://nomads.ncdc.noaa.gov/modeldata/cmd_flxf/2006/200606/20060601/flxf06.gdas.2006060100.grb2 what's in flx vs spl...
  #http://nomads.ncdc.noaa.gov/modeldata/cmd_flxf/2006/200606/20060601/splf06.gdas.2006060100.grb2 for surface data???
  #don't use the pgbh00 files since it's a model spinup. 
  #We use 6hr forecast from nl since all variables are available, as opposed to the nl (analysis) files (without at least the flux/heating vars i think)
  
  #we don't need to rename the gribfile since we relink it to,say, GRIBFILE.AAA for wps
  #rename the h file to be of the time that it's valid
  rightNameh = 'pgbh06.gdas.'+ymdhh+'.grib2'
  os.system('mv ./fnameh ./'+rightNameh)
'''

