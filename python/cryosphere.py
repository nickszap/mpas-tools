#download images from cryosphere today

import os
import datetime as dt
from time import strftime, sleep

def url_time(t0):
  urlBase = 'http://arctic.atmos.uiuc.edu/cryosphere/IMAGES/ARCHIVE/'
  t0 = dt.datetime.timetuple(t0)
  ymd = strftime('%Y%m%d',t0)
  return urlBase+ymd+'.jpg'

def getFiles():
  t0 = dt.datetime(2006,7,7,0)
  tf = dt.datetime(2006,9,30,0)
  deltaT = dt.timedelta(days=1)
  i=0; t=t0;
  while (t<=tf):
    #t = t0+i*deltaT
    s = url_time(t)
    cmd = 'wget '+ s
    os.system(cmd)

    i=i+1
    t = t0+i*deltaT
    sleep(3)

if __name__=='__main__':
  getFiles()

