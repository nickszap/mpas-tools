import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import netCDF4

def demo_scatter():
  fig = plt.figure()
  ax = fig.add_subplot(111)
  ax.set_title('click on points')

  line, = ax.plot(np.random.rand(100), 'o', picker=5)  # 5 points tolerance

  def onpick(event):
      thisline = event.artist
      xdata = thisline.get_xdata()
      ydata = thisline.get_ydata()
      ind = event.ind
      print 'onpick points:', zip(xdata[ind], ydata[ind])

  fig.canvas.mpl_connect('pick_event', onpick)

  plt.show()

def onpick_map(event):
  thisline = event.artist
  #xdata = thisline.get_xdata()
  #ydata = thisline.get_ydata()
  ind = event.ind
  #print 'onpick points:', zip(xdata[ind], ydata[ind])
  print 'index of pick: ', ind

def demo_map(lat, lon):
  fig = plt.figure()
  ax = fig.add_subplot(111)
  ax.set_title('click on points')
  
  map = Basemap(projection='ortho',lon_0=0,lat_0=90, resolution='l')
  x,y = map(lon, lat)
  
  map.drawcoastlines()
  map.drawmapboundary()
  #map.scatter(x,y, picker=5)
  map.scatter(x,y,s=0.01, picker=1) #s gives area of marker. use small area so can see boundaries
  
  fig.canvas.mpl_connect('pick_event', onpick_map)
  
  plt.show()

def run_map():
  #fName = '/arctic1/nick/cases/v1.0/x4/longer/x4.kf.output.2006-08-15_00.00.00.nc'
  fName = '/arctic1/nick/cases/cfsr/output/x4.cfsr.output.2006-08-01_00.00.00.nc'
  #fName = '/arctic1/mduda/60km/x1.163842.output.2006-07-08_00.00.00.nc'
  data = netCDF4.Dataset(fName,'r')
  
  r2d = 180./np.pi
  latCell = data.variables['latCell'][:]*r2d
  lonCell = data.variables['lonCell'][:]*r2d
  
  demo_map(latCell, lonCell)

if __name__=='__main__':
  #demo()
  run_map()
