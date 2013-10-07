import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import datetime as dt

import compare
import output_data
import sys; sys.path.append("/data01/tracks/");
import cfsrStats
import mpasStats

def calc_distSphere_multiple(r, lat1, lon1, lat2, lon2):
  #input angles in what sine, cosine take...radians apparently
  
  dlat = lat2-lat1
  dlon = lon2-lon1
  
  latTerm = np.sin(.5*dlat)
  latTerm *= latTerm
  lonTerm = np.sin(.5*dlon); lonTerm = lonTerm*lonTerm*np.cos(lat1)*np.cos(lat2)
  
  dAngle = np.sqrt(latTerm+lonTerm)
  dVals = 2.*r*np.arcsin(dAngle)
  
  return(dVals)

def listStringsToFloat(s):
  a = [float(i) for i in s]
  return a

def get_commonIndices_time(dateList1, dateList2):
  #indices and values of common time between 2 datelists
  ind1 = [i for i,d in enumerate(dateList1) if d in dateList2]
  ind2 = [i for i,d in enumerate(dateList2) if d in dateList1]
  dates = [dateList1[i] for i in ind1]
  return (ind1, ind2, dates)

def plotTracks_dvals_time():
  
  tBase = dt.datetime(2006,8,1,0)
  
  cellInds_cfsr = cfsrStats.x4_cfsr_cell
  
  nSteps = len(cellInds_cfsr)
  dateList_cfsr = [ tBase + dt.timedelta(hours=12*x) for x in range(0,nSteps) ]
  
  vals_cfsr = np.array(cfsrStats.x4_cfsr_minVal)
  vals_t = np.array(mpasStats.x4_t_vals)
  vals_kf = np.array(mpasStats.x4_kf_vals)
  
  cellInds_t = mpasStats.x4_t_cells
  cellInds_kf = mpasStats.x4_kf_cells
  
  nSteps = len(cellInds_t)
  dateList_x4 = [ tBase + dt.timedelta(hours=6*x) for x in range(0,nSteps) ]
  
  #
  ind_cfsr = [i for i,d in enumerate(dateList_cfsr) if d in dateList_x4]
  ind_x4 = [i for i,d in enumerate(dateList_x4) if d in dateList_cfsr]
  dates_common_x4 = [dateList_x4[i] for i in ind_x4]
  
  vals_x7_t = np.array(mpasStats.x7_t_vals)
  vals_x7_kf = np.array(mpasStats.x7_kf_vals)
  nSteps = len(vals_x7_t)
  dateList_x7 = [ tBase + dt.timedelta(hours=6*x) for x in range(0,nSteps) ]
  nSteps = len(vals_x7_kf)
  dateList_x7_kf = [ tBase + dt.timedelta(hours=6*x) for x in range(0,nSteps) ]
  
  ind_cfsr_x7 = [i for i,d in enumerate(dateList_cfsr) if d in dateList_x7]
  ind_x7 = [i for i,d in enumerate(dateList_x7) if d in dateList_cfsr]
  dates_common_x7_t = [dateList_x7[i] for i in ind_x7]
  
  ind_cfsr_x7_kf = [i for i,d in enumerate(dateList_cfsr) if d in dateList_x7_kf]
  ind_x7_kf = [i for i,d in enumerate(dateList_x7_kf) if d in dateList_cfsr]
  dates_common_x7_kf = [dateList_x7_kf[i] for i in ind_x7_kf]
  
  plt.figure() #---------------------------
  plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%m/%d')) #('%Y/%m/%d/%H'))
  plt.gca().xaxis.set_major_locator(mdates.DayLocator(interval=4))
  plt.gca().xaxis.set_minor_locator(mdates.DayLocator(interval=1))
  
  #use the same line properties for scheme as set in this section
  #x4 2006-08-01_00
  plt.plot(dates_common_x4,vals_t[ind_x4]-vals_cfsr[ind_cfsr], 'b', label='x4_t')
  plt.plot(dates_common_x4,vals_kf[ind_x4]-vals_cfsr[ind_cfsr], 'g', label='x4_kf')
  #x7 2006-08-01_00
  plt.plot(dates_common_x7_t,vals_x7_t[ind_x7]-vals_cfsr[ind_cfsr_x7], 'r.-', label='x7_t')
  plt.plot(dates_common_x7_kf,vals_x7_kf[ind_x7_kf]-vals_cfsr[ind_cfsr_x7_kf], 'k--', label='x7_kf')
  
  #x7 2006-08-07_18
  x7_kf_0807_vals = np.array(mpasStats.x7_kf_0807_vals)
  nSteps = len(x7_kf_0807_vals)
  tBase = dt.datetime(2006,8,7,18)
  dateList_x = [ tBase + dt.timedelta(hours=6*x) for x in range(0,nSteps) ]
  
  ind_cfsr = [i for i,d in enumerate(dateList_cfsr) if d in dateList_x]
  ind_x = [i for i,d in enumerate(dateList_x) if d in dateList_cfsr]
  dates_common = [dateList_x[i] for i in ind_x]
  plt.plot(dates_common,x7_kf_0807_vals[ind_x]-vals_cfsr[ind_cfsr], 'k--')
  
  #x7 2006-08-03_12
  x_vals = np.array(mpasStats.x7_t_0803_vals)
  nSteps = len(x_vals)
  tBase = dt.datetime(2006,8,3,12)
  dateList_x = [ tBase + dt.timedelta(hours=6*x) for x in range(0,nSteps) ]
  
  ind_cfsr = [i for i,d in enumerate(dateList_cfsr) if d in dateList_x]
  ind_x = [i for i,d in enumerate(dateList_x) if d in dateList_cfsr]
  dates_common = [dateList_x[i] for i in ind_x]
  plt.plot(dates_common,x_vals[ind_x]-vals_cfsr[ind_cfsr], 'r.-')
  
  #x4_t 2006081500
  x_vals = np.array(mpasStats.x4_t_0815_vals)
  nSteps = len(x_vals)
  tBase = dt.datetime(2006,8,15,0)
  dateList_x = [ tBase + dt.timedelta(hours=6*x) for x in range(0,nSteps) ]
  
  ind_cfsr = [i for i,d in enumerate(dateList_cfsr) if d in dateList_x]
  ind_x = [i for i,d in enumerate(dateList_x) if d in dateList_cfsr]
  dates_common = [dateList_x[i] for i in ind_x]
  
  plt.plot(dates_common,x_vals[ind_x]-vals_cfsr[ind_cfsr], 'b')
  
  #all cfsr ------------------------
  #cells
  c_all_cells = np.array(cfsrStats.x4_cfsr_july_cells+cfsrStats.x4_cfsr_cell+cfsrStats.x4_cfsr_0919_cells)
  tBase = dt.datetime(2006,7,24,12)
  nSteps = len(c_all_cells)-len(cfsrStats.x4_cfsr_0919_cells)
  dateList_c_all = [ tBase + dt.timedelta(hours=12*x) for x in range(0,nSteps) ]
  tBase = dt.datetime(2006,9,19,0)
  nSteps = len(cfsrStats.x4_cfsr_0919_cells)
  dateList_c_all = dateList_c_all + [ tBase + dt.timedelta(hours=12*x) for x in range(0,nSteps) ]
  #values
  c_all_vals = cfsrStats.x4_cfsr_july_vals+cfsrStats.x4_cfsr_minVal+cfsrStats.x4_cfsr_0919_vals
  c_all_vals = np.array(c_all_vals)
  
  #x4_t 2006072412
  x_vals = np.array(mpasStats.x4_t_0724_vals)
  nSteps = len(x_vals)
  tBase = dt.datetime(2006,7,24,12)
  dateList_x = [ tBase + dt.timedelta(hours=6*x) for x in range(0,nSteps) ]
  
  ind_cfsr = [i for i,d in enumerate(dateList_c_all) if d in dateList_x]
  ind_x = [i for i,d in enumerate(dateList_x) if d in dateList_c_all]
  dates_common = [dateList_x[i] for i in ind_x]
  print x_vals[ind_x]
  print c_all_vals[ind_cfsr]
  print dates_common
  plt.plot(dates_common,x_vals[ind_x]-c_all_vals[ind_cfsr], 'b')
  
  #x4_t_0919_vals
  x_vals = np.array(mpasStats.x4_t_0919_vals)
  nSteps = len(x_vals)
  tBase = dt.datetime(2006,9,19,0)
  dateList_x = [ tBase + dt.timedelta(hours=6*x) for x in range(0,nSteps) ]
  
  ind_cfsr = [i for i,d in enumerate(dateList_c_all) if d in dateList_x]
  ind_x = [i for i,d in enumerate(dateList_x) if d in dateList_c_all]
  dates_common = [dateList_x[i] for i in ind_x]
  
  plt.plot(dates_common,x_vals[ind_x]-c_all_vals[ind_cfsr], 'b')
  
  plt.ylabel('$\Delta\Theta_{min}$ , K')
  
  plt.gcf().autofmt_xdate()
  
  #plotOrder = ['x4_t', 'x4_kf', 'x7_t', 'x7_kf']
  #plt.legend(plotOrder)
  plt.legend()
  
  plt.show()

def plotTracks_lat_time():
  #how to show longitude is a vexing question since it wraps -180,180 or 360,0...
  
  fpath = '/arctic1/nick/cases/cfsr/output/x4.cfsr.output.2006-08-01_00.00.00.nc'
  data = output_data.open_netcdf_data(fpath)
  latCell_x4 = data.variables['latCell'][:]*180./np.pi
  data.close()
  fpath = '/arctic1/nick/cases/vduda/x7/x7.kf.output.2006-08-07_18.00.00.nc'
  data = output_data.open_netcdf_data(fpath)
  latCell_x7 = data.variables['latCell'][:]*180./np.pi
  data.close()
  
  #cfsr
  cells_cfsr = cfsrStats.x4_cfsr_cell
  lat_cfsr = latCell_x4[cells_cfsr]
  tBase = dt.datetime(2006,8,1,0)
  nSteps = len(cells_cfsr)
  dateList_cfsr = [ tBase + dt.timedelta(hours=12*x) for x in range(0,nSteps) ]
  
  plt.figure()
  plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%m/%d'))
  plt.gca().xaxis.set_major_locator(mdates.DayLocator(interval=2))
  plt.gca().xaxis.set_minor_locator(mdates.DayLocator(interval=1))
  
  #2006-08-01
  cells_x = mpasStats.x4_t_cells
  lat_x = latCell_x4[cells_x]
  tBase = dt.datetime(2006,8,1,0)
  nSteps = len(cells_x)
  dateList_x = [ tBase + dt.timedelta(hours=6*x) for x in range(0,nSteps) ]
  
  ind_cfsr, ind_x, dates_common = get_commonIndices_time(dateList_cfsr, dateList_x)
  plt.plot(dates_common,lat_x[ind_x]-lat_cfsr[ind_cfsr], 'b', label='x4_t')
  
  cells_x = mpasStats.x4_kf_cells
  lat_x = latCell_x4[cells_x]
  tBase = dt.datetime(2006,8,1,0)
  nSteps = len(cells_x)
  dateList_x = [ tBase + dt.timedelta(hours=6*x) for x in range(0,nSteps) ]
  
  ind_cfsr, ind_x, dates_common = get_commonIndices_time(dateList_cfsr, dateList_x)
  plt.plot(dates_common,lat_x[ind_x]-lat_cfsr[ind_cfsr], 'g', label='x4_kf')
  
  cells_x = mpasStats.x7_t_cells
  lat_x = latCell_x7[cells_x]
  tBase = dt.datetime(2006,8,1,0)
  nSteps = len(cells_x)
  dateList_x = [ tBase + dt.timedelta(hours=6*x) for x in range(0,nSteps) ]
  
  ind_cfsr, ind_x, dates_common = get_commonIndices_time(dateList_cfsr, dateList_x)
  plt.plot(dates_common,lat_x[ind_x]-lat_cfsr[ind_cfsr], 'r.-', label='x7_t')
  
  cells_x = mpasStats.x7_kf_cells
  lat_x = latCell_x7[cells_x]
  tBase = dt.datetime(2006,8,1,0)
  nSteps = len(cells_x)
  dateList_x = [ tBase + dt.timedelta(hours=6*x) for x in range(0,nSteps) ]
  
  ind_cfsr, ind_x, dates_common = get_commonIndices_time(dateList_cfsr, dateList_x)
  plt.plot(dates_common,lat_x[ind_x]-lat_cfsr[ind_cfsr], 'k--', label='x7_kf')
  
  #x7 2006-08-07_18
  cells_x = mpasStats.x7_kf_0807_cells
  lat_x = latCell_x7[cells_x]
  tBase = dt.datetime(2006,8,7,18)
  nSteps = len(cells_x)
  dateList_x = [ tBase + dt.timedelta(hours=6*x) for x in range(0,nSteps) ]
  
  ind_cfsr, ind_x, dates_common = get_commonIndices_time(dateList_cfsr, dateList_x)
  plt.plot(dates_common,lat_x[ind_x]-lat_cfsr[ind_cfsr], 'k--')
  
  #x7 2006-08-03_12
  cells_x = mpasStats.x7_t_0803_cells
  lat_x = latCell_x7[cells_x]
  tBase = dt.datetime(2006,8,3,12)
  nSteps = len(cells_x)
  dateList_x = [ tBase + dt.timedelta(hours=6*x) for x in range(0,nSteps) ]
  
  ind_cfsr, ind_x, dates_common = get_commonIndices_time(dateList_cfsr, dateList_x)
  plt.plot(dates_common,lat_x[ind_x]-lat_cfsr[ind_cfsr], 'r.-')
  
  #x4_t 2006081500
  cells_x = mpasStats.x4_t_0815_cells
  lat_x = latCell_x4[cells_x]
  tBase = dt.datetime(2006,8,15,0)
  nSteps = len(cells_x)
  dateList_x = [ tBase + dt.timedelta(hours=6*x) for x in range(0,nSteps) ]
  
  ind_cfsr, ind_x, dates_common = get_commonIndices_time(dateList_cfsr, dateList_x)
  plt.plot(dates_common,lat_x[ind_x]-lat_cfsr[ind_cfsr], 'b')
  
  
  plt.ylabel('$\Delta\phi_{min}$'+', ' +'$\deg$')
  
  plt.gcf().autofmt_xdate()
  
  #plotOrder = ['x4_t', 'x4_kf', 'x7_t', 'x7_kf']
  #plt.legend(plotOrder)
  plt.legend(loc=3)
  
  plt.show()

def plotTracks_dist_time():
  #how to show longitude is a vexing question since it wraps -180,180 or 360,0...
  
  rEarth = 6371. #km
  
  fpath = '/arctic1/nick/cases/cfsr/output/x4.cfsr.output.2006-08-01_00.00.00.nc'
  data = output_data.open_netcdf_data(fpath)
  latCell_x4 = data.variables['latCell'][:]
  lonCell_x4 = data.variables['lonCell'][:]
  data.close()
  fpath = '/arctic1/nick/cases/vduda/x7/x7.kf.output.2006-08-07_18.00.00.nc'
  data = output_data.open_netcdf_data(fpath)
  latCell_x7 = data.variables['latCell'][:]
  lonCell_x7 = data.variables['lonCell'][:]
  data.close()
  
  #cfsr
  cells_cfsr = cfsrStats.x4_cfsr_cell
  lat_cfsr = latCell_x4[cells_cfsr]
  lon_cfsr = lonCell_x4[cells_cfsr]
  tBase = dt.datetime(2006,8,1,0)
  nSteps = len(cells_cfsr)
  dateList_cfsr = [ tBase + dt.timedelta(hours=12*x) for x in range(0,nSteps) ]
  
  plt.figure()
  plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%m/%d'))
  plt.gca().xaxis.set_major_locator(mdates.DayLocator(interval=4))
  plt.gca().xaxis.set_minor_locator(mdates.DayLocator(interval=1))
  
  #2006-08-01
  cells_x = mpasStats.x4_t_cells
  lat_x = latCell_x4[cells_x]
  lon_x = lonCell_x4[cells_x]
  tBase = dt.datetime(2006,8,1,0)
  nSteps = len(cells_x)
  dateList_x = [ tBase + dt.timedelta(hours=6*x) for x in range(0,nSteps) ]
  
  ind_cfsr, ind_x, dates_common = get_commonIndices_time(dateList_cfsr, dateList_x)
  
  vals = calc_distSphere_multiple(rEarth, lat_x[ind_x], lon_x[ind_x], lat_cfsr[ind_cfsr], lon_cfsr[ind_cfsr])
  plt.plot(dates_common,vals, 'b', label='x4_t')
  
  cells_x = mpasStats.x4_kf_cells
  lat_x = latCell_x4[cells_x]
  lon_x = lonCell_x4[cells_x]
  tBase = dt.datetime(2006,8,1,0)
  nSteps = len(cells_x)
  dateList_x = [ tBase + dt.timedelta(hours=6*x) for x in range(0,nSteps) ]
  
  ind_cfsr, ind_x, dates_common = get_commonIndices_time(dateList_cfsr, dateList_x)
  vals = calc_distSphere_multiple(rEarth, lat_x[ind_x], lon_x[ind_x], lat_cfsr[ind_cfsr], lon_cfsr[ind_cfsr])
  plt.plot(dates_common,vals, 'g', label='x4_kf')
  
  cells_x = mpasStats.x7_t_cells
  lat_x = latCell_x7[cells_x]
  lon_x = lonCell_x7[cells_x]
  tBase = dt.datetime(2006,8,1,0)
  nSteps = len(cells_x)
  dateList_x = [ tBase + dt.timedelta(hours=6*x) for x in range(0,nSteps) ]
  
  ind_cfsr, ind_x, dates_common = get_commonIndices_time(dateList_cfsr, dateList_x)
  vals = calc_distSphere_multiple(rEarth, lat_x[ind_x], lon_x[ind_x], lat_cfsr[ind_cfsr], lon_cfsr[ind_cfsr])
  plt.plot(dates_common,vals, 'r.-', label='x7_t')
  
  cells_x = mpasStats.x7_kf_cells
  lat_x = latCell_x7[cells_x]
  lon_x = lonCell_x7[cells_x]
  tBase = dt.datetime(2006,8,1,0)
  nSteps = len(cells_x)
  dateList_x = [ tBase + dt.timedelta(hours=6*x) for x in range(0,nSteps) ]
  
  ind_cfsr, ind_x, dates_common = get_commonIndices_time(dateList_cfsr, dateList_x)
  vals = calc_distSphere_multiple(rEarth, lat_x[ind_x], lon_x[ind_x], lat_cfsr[ind_cfsr], lon_cfsr[ind_cfsr])
  plt.plot(dates_common,vals, 'k--', label='x7_kf')
  
  #x7 2006-08-07_18
  cells_x = mpasStats.x7_kf_0807_cells
  lat_x = latCell_x7[cells_x]
  lon_x = lonCell_x7[cells_x]
  tBase = dt.datetime(2006,8,7,18)
  nSteps = len(cells_x)
  dateList_x = [ tBase + dt.timedelta(hours=6*x) for x in range(0,nSteps) ]
  
  ind_cfsr, ind_x, dates_common = get_commonIndices_time(dateList_cfsr, dateList_x)
  vals = calc_distSphere_multiple(rEarth, lat_x[ind_x], lon_x[ind_x], lat_cfsr[ind_cfsr], lon_cfsr[ind_cfsr])
  plt.plot(dates_common,vals, 'k--')
  
  #x7 2006-08-03_12
  cells_x = mpasStats.x7_t_0803_cells
  lat_x = latCell_x7[cells_x]
  lon_x = lonCell_x7[cells_x]
  tBase = dt.datetime(2006,8,3,12)
  nSteps = len(cells_x)
  dateList_x = [ tBase + dt.timedelta(hours=6*x) for x in range(0,nSteps) ]
  
  ind_cfsr, ind_x, dates_common = get_commonIndices_time(dateList_cfsr, dateList_x)
  vals = calc_distSphere_multiple(rEarth, lat_x[ind_x], lon_x[ind_x], lat_cfsr[ind_cfsr], lon_cfsr[ind_cfsr])
  plt.plot(dates_common,vals, 'r.-')
  
  #x4_t 2006081500
  cells_x = mpasStats.x4_t_0815_cells
  lat_x = latCell_x4[cells_x]
  lon_x = lonCell_x4[cells_x]
  tBase = dt.datetime(2006,8,15,0)
  nSteps = len(cells_x)
  dateList_x = [ tBase + dt.timedelta(hours=6*x) for x in range(0,nSteps) ]
  
  ind_cfsr, ind_x, dates_common = get_commonIndices_time(dateList_cfsr, dateList_x)
  vals = calc_distSphere_multiple(rEarth, lat_x[ind_x], lon_x[ind_x], lat_cfsr[ind_cfsr], lon_cfsr[ind_cfsr])
  plt.plot(dates_common,vals, 'b')
  
  #all cfsr ------------------------
  #cells
  cells_cfsr = np.array(cfsrStats.x4_cfsr_july_cells+cfsrStats.x4_cfsr_cell+cfsrStats.x4_cfsr_0919_cells)
  lat_cfsr = latCell_x4[cells_cfsr]
  lon_cfsr = lonCell_x4[cells_cfsr]
  
  tBase = dt.datetime(2006,7,24,12)
  nSteps = len(cells_cfsr)-len(cfsrStats.x4_cfsr_0919_cells)
  dateList_c_all = [ tBase + dt.timedelta(hours=12*x) for x in range(0,nSteps) ]
  tBase = dt.datetime(2006,9,19,0)
  nSteps = len(cfsrStats.x4_cfsr_0919_cells)
  
  dateList_cfsr = dateList_c_all + [ tBase + dt.timedelta(hours=12*x) for x in range(0,nSteps) ]
  
  #x4_t 2006072412
  cells_x = mpasStats.x4_t_0724_cells
  lat_x = latCell_x4[cells_x]
  lon_x = lonCell_x4[cells_x]
  nSteps = len(cells_x)
  tBase = dt.datetime(2006,7,24,12)
  dateList_x = [ tBase + dt.timedelta(hours=6*x) for x in range(0,nSteps) ]
  
  ind_cfsr, ind_x, dates_common = get_commonIndices_time(dateList_cfsr, dateList_x)
  vals = calc_distSphere_multiple(rEarth, lat_x[ind_x], lon_x[ind_x], lat_cfsr[ind_cfsr], lon_cfsr[ind_cfsr])
  plt.plot(dates_common,vals, 'b')
  
  #x4_t_0919_vals
  cells_x = mpasStats.x4_t_0919_cells
  lat_x = latCell_x4[cells_x]
  lon_x = lonCell_x4[cells_x]
  nSteps = len(cells_x)
  tBase = dt.datetime(2006,9,19,0)
  dateList_x = [ tBase + dt.timedelta(hours=6*x) for x in range(0,nSteps) ]
  
  ind_cfsr, ind_x, dates_common = get_commonIndices_time(dateList_cfsr, dateList_x)
  vals = calc_distSphere_multiple(rEarth, lat_x[ind_x], lon_x[ind_x], lat_cfsr[ind_cfsr], lon_cfsr[ind_cfsr])
  plt.plot(dates_common,vals, 'b')
  
  plt.ylabel('$|\Delta x|$'+', km')
  
  plt.gcf().autofmt_xdate()
  
  #plotOrder = ['x4_t', 'x4_kf', 'x7_t', 'x7_kf']
  #plt.legend(plotOrder)
  #plt.legend(loc=3) #bottom left
  plt.legend()
  
  plt.show()

def dataOutput_augTheta():
  x4_t_vals = '''
301.4533597
299.3297445
297.3610777
295.6199052
294.1198904
292.7884307
291.5967839
290.2545776
289.4173524
301.5380895
309.5461473
306.6020387
308.7006765
316.4247718
'''.split()

  x4_kf_vals = '''
301.4533597
299.3231661
297.3670824
295.614679
294.1068815
292.8219867
291.5510108
290.542331
292.0127833
288.2663763
287.3617802
292.8418938
288.2731879
293.1016675
'''.split()

  x4_cfsr_vals = '''
301.4533597
300.9854183
299.7418007
298.3698262
297.3424742
297.1961196
296.5385434
297.8222468
297.5584998
293.9392789
299.7709242
291.424318
293.2356347
291.5369664
'''.split()

  return (np.array(listStringsToFloat(x4_t_vals)),
          np.array(listStringsToFloat(x4_kf_vals)),
          np.array(listStringsToFloat(x4_cfsr_vals)))

def latLonAugust():
  test = '''
1.54316421	5.1796831
1.51601256	5.81537576
1.48528808	0.01336435
1.51747693	0.07884841
1.5072933	5.24327442
1.48272698	4.93401599
1.50019658	4.46295316
1.45248069	3.77505923
1.4403253	3.8987524
1.42551189	3.33960653
1.50103199	3.28630323
1.44294467	3.5556655
1.45912444	3.86952366
1.50303516	4.0109568
'''.split()
  nVals = len(test)/2
  ll_x4t = np.empty((nVals,2))
  for i in xrange(nVals):
    ll_x4t[i,0] = float(test[i*2+0])
    ll_x4t[i,1] = float(test[i*2+1])

  test = '''
1.54316421	5.1796831
1.51601256	5.81537576
1.48812449	6.27602155
1.51729479	0.01138975
1.50190743	5.17119768
1.4591825	4.79634907
1.39362819	4.62856897
1.38151196	5.05074463
1.44073973	4.72432944
1.42630467	4.46591225
1.49086623	4.86448714
1.42414033	4.48639981
1.48149108	4.57367372
1.41043434	4.82314659
  '''.split()

  nVals = len(test)/2
  ll_x4kf = np.empty((nVals,2))
  for i in xrange(nVals):
    ll_x4kf[i,0] = float(test[i*2+0])
    ll_x4kf[i,1] = float(test[i*2+1])

  test = '''
  1.54316421	5.1796831
1.51884694	6.12215272
1.49367543	6.23067401
1.50865698	0.09336126
1.54236896	5.51816595
1.50192029	4.85472723
1.55342346	0.01263294
1.50820135	2.46813339
1.45092428	3.06426216
1.41299368	3.49703965
1.48935102	3.6475312
1.52895676	2.69768265
1.55411551	0.75387787
1.48483641	1.00666332
1.44186464	0.70839626
  '''.split()

  nVals = len(test)/2 #-1
  ll_x4cfsr = np.empty((nVals,2))
  for i in xrange(nVals):
    ll_x4cfsr[i,0] = float(test[i*2+0])
    ll_x4cfsr[i,1] = float(test[i*2+1])

  return (ll_x4t, ll_x4kf, ll_x4cfsr)

def latLonJuly():
  test = '''
1.212997771	2.845908187
1.16993311	2.941363731
1.195835449	2.966504771
1.258090684	3.11303472
1.445794019	3.599530112
1.510427886	2.732732409
1.53737124	3.601857
1.548991831	5.121989569
1.520217946	6.00979144
1.499463612	6.2305855
1.508350857	6.259060315
1.543715349	5.712447387
1.509549321	4.756857031
1.562092557	6.02665685
1.490738554	2.497406346
'''.split()
  nVals = len(test)/2
  ll_cfsr = np.empty((nVals,2))
  for i in xrange(nVals):
    ll_cfsr[i,0] = float(test[i*2+0])
    ll_cfsr[i,1] = float(test[i*2+1])

  test = '''
1.217388708	2.84812407
1.196031175	2.918323848
1.19358669	2.944479287
1.272006926	3.053166596
1.338225013	3.043120919
1.408908574	2.130969922
1.392540009	2.629249335
1.357181908	2.861174418
1.332163461	2.664581546
1.323544993	2.723711201
1.354443039	2.694047911
1.406944846	2.692151811
1.472060841	2.798806627
1.502971049	2.093953831
  '''.split()

  nVals = len(test)/2
  ll_x1kf = np.empty((nVals,2))
  for i in xrange(nVals):
    ll_x1kf[i,0] = float(test[i*2+0])
    ll_x1kf[i,1] = float(test[i*2+1])

  test = '''
1.217388708	2.84812407
1.196031175	2.918323848
1.199675892	2.96454585
1.272006926	3.053166596
1.375618006	3.185763156
1.395103118	2.168086216
1.347238611	2.668686803
1.347584442	2.623826232
1.263430148	2.667370143
1.300407551	2.732333609
1.19358669	2.944479287
1.340858492	3.192108272
1.369424235	3.05861464
1.352619493	3.508640038

  '''.split()

  nVals = len(test)/2 #-1
  ll_x1t = np.empty((nVals,2))
  for i in xrange(nVals):
    ll_x1t[i,0] = float(test[i*2+0])
    ll_x1t[i,1] = float(test[i*2+1])

  test = '''
1.212997771	2.845908187
1.196570791	2.900211515
1.191673736	2.9388886
1.261087794	3.04748311
1.382640281	3.266787733
1.473947153	4.315275514
1.418792365	5.323086978
1.338169326	5.414997452
1.317451217	5.206299118
1.273662041	5.115835832
1.263895226	5.056073121
1.271305982	4.945580364
1.223976209	4.910642253
1.231420445	4.839901418
  '''.split()

  nVals = len(test)/2 #-1
  ll_x4kf = np.empty((nVals,2))
  for i in xrange(nVals):
    ll_x4kf[i,0] = float(test[i*2+0])
    ll_x4kf[i,1] = float(test[i*2+1])

  test = '''
1.216457828	2.844462006
1.196570791	2.900211515
1.199463713	2.966212767
1.265498372	3.025425506
1.385609849	3.19624907
1.501141506	4.239730212
1.441810769	5.692487932
1.354588004	6.021231672
1.324649654	5.970031265
1.274748275	6.145800806
1.311877473	6.211914166
1.283329736	0.171220193
1.27442924	0.219920789
1.332937845	0.003042817
  '''.split()

  nVals = len(test)/2 #-1
  ll_x4t = np.empty((nVals,2))
  for i in xrange(nVals):
    ll_x4t[i,0] = float(test[i*2+0])
    ll_x4t[i,1] = float(test[i*2+1])

  return (ll_x1t,ll_x1kf,ll_x4t, ll_x4kf, ll_cfsr)

def dataOutputJuly():
  cfsr = '''
  245.3551355
243.9758663
245.0400394
245.2595964
244.8970497
246.0079467
246.4601132
244.7694969
244.6749876
243.5938506
242.8745024
242.0907957
241.7707244
241.227902
  '''.split()

  x1kf = '''
  246.3225296
244.9111884
243.6969094
243.1803623
244.632451
242.4324238
242.3478019
242.3436291
241.7382206
241.1696924
241.3824351
241.3916151
239.9552873
239.1590827
  '''.split()
  
  x1t = '''
  246.3129951
244.8177719
243.4413255
242.958061
244.0290964
241.675864
240.9991272
241.5651499
240.82667
240.5160811
241.8815674
244.7844279
244.7624418
244.4322401
  '''.split()
  
  x4kf = '''
  246.1049343
244.3654778
243.154677
242.3692446
242.2681844
242.2380734
241.5894303
242.475133
243.6147359
242.5476453
241.8654697
240.5047682
240.8522627
241.5634222
  '''.split()

  x4t = '''
  246.0947353
244.2939771
243.027821
242.3491073
242.7045674
241.9908339
240.8363903
240.2901002
239.3016135
239.8764337
241.0056306
240.0258912
240.0636049
240.3516423
  '''.split()

  return (np.array(listStringsToFloat(x1kf)),
          np.array(listStringsToFloat(x1t)),
          np.array(listStringsToFloat(x4kf)),
          np.array(listStringsToFloat(x4t)),
          np.array(listStringsToFloat(cfsr)))

def runPlot_theta():
  #(x4_t, x4_kf, x4_cfsr) = dataOutput()
  (x1kf, x1t, x4kf, x4t, cfsr) = dataOutputJuly()

  dataPlot = [None]*4
  dataPlot[0] = x1kf-cfsr
  dataPlot[1] = x1t-cfsr
  dataPlot[2] = x4kf-cfsr
  dataPlot[3] = x4t-cfsr

  compare.plot_save_lists(dataPlot, 4, ['x1_kf','x1_t', 'x4_kf', 'x4_t'], ' ', 'dtheta_2007-24-12')

def runPlot_ll():

  #ll_x4t, ll_x4kf, ll_x4cfsr = latLonAugust()
  (ll_x1t,ll_x1kf,ll_x4t, ll_x4kf, ll_cfsr) = latLonJuly()
  #print ll_x4t
  #print ll_x4kf
  #print ll_x4cfsr
  nVals = ll_x4t.shape[0]
  dVal = np.empty((4,nVals))
  for i in xrange(nVals):
    dVal[0,i] = conn.calc_distSphere(6371000., ll_x1t[i,:], ll_cfsr[i,:])
    dVal[1,i] = conn.calc_distSphere(6371000., ll_x1kf[i,:], ll_cfsr[i,:])
    dVal[2,i] = conn.calc_distSphere(6371000., ll_x4t[i,:], ll_cfsr[i,:])
    dVal[3,i] = conn.calc_distSphere(6371000., ll_x4kf[i,:], ll_cfsr[i,:])

  dataPlot = [None]*4
  for i in xrange(4):
    dataPlot[i] = dVal[i,:]/1000.
  compare.plot_save_lists(dataPlot, 4, ['x1_t', 'x1_kf','x4_t', 'x4_kf'], ' ', 'dloc_2006-07-24_12')

if __name__ == '__main__':
  #runPlot_theta()
  #runPlot_ll()
  plotTracks_dist_time()
  plotTracks_dvals_time()