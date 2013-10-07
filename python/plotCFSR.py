#here we're working on the CFSR lat-lon grid

import numpy as np
import os
import matplotlib.pyplot as plt

def example_trop():
  '''
  tropopause variables are:
  343:59438456:d=1990050500:PRES:tropopause:anl:  looks like Pa
  344:59712381:d=1990050500:HGT:tropopause:anl:   looks like m
  345:59999169:d=1990050500:TMP:tropopause:anl:   in K, looks like T (not theta)
  346.1:60085837:d=1990050500:UGRD:tropopause:anl: looks like m/s
  346.2:60085837:d=1990050500:VGRD:tropopause:anl: looks like m/s
  347:60308670:d=1990050500:VWSH:tropopause:anl:  vertical wind shear
  '''
  fpath = '/home/nickszap/Downloads/pgbhnl.gdas.1990050500.grb2'
  #os.system("wgrib2 "+fpath)
  #os.system('wgrib2 '+fpath+' -match "(TMP:tropopause)"')
  #wgrib2 pgbhnl.gdas.1990050500.grb2 -match "(TMP:tropopause)" -grib pgbhnl.gdas.1990050500.small.grb2
  #matchString = '"(TMP:tropopause)"'
  matchString = '"(tropopause)"'
  outName = fpath+'.trop.grb2'
  cmd = 'wgrib2 '+fpath+' -match '+matchString+' -grib '+outName
  os.system(cmd)
  
  cmd = 'ncl_convert2nc '+outName
  os.system(cmd) #can plot this with Paraview using generic netcdf reader
