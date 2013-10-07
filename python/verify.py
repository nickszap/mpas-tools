# -*- coding: utf-8 -*-

#imports
import netCDF4
import numpy as np
import output_data

def example():
  #file properties
  ana_fname = '/home/nickszap/research/mpas/output.2010-10-23_00:00:00.nc' #analysis file
  sim_fname = '/home/nickszap/research/mpas/output.2010-10-23_00:00:00.nc' #simulation file
  
  vtkfname = 'fieldDiff.vtk' #output file
  
  #read in reference and simulation files
  anaData = output_data.open_netcdf_data(ana_fname)
  simData = output_data.open_netcdf_data(sim_fname)
  
  #assuming same mesh between ana and sim
  nTimes = len(data.dimensions['Time'])
  nCells = len(data.dimensions['nCells']);
  nVert = len(data.dimensions['nVertLevels'])
  
  #I don't know how to write a 3 or 4d field into the .nc format so we're going to have
  #to hold off on that.
  
  #take differences of individual fields (3d, height levels, pressure levels,...)
  #and generate output (surfaces-horiz and vertical, field statistics,...).
  tSim = 0 #index into time steps of simulation corresponding to analysis
  tAna = 0
  fieldKeys = ['pv_cell','theta','rho']
  for k in fieldKeys:
    #index is [time,cell,vLevel]
    varAna = data.variables[k][tAna,:,:]
    varSim = data.variables[k][tSim,:,:]
    diff = varSim-varAna
    
    avg = sum(diff)/(nCells*nVert)
    str = 'Variable %s has sim-ana mean= %s for simInd %s\n' %(k, avg, tSim)
    print str
    
    rms = sum(diff*diff)/(nCells*nVert) #square and mean
    rms = np.sqrt(rms) #root
    str = 'Variable %s has sim-ana RMS= %s for simInd %s\n' %(k, rms, tSim)
    print str
  
  #take differences of derived fields
  #pressure: 'pressure_p' + 'pressure_base'
  
  #velocity: vector addition of 'uReconstructX' + 'uReconstructY' + 'uReconstructZ'
  
def calc_rms(var, nCells, nVert):
  #calc the rms of a field
  diff=var
  rms = np.sqrt(sum(diff*diff)/(nCells*nVert)) #square, mean, root
  return rms
