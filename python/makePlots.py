#module for examples of generating fields from mpas output.nc fields

#python libraries
import os
import datetime as dt

#my libraries
import sys; sys.path.append("/home/nickszap/Dropbox/pythonScripts");
import output_data
import cfsr
import vars

def example():
  
  #file properties
  ncfname = '/arctic1/nick/cases/cfsr/output.2006-08-07_12.00.00.nc'
  #ncfname = '/home/nickszap/research/mpas/output.2010-10-23_00:00:00.nc' #input file
  vtkfname = 'plotTest.vtk' #output file
  
  data = output_data.open_netcdf_data(ncfname)
  
  #open the output vtk file and write header. -------------------
  fvtk = output_data.write_vtk_header_polydata(vtkfname, ncfname)
  
  #write nodes and cells
  nNodes =  output_data.write_vtk_xyzNodes(fvtk, data)
  nCells = output_data.write_vtk_polygons(fvtk, data)
  nLevels =  len(data.dimensions['nVertLevels'])
  
  #write some cell data --------------------
  fvtk.write('\nCELL_DATA '+str(nCells)+'\n')
  
  output_data.write_vtk_staticGeoFields(fvtk,data,nCells)
  
  #time dependent stuff goes in different files
  time = 0
  output_data.write_vtk_pressureHeights(fvtk, data, nCells, time, nLevels, 50000) #500mb = 50,000Pa
  
  #write some node data
  
  #close the .nc and vtk files
  data.close()
  fvtk.close()

def example_vars():
  #for every 6 hours of the IC CFSR data,
  #write a vtk file with slp, 500mb heights, and theta on 2 pv surface
  
  t0 = dt.datetime(2006,6,1,0)
  #tf = dt.datetime(2006,9,29,18)
  tf = dt.datetime(2006,6,1,13)
  h6 = dt.timedelta(hours=6)
  
  cfsrPath = '/arctic1/nick/cases/cfsr/'
  vtkBase = 'cfsrIC.'
  
  i=1; t=t0;
  while (t<=tf): #since increment t after check, don't do
    #open the .nc data file for this datetime
    tString = cfsr.form_cfsrTimeString(t)
    ncName = 'vert_sfc.'+tString+'.nc' #initial condition netcdf file
    ncName = cfsrPath+ncName    
    data = output_data.open_netcdf_data(ncName)
  
    #open the output vtk file and write header. -------------------
    vtkfname = vtkBase+str(i-1)+'.vtk'
    fvtk = output_data.write_vtk_header_polydata(vtkfname, ncName)

    #write nodes and cells
    nNodes =  output_data.write_vtk_xyzNodes(fvtk, data)
    nCells = output_data.write_vtk_polygons(fvtk, data)
    nLevels =  len(data.dimensions['nVertLevels'])

    #write some cell data --------------------
    fvtk.write('\nCELL_DATA '+str(nCells)+'\n')
    
    timeInd = 0
    
    #
    data.close()
    fvtk.close()
    
    #increment day
    t = t0+i*h6; i = i+1;

if __name__=='__main__':
  example()


