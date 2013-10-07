# -*- coding: utf-8 -*-

# imports
import netCDF4
import numpy as np
import output_data

def driver_ncGeo2vtk(ncfname):
  '''
  given .nc (MPAS NETCDF) file, output the surface geography data into classic vtk format.
  The MPAS NetCDF reader in Paraview loads only the "critical" netcdf variables, ie with time and vertLevels.
  Rather than edit that reader, we can create a file of the vars we care about on the scvt mesh rather than the dual.
  I think they have to use the dual for the volume mesh since the elements need to be of supported VTK type (eg prism, hex,...)
  '''
  
  #file properties
  #ncfname = '/home/nickszap/research/mpas/output.2010-10-23_00:00:00.nc' #input file
  vtkfname = 'geoTest.vtk' #output file
  
  data = output_data.open_netcdf_data(ncfname)
  
  #open the output vtk file and write header.
  fvtk = output_data.write_vtk_header_polydata(vtkfname, ncfname)
  
  #write nodes and cells
  nNodes =  output_data.write_vtk_xyzNodes(fvtk, data)
  nCells = output_data.write_vtk_polygons(fvtk, data)
  vLevels =  len(data.dimensions['nVertLevels'])
  nTimes = len(data.dimensions['Time'])
  
  #write some geographic cell data
  fvtk.write('\nCELL_DATA '+str(nCells)+'\n')
  write_vtk_staticGeoFields(fvtk,data,nCells)
  write_vtk_timeGeoFields(fvtk,data,nCells, nTimes)
  
  #write node data (none for geog)
  
  #close the .nc and vtk files
  data.close()
  fvtk.close()
#end driver

def write_vtk_staticGeoFields(f,data,nCells):
  #write some surface fields of the cells in our scvt mesh that don't change over time
  
  #could use a dict if want titles of vars different from keys
  vars = ['ter','landmask','snoalb']
  for var in vars:
    f.write('\nSCALARS '+var+' float 1\nLOOKUP_TABLE default\n')
    vals = data.variables[var][:]
    for cell in range(nCells):
      string  = str(float(vals[cell]))+'\n'
      f.write(string)
#end

def write_vtk_timeGeoFields(f,data,nCells, nTimes):
  #write some surface fields of the cells in our scvt mesh that do change over time
  vars = ['vegfra','xice','sfc_albbck', 'skintemp', 'sst', 'tmn'] #use seaice if want xice binary
  for t in range(nTimes):
    for var in vars:
      f.write('\nSCALARS '+var+str(t)+' float 1\nLOOKUP_TABLE default\n')
      vals = data.variables[var][t,:]
      for cell in range(nCells):
        string  = str(float(vals[cell]))+'\n'
        f.write(string)
#end
  
def write_vtk_cellLandType(f, data, nCells):
  #string = '\nCELL_DATA '+str(nCells)+'\n'
  #f.write(string)
  f.write('\nSCALARS waterLand int 1\nLOOKUP_TABLE default\n')
  wl = data.variables['landmask'][:]
  for i in range(nCells):
    string  = str(wl[i])+'\n'
    f.write(string)



