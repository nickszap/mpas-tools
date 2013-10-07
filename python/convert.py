
import output_data
import vars

def derivedSfcs(ncfname, vtkfname):
  #write some derived surfaces to file
  
  data = output_data.open_netcdf_data(ncfname)
  
  #header info
  fvtk = output_data.write_vtk_header_polydata(vtkfname, ncfname)
  nNodes =  output_data.write_vtk_xyzNodes(fvtk, data)
  nCells = output_data.write_vtk_polygons(fvtk, data)
  nLevels =  len(data.dimensions['nVertLevels'])
  
  #cell data
  fvtk.write('\nCELL_DATA '+str(nCells)+'\n')
  
  #geo for reference
  output_data.write_vtk_staticGeoFields(f,data,nCells)
  
  time = 0
  #500 mb
  output_data.write_vtk_pressureHeights(fvtk, data, nCells, time, vLevels, 50000.)
  
  #theta on dynamic tropopause
  pv = np.empty((nCells,nLevels), dtype=float)
  for hcell in range(nCells):
    for l in range(nLevels):
      pv[hcell,l] = vars.calc_ertelPV(data, 'theta', time, hcell, l, nLevels)
  #
  
  pvuVal = 2.
  thetaVal = np.empty(nCells)
  for hcell in range(nCells):
    (l,dl) = output_data.calcIndexOfValue(pvuVal,pv[hcell,:], nLevels)
    thetaVal[hcell] = output_data.calcValueOfIndex(l,dl,data.variables['theta'][time,hcell,:])
  output_data.write_levelData_float('theta_pv', fvtk, thetaVal, nCells)
  
  #slp
  slp = np.empty(nCells)
  for hcell in range(nCells):
    slp[hcell] = vars.calc_slp(data, hcell, nLevels, time)
  output_data.write_levelData_float('slp', fvtk, slp, nCells)
  
  #close da files
  fvtk.close()
  data.close()

def make_vtkName(base, tInd, lInd):
  #visit and paraview have an expected format for file names for animation.
  #These are listed in different formats but are generally name_N.vtk for index N.
  #Inputs are base of string, time and level indices
  
  s = base+'_l'+str(lInd)+'_t'+str(tInd)+'.vtk'
  return s

def example(ncNameFile):
  #Input file is a list of each output file on its own line.
  #Paraview has a hard time animating the 0 time step .nc files. Maybe it's an issue with the MPAS reader?
  #We can get around it by creating vtk files of the fields of interest named fnameN.vtk where N is an integer that indicates time.
  
  #store the names of the files we want in order in a file
  #ncNameFile = 'ncNames.txt'
  fp = open(ncNameFile,'r')

  #output the files to
  outNameBase = 'june'
  
  #
  i = 0
  for line in fp:
    ncfname = line.rstrip('\n') #do we have to strip off any characters like \n?
    
    if (i==0):
      #reference field for orientation
      vtkfname = outNameBase+'_ref'+'.vtk'
      data = output_data.open_netcdf_data(ncfname)

      #header and mesh info
      fvtk = output_data.write_vtk_header_polydata(vtkfname, ncfname)
      nNodes =  output_data.write_vtk_xyzNodes(fvtk, data)
      nCells = output_data.write_vtk_polygons(fvtk, data)
      fvtk.write('\nCELL_DATA '+str(nCells)+'\n')
      output_data.write_vtk_staticGeoFields(f,data,nCells)
      fvtk.close()
      data.close()
    #
    
    vtkfname = outNameBase+str(i)+'.vtk'

    data = output_data.open_netcdf_data(ncfname)
    
    #header and mesh info
    fvtk = output_data.write_vtk_header_polydata(vtkfname, ncfname)
    nNodes =  output_data.write_vtk_xyzNodes(fvtk, data)
    nCells = output_data.write_vtk_polygons(fvtk, data)
    
    #write cell data
    fvtk.write('\nCELL_DATA '+str(nCells)+'\n')
    time = 0
    vLevel = 18
    output_data.write_vtk_cellCenterVelocity(fvtk, data, time, vLevel, nCells)
    output_data.write_vtk_var_timeLevelCells(fvtk, 'pv_cell', data, vLevel, time, nCells)
    
    i = i+1
    #close files
    data.close()
    fvtk.close()

  fp.close()

if __name__=='__main__':
  example('ncNames.txt')
  
#
