#code namelist sections as functions so that we can loop through options later.
#for init, the preproc_stages define the fields that get worked on.
#for later, s = '''testing %s''' %("\'?\'") gives s= testing '?' (as a string). %("'?'") also works i think
#Apparently, str.format() is the now preferred 

#format of MPAS time is '2006-07-06_06:00:00'. For example, do:
#time.strftime('%Y-%m-%d_%H:%M:%S', datetime.datetime.timetuple(datetime.datetime(2007,12,28,0)))

import arctic

def addQuotes2String(text):
  #many namelist options need to be written as 'text' (ie, within single quotes)
  text = "'%s'" %(text)
  return text
  
def writeNamelist_geo(f):
  s = '''&nhyd_model
   config_test_case = 7
/
'''
  s = s+formString_dimensionsGeo()
  f.write(s) #model and dims
  s = formString_data()
  f.write(s)
  #ts = "'.true.'"; fs = "'.false.'";
  ts = '.true.'; fs = '.false.'; #they don't have quotes!
  s = formString_preproc_vars(ts,fs,fs,fs,fs)
  f.write(s)
  s = formString_initIO_vars("'grid.nc'", "'static.nc'")
  f.write(s)
#end

def writeNamelist_vertSfc(f, tStart,nNodes, nProcsPerNode):
  #there's no harm in initializing both together as laura does.
  #since it's fast enough (~1 minute), easy-peezy to just run together.
  #Note that WPS intermediate file needs to have pressure and surface data concatenated!
  
  #tStart = "'2006-07-01_00:00:00'"
  s = formString_initModel_vars(7, tStart, tStart)
  f.write(s)
  s = formString_dimensions()
  f.write(s)
  s = formString_data()
  f.write(s)
  s = formString_vertical()
  f.write(s)
  ts = '.true.'; fs = '.false.'; #they don't have quotes!
  s = formString_preproc_vars(fs,ts,ts,ts,ts)
  f.write(s)
  s = formString_initIO_vars("'static.163842.nc'", "'vert_sfc.nc'",nNodes,nProcsPerNode)
  f.write(s)
  s = formString_decomposition(-1)
  f.write(s)
#end

#def writeNamelist_vertSfc_cfsrCase(f, tStart,nNodes, nProcsPerNode):
def writeNamelist_vertSfc_cfsrCase(f, tStart):
  #there's no harm in initializing both together as laura does.
  #since it's fast enough (~1 minute), easy-peezy to just run together.
  #Note that WPS intermediate file needs to have pressure and surface data concatenated!
  
  #tStart = "'2006-07-01_00:00:00'"
  s = formString_initModel_vars(7, tStart, tStart)
  f.write(s)
  s = formString_dimensions()
  f.write(s)
  s = formString_data()
  f.write(s)
  s = formString_vertical()
  f.write(s)
  ts = '.true.'; fs = '.false.'; #they don't have quotes!
  s = formString_preproc_vars(fs,ts,ts,ts,ts)
  f.write(s)
  #s = formString_initIO_vars("'static.163842.nc'", "'vert_sfc."+tStart[0:13]+".nc'",nNodes,nProcsPerNode)
  s = formString_initIO_vars("'static.163842.nc'", "'vert_sfc."+tStart[0:13]+".nc'",0,1)
  f.write(s)
  s = formString_decomposition(163842)
  f.write(s)
#end

def writeNamelist_sfcUpdate(f,tStart,tEnd):
  
  s = formString_initModel_vars(8, tStart, tEnd)
  f.write(s)
#end

def writeNamelist_integrate(f, tStart, tEnd, nNodes, nProcsPerNode):
  #Namelist for integration of model: model, damping, dims, io, decomp, restart, physics
  
  s = formString_model(tStart, tEnd)
  f.write(s)
  s = formString_damping()
  f.write(s)
  s = formString_runDimensions()
  f.write(s)
  s = formString_runIO(nNodes, nProcsPerNode)
  f.write(s)
  s = formString_decomposition()
  f.write(s)
  s = formString_restart()
  f.write(s)
  s = formString_physics()
  f.write(s)
#end

#Initialization options ------------------------------------
#placing text on the same line as the ''' gives no newline
def formString_initModel():
  s = '''&nhyd_model
   config_test_case = 7
   config_theta_adv_order = 3
   config_start_time      = '2006-07-01_00:00:00'
   config_stop_time       = '2006-07-01_00:00:00'
/
'''
  return s

def formString_initModel_vars(case, tStart, tStop):
  s = '''&nhyd_model
   config_test_case = %d
   config_theta_adv_order = 3
   config_start_time      = %s
   config_stop_time       = %s
/
''' %(case, tStart, tStop)
  return s

def formString_dimensions():
  s = '''
&dimensions
   config_nvertlevels     = 41
   config_nsoillevels     = 4
   config_nfglevels       = 38
   config_nfgsoillevels   = 4
/
'''
  return s

def formString_dimensionsGeo():
  s = '''
&dimensions
   config_nvertlevels     = 3
   config_nsoillevels     = 4
   config_nfglevels       = 3
   config_nfgsoillevels   = 4
/
'''
  return s

def formString_data():
  s = '''
&data_sources
   config_geog_data_path  = '/raid1/comm/geog/'
   config_met_prefix      = 'CFSR'
   config_sfc_prefix      = 'CFSR'
   config_fg_interval     = 21600
/
'''
  return s

def formString_vertical():
  s = '''
&vertical_grid
   config_ztop            = 30000.0
   config_nsmterrain      = 2
   config_smooth_surfaces = .true.
/
'''
  return s

def formString_preproc():
  s = '''
&preproc_stages 
   config_static_interp   = .false.
   config_vertical_grid   = .true.
   config_met_interp      = .true.
   config_input_sst       = .true.
   config_frac_seaice     = .true.
/
'''
  return s

def formString_preproc_vars(static, vertical, met, sst, seaIce):
  #call as something like ...vars('.true.','.false.',...)
  s = '''
&preproc_stages
   config_static_interp   = %s
   config_vertical_grid   = %s
   config_met_interp      = %s
   config_input_sst       = %s
   config_frac_seaice     = %s
/
''' %(static, vertical, met, sst, seaIce)
  return s

def formString_initIO():
  #we think num_iotasks*pio_stride needs to equal # processes
  s = '''
&io
   config_input_name         = 'x1.40962.static.nc'
   config_output_name        = 'x1.40962.initFLXF06.2006-07-01_00.nc'
   config_pio_num_iotasks    = 10
   config_pio_stride         = 12
/
'''
  return s

def formString_initIO_vars(inputName, outputName, nNodes, nProcPerNode):
  #call as something like ...vars('x1.40962.static.nc','x1.40962.initFLXF06.2006-07-01_00.nc',10,12)
  
  #From cavallo's email after meeting at ncar on 2/14/2013,
  #config_pio_num_iotasks must be set to the number of machines that will be used 
  #(up to 10 on Arctic, or 9 if you don't want to use node 6), 
  #and config_pio_stride is the number of nodes per machines (up to 12 on Arctic)
  
  #there are issues with num_iotasks=1 so avoid that.
  #we'll assume that nProcsPerNode is even
  if (nNodes==1):
    nNodes = 2;
    nProcPerNode = nProcPerNode/2;
  
  s = '''
&io
   config_input_name         = %s
   config_output_name        = %s
   config_pio_num_iotasks    = %s
   config_pio_stride         = %s
/
''' %(inputName, outputName,nNodes, nProcPerNode)
  return s

#nhyd.exe options
def formString_decomposition(nCells):
  #graphFile = 'x1.{0}.graph.info.part.'.format(nCells)
  s = '''
&decomposition
   config_number_of_blocks = 0
   config_block_decomp_file_prefix = 'x1.{0}.graph.info.part.' 
   config_explicit_proc_decomp = .false.
   config_proc_decomp_file_prefix = 'x1.{0}.graph.info.part.'
/
'''.format(nCells)
  return s

def formString_restart():

  s = '''
&restart
   config_restart_interval   = '01_00:00:00'
   config_do_restart         = .false.
/
'''
  return s

def formString_model(tStart, tEnd):
  s = '''&nhyd_model
   config_time_integration = 'SRK3'
   config_dt = 450.0
   config_start_time = %s
   config_stop_time  = %s
   config_sfc_update_interval = 'none'
   config_number_of_sub_steps = 6
   config_h_mom_eddy_visc2    = 0.0e+04
   config_h_mom_eddy_visc4    = 0.0e+13
   config_v_mom_eddy_visc2    = 0.0e+04
   config_h_theta_eddy_visc2  = 0.0e+04
   config_h_theta_eddy_visc4  = 0.0e+13
   config_v_theta_eddy_visc2  = 0.0e+04
   config_horiz_mixing        = '2d_smagorinsky'
   config_len_disp            = 120000
   config_u_vadv_order        = 3
   config_w_vadv_order        = 3
   config_theta_vadv_order    = 3
   config_scalar_vadv_order   = 3
   config_theta_adv_order     = 3
   config_scalar_adv_order    = 3
   config_scalar_advection    = .true.
   config_positive_definite   = .false.
   config_monotonic           = .true.
   config_coef_3rd_order      = 0.25
   config_epssm               = 0.1
   config_smdiv               = 0.1
   config_h_ScaleWithMesh     = .false.
   config_newpx               = .false.
/
''' %(tStart, tEnd)
  return s

def formString_runIO(nNodes, nProcsPerNode):
  #there are issues with num_iotasks=1 so avoid that.
  #we'll assume that nProcsPerNode is even
  if (nNodes==1):
    nNodes = 2;
    nProcsPerNode = nProcsPerNode/2;
    
  s = '''
&io
   config_input_name         = 'vert_sfc.nc'
   config_output_name        = 'output.nc'
   config_restart_name       = 'restart.nc'
   config_sfc_update_name    = 'sfc_update.nc'
   config_output_interval    = '00_06:00:00'
   config_frames_per_outfile = 25
   config_pio_num_iotasks    = {0}
   config_pio_stride         = {1}
/
'''.format(nNodes,nProcsPerNode)
  return s

def formString_runDimensions():
  s = '''
&dimensions
   config_nvertlevels     = 41
/
'''
  return s

def formString_damping():
  s = '''
&damping
   config_zd = 26000.0
   config_xnutr = 0.2
/
'''
  return s

def formString_physics():
  s = '''
&physics
 config_frac_seaice          = .false.
 config_sfc_albedo           = .true.
 config_sfc_snowalbedo       = .true.
 config_sst_update           = .false.
 config_sstdiurn_update      = .false.
 config_deepsoiltemp_update  = .false.
 config_bucket_update        = '06:00:00'
 config_bucket_rainc         = 100.0
 config_bucket_rainnc        = 100.0
 config_bucket_radt          = 1.0e9
 config_radtlw_interval      = '00:30:00'
 config_radtsw_interval      = '00:30:00'
 config_conv_interval        = 'none'
 config_n_microp = 5
 config_microp_scheme        = 'wsm6'
 config_conv_shallow_scheme  = 'off'
 config_conv_deep_scheme     = 'kain_fritsch'
 config_eddy_scheme          = 'off'
 config_lsm_scheme           = 'noah'
 config_pbl_scheme           = 'ysu'
 config_radt_lw_scheme       = 'rrtmg_lw'
 config_radt_sw_scheme       = 'rrtmg_sw'
 config_sfclayer_scheme      = 'monin_obukhov'
/
'''
  return s

if __name__ == '__main__':
  print("Nope. Import and call functions\n")
  
  
