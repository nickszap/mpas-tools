import os

def main():
  dates = ['2006-08-01_06:00:00','2006-07-20_00:00:00','2006-07-17_18:00:00', '2006-07-12_18:00:00','2006-07-10_00:00:00']

  fname_ic = 'namelist.ic1'
  fname_integrate = 'namelist.int1'

  for d in dates:
    #run ic----------------------
    #write namelist    
    f = open(fname_ic,'w')
    s = namelist_ic(d)
    f.write(s)
    f.close()

    #link to proper name
    cmd = 'rm namelist.input'
    os.system(cmd)
    cmd = 'ln -s '+fname_ic+' namelist.input'
    os.system(cmd)

    #generate ic
    cmd = 'time mpiexec -n 1 ./init_atmosphere_model'
    os.system(cmd)

    #run integration-------------------------
    #write namelist
    f = open(fname_integrate,'w')
    s = namelist_integrate(d)
    f.write(s)
    f.close()

    #link to proper name
    cmd = 'rm namelist.input'
    os.system(cmd)
    cmd = 'ln -s '+fname_integrate+' namelist.input'
    os.system(cmd)

    #integrate
    cmd = 'time mpiexec -f ~/machines -n 108 ./atmosphere_model'
    os.system(cmd)

def namelist_ic(startTime):
  #startTime = '2006-07-16_06:00:00'
  s = '''&nhyd_model
   config_test_case = 7
   config_theta_adv_order = 3
   config_start_time      = {0}
   config_stop_time       = {0}
/

&dcmip
   config_dcmip_case          = '2-0-0'
   config_planet_scale        = 1.0
   config_rotation_rate_scale = 1.0
/

&dimensions
   config_nvertlevels     = 41
   config_nsoillevels     = 4
   config_nfglevels       = 38
   config_nfgsoillevels   = 4
/

&data_sources
   config_geog_data_path  = '/raid1/comm/geog/'
   config_met_prefix      = 'CFSR'
   config_sfc_prefix      = 'CFSR'
   config_fg_interval     = 21600
/

&vertical_grid
   config_ztop            = 30000.0
   config_nsmterrain      = 2
   config_smooth_surfaces = .true.
/

&preproc_stages 
   config_static_interp   = .false.
   config_vertical_grid   = .true.
   config_met_interp      = .true.
   config_input_sst       = .true.
/

&io
   config_input_name         = 'x1.163842.static.nc'
   config_output_name        = 'x1.163842.init.nc'
   config_sfc_update_name    = 'x1.163842.surface.nc'
   config_pio_num_iotasks    = 0
   config_pio_stride         = 1
/

&decomposition
   config_number_of_blocks = 0
   config_block_decomp_file_prefix = 'x1.163842.graph.info.part.' 
   config_explicit_proc_decomp = .false.
   config_proc_decomp_file_prefix = 'graph.info.part.'
/

&restart
/
'''.format(startTime)
  return s

def namelist_integrate(startTime):
  #startTime = '2006-07-16_06:00:00'
  s = '''&nhyd_model
   config_time_integration = 'SRK3'
   config_dt = 180.0
   config_start_time   = {0}
   config_run_duration = '14_00:00:00'
   config_number_of_sub_steps = 6
   config_h_mom_eddy_visc2    = 0.0
   config_h_mom_eddy_visc4    = 0.0
   config_v_mom_eddy_visc2    = 0.0
   config_h_theta_eddy_visc2  = 0.0
   config_h_theta_eddy_visc4  = 0.0
   config_v_theta_eddy_visc2  = 0.0
   config_horiz_mixing        = '2d_smagorinsky'
   config_len_disp            = 60000.0
   config_visc4_2dsmag        = 0.05
   config_u_vadv_order        = 3
   config_w_vadv_order        = 3
   config_theta_vadv_order    = 3
   config_scalar_vadv_order   = 3
   config_w_adv_order         = 3
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
   config_sfc_update_interval = '6:00:00'
/

&damping
   config_zd = 22000.0
   config_xnutr = 0.2
/

&io
   config_input_name         = 'x1.163842.init.nc'
   config_output_name        = 'x1.163842.output.nc'
   config_restart_name       = 'x1.163842.restart.nc'
   config_sfc_update_name    = 'x1.163842.surface.nc'
   config_output_interval    = '6:00:00'
   config_frames_per_outfile = 28
   config_pio_num_iotasks    = 9
   config_pio_stride         = 12
/

&decomposition
   config_number_of_blocks         = 0
   config_block_decomp_file_prefix = 'x1.163842.graph.info.part.'
/

&restart
   config_restart_interval = '3_12:00:00'
   config_do_restart = .false.
/

&physics
   config_frac_seaice         = .true.
   config_sfc_albedo          = .true.
   config_sfc_snowalbedo      = .true.
   config_sst_update          = .true.
   config_sstdiurn_update     = .false.
   config_deepsoiltemp_update = .false.
   config_bucket_update       = 'none'
   config_bucket_rainc        = 100.0
   config_bucket_rainnc       = 100.0
   config_bucket_radt         = 1.0e9
   config_radtlw_interval     = '00:30:00'
   config_radtsw_interval     = '00:30:00'
   config_conv_interval       = 'none'
   config_pbl_interval        = 'none'
   config_n_microp            = 3
   config_microp_scheme       = 'wsm6'
   config_conv_shallow_scheme = 'off'
   config_conv_deep_scheme    = 'kain_fritsch'
   config_eddy_scheme         = 'off'
   config_lsm_scheme          = 'noah'
   config_pbl_scheme          = 'ysu'
   config_gwdo_scheme         = 'off'
   config_radt_cld_scheme     = 'off'
   config_radt_lw_scheme      = 'rrtmg_lw'
   config_radt_sw_scheme      = 'rrtmg_sw'
   config_sfclayer_scheme     = 'monin_obukhov'
/
'''.format(startTime)
  return s


if __name__=='__main__':
  main()
