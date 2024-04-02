      integer, pointer :: config_init_case
      character (len=StrKIND), pointer :: config_calendar_type
      character (len=StrKIND), pointer :: config_start_time
      character (len=StrKIND), pointer :: config_stop_time
      integer, pointer :: config_theta_adv_order
      real (kind=RKIND), pointer :: config_coef_3rd_order
      integer, pointer :: config_num_halos
      character (len=StrKIND), pointer :: config_interface_projection

      integer, pointer :: config_nvertlevels
      integer, pointer :: config_nsoillevels
      integer, pointer :: config_nfglevels
      integer, pointer :: config_nfgsoillevels
      integer, pointer :: config_months

      character (len=StrKIND), pointer :: config_geog_data_path
      character (len=StrKIND), pointer :: config_met_prefix
      character (len=StrKIND), pointer :: config_sfc_prefix
      integer, pointer :: config_fg_interval
      character (len=StrKIND), pointer :: config_landuse_data
      character (len=StrKIND), pointer :: config_topo_data
      character (len=StrKIND), pointer :: config_vegfrac_data
      character (len=StrKIND), pointer :: config_albedo_data
      character (len=StrKIND), pointer :: config_maxsnowalbedo_data
      integer, pointer :: config_supersample_factor
      logical, pointer :: config_use_spechumd

      real (kind=RKIND), pointer :: config_ztop
      integer, pointer :: config_nsmterrain
      logical, pointer :: config_smooth_surfaces
      real (kind=RKIND), pointer :: config_dzmin
      integer, pointer :: config_nsm
      logical, pointer :: config_tc_vertical_grid
      character (len=StrKIND), pointer :: config_specified_zeta_levels
      logical, pointer :: config_blend_bdy_terrain

      character (len=StrKIND), pointer :: config_extrap_airtemp

      logical, pointer :: config_static_interp
      logical, pointer :: config_native_gwd_static
      real (kind=RKIND), pointer :: config_gwd_cell_scaling
      logical, pointer :: config_vertical_grid
      logical, pointer :: config_met_interp
      logical, pointer :: config_input_sst
      logical, pointer :: config_frac_seaice

      real (kind=RKIND), pointer :: config_tsk_seaice_threshold

      integer, pointer :: config_pio_num_iotasks
      integer, pointer :: config_pio_stride

      character (len=StrKIND), pointer :: config_block_decomp_file_prefix
      integer, pointer :: config_number_of_blocks
      logical, pointer :: config_explicit_proc_decomp
      character (len=StrKIND), pointer :: config_proc_decomp_file_prefix

