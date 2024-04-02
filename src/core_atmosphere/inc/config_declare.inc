      character (len=StrKIND), pointer :: config_time_integration
      integer, pointer :: config_time_integration_order
      real (kind=RKIND), pointer :: config_dt
      character (len=StrKIND), pointer :: config_calendar_type
      character (len=StrKIND), pointer :: config_start_time
      character (len=StrKIND), pointer :: config_stop_time
      character (len=StrKIND), pointer :: config_run_duration
      logical, pointer :: config_split_dynamics_transport
      integer, pointer :: config_number_of_sub_steps
      integer, pointer :: config_dynamics_split_steps
      real (kind=RKIND), pointer :: config_h_mom_eddy_visc2
      real (kind=RKIND), pointer :: config_h_mom_eddy_visc4
      real (kind=RKIND), pointer :: config_v_mom_eddy_visc2
      real (kind=RKIND), pointer :: config_h_theta_eddy_visc2
      real (kind=RKIND), pointer :: config_h_theta_eddy_visc4
      real (kind=RKIND), pointer :: config_v_theta_eddy_visc2
      character (len=StrKIND), pointer :: config_horiz_mixing
      real (kind=RKIND), pointer :: config_len_disp
      real (kind=RKIND), pointer :: config_visc4_2dsmag
      real (kind=RKIND), pointer :: config_del4u_div_factor
      integer, pointer :: config_w_adv_order
      integer, pointer :: config_theta_adv_order
      integer, pointer :: config_scalar_adv_order
      integer, pointer :: config_u_vadv_order
      integer, pointer :: config_w_vadv_order
      integer, pointer :: config_theta_vadv_order
      integer, pointer :: config_scalar_vadv_order
      logical, pointer :: config_scalar_advection
      logical, pointer :: config_positive_definite
      logical, pointer :: config_monotonic
      real (kind=RKIND), pointer :: config_coef_3rd_order
      real (kind=RKIND), pointer :: config_smagorinsky_coef
      logical, pointer :: config_mix_full
      real (kind=RKIND), pointer :: config_epssm
      real (kind=RKIND), pointer :: config_smdiv
      real (kind=RKIND), pointer :: config_apvm_upwinding
      logical, pointer :: config_h_ScaleWithMesh
      integer, pointer :: config_num_halos
      real (kind=RKIND), pointer :: config_relax_zone_divdamp_coef

      real (kind=RKIND), pointer :: config_zd
      real (kind=RKIND), pointer :: config_xnutr
      real (kind=RKIND), pointer :: config_mpas_cam_coef
      integer, pointer :: config_number_cam_damping_levels
      logical, pointer :: config_rayleigh_damp_u
      real (kind=RKIND), pointer :: config_rayleigh_damp_u_timescale_days
      integer, pointer :: config_number_rayleigh_damp_u_levels

      logical, pointer :: config_apply_lbcs

      character (len=StrKIND), pointer :: config_restart_timestamp_name
      integer, pointer :: config_pio_num_iotasks
      integer, pointer :: config_pio_stride

      character (len=StrKIND), pointer :: config_block_decomp_file_prefix
      integer, pointer :: config_number_of_blocks
      logical, pointer :: config_explicit_proc_decomp
      character (len=StrKIND), pointer :: config_proc_decomp_file_prefix

      logical, pointer :: config_do_restart
      logical, pointer :: config_do_DAcycling

      logical, pointer :: config_print_global_minmax_vel
      logical, pointer :: config_print_detailed_minmax_vel
      logical, pointer :: config_print_global_minmax_sca

      character (len=StrKIND), pointer :: config_IAU_option
      real (kind=RKIND), pointer :: config_IAU_window_length_s

      logical, pointer :: config_jedi_da

      character (len=StrKIND), pointer :: config_halo_exch_method

      character (len=StrKIND), pointer :: input_soil_data
      integer, pointer :: input_soil_temperature_lag
      integer, pointer :: num_soil_layers
      integer, pointer :: months
      integer, pointer :: noznlev
      integer, pointer :: naerlev
      integer, pointer :: camdim1
      logical, pointer :: config_frac_seaice
      logical, pointer :: config_sfc_albedo
      logical, pointer :: config_sfc_snowalbedo
      logical, pointer :: config_sst_update
      logical, pointer :: config_sstdiurn_update
      logical, pointer :: config_deepsoiltemp_update
      logical, pointer :: config_o3climatology
      logical, pointer :: config_microp_re
      logical, pointer :: config_ysu_pblmix
      integer, pointer :: config_n_microp
      character (len=StrKIND), pointer :: config_radtlw_interval
      character (len=StrKIND), pointer :: config_radtsw_interval
      character (len=StrKIND), pointer :: config_conv_interval
      character (len=StrKIND), pointer :: config_pbl_interval
      character (len=StrKIND), pointer :: config_camrad_abs_update
      character (len=StrKIND), pointer :: config_greeness_update
      character (len=StrKIND), pointer :: config_bucket_update
      character (len=StrKIND), pointer :: config_physics_suite
      character (len=StrKIND), pointer :: config_microp_scheme
      character (len=StrKIND), pointer :: config_convection_scheme
      character (len=StrKIND), pointer :: config_lsm_scheme
      character (len=StrKIND), pointer :: config_pbl_scheme
      character (len=StrKIND), pointer :: config_gwdo_scheme
      character (len=StrKIND), pointer :: config_radt_cld_scheme
      character (len=StrKIND), pointer :: config_radt_lw_scheme
      character (len=StrKIND), pointer :: config_radt_sw_scheme
      character (len=StrKIND), pointer :: config_sfclayer_scheme
      integer, pointer :: config_gfconv_closure_deep
      integer, pointer :: config_gfconv_closure_shallow
      real (kind=RKIND), pointer :: config_bucket_radt
      real (kind=RKIND), pointer :: config_bucket_rainc
      real (kind=RKIND), pointer :: config_bucket_rainnc
      logical, pointer :: config_oml1d
      real (kind=RKIND), pointer :: config_oml_hml0
      real (kind=RKIND), pointer :: config_oml_gamma
      real (kind=RKIND), pointer :: config_oml_relaxation_time

      character (len=StrKIND), pointer :: config_sounding_interval

