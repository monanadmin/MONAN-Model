   subroutine init_atm_setup_nmlrec_nhyd_model(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      integer :: config_init_case = 7
      character (len=StrKIND) :: config_calendar_type = 'gregorian'
      character (len=StrKIND) :: config_start_time = '2010-10-23_00:00:00'
      character (len=StrKIND) :: config_stop_time = '2010-10-23_00:00:00'
      integer :: config_theta_adv_order = 3
      real (kind=RKIND) :: config_coef_3rd_order = 0.250000
      integer :: config_num_halos = 2
      character (len=StrKIND) :: config_interface_projection = 'linear_interpolation'

      namelist /nhyd_model/ &
         config_init_case, &
         config_calendar_type, &
         config_start_time, &
         config_stop_time, &
         config_theta_adv_order, &
         config_coef_3rd_order, &
         config_num_halos, &
         config_interface_projection
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, nhyd_model, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_int(dminfo, config_init_case)
         call mpas_dmpar_bcast_char(dminfo, config_calendar_type)
         call mpas_dmpar_bcast_char(dminfo, config_start_time)
         call mpas_dmpar_bcast_char(dminfo, config_stop_time)
         call mpas_dmpar_bcast_int(dminfo, config_theta_adv_order)
         call mpas_dmpar_bcast_real(dminfo, config_coef_3rd_order)
         call mpas_dmpar_bcast_int(dminfo, config_num_halos)
         call mpas_dmpar_bcast_char(dminfo, config_interface_projection)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record nhyd_model')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_init_case = $i', intArgs=(/config_init_case/))
            call mpas_log_write('        config_calendar_type = '//mpas_log_escape_dollars(config_calendar_type))
            call mpas_log_write('        config_start_time = '//mpas_log_escape_dollars(config_start_time))
            call mpas_log_write('        config_stop_time = '//mpas_log_escape_dollars(config_stop_time))
            call mpas_log_write('        config_theta_adv_order = $i', intArgs=(/config_theta_adv_order/))
            call mpas_log_write('        config_coef_3rd_order = $r', realArgs=(/config_coef_3rd_order/))
            call mpas_log_write('        config_num_halos = $i', intArgs=(/config_num_halos/))
            call mpas_log_write('        config_interface_projection = '//mpas_log_escape_dollars(config_interface_projection))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record nhyd_model.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_init_case', config_init_case)
      call mpas_pool_add_config(configPool, 'config_calendar_type', config_calendar_type)
      call mpas_pool_add_config(configPool, 'config_start_time', config_start_time)
      call mpas_pool_add_config(configPool, 'config_stop_time', config_stop_time)
      call mpas_pool_add_config(configPool, 'config_theta_adv_order', config_theta_adv_order)
      call mpas_pool_add_config(configPool, 'config_coef_3rd_order', config_coef_3rd_order)
      call mpas_pool_add_config(configPool, 'config_num_halos', config_num_halos)
      call mpas_pool_add_config(configPool, 'config_interface_projection', config_interface_projection)

   end subroutine init_atm_setup_nmlrec_nhyd_model


   subroutine init_atm_setup_nmlrec_dimensions(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      integer :: config_nvertlevels = 55
      integer :: config_nsoillevels = 4
      integer :: config_nfglevels = 38
      integer :: config_nfgsoillevels = 4
      integer :: config_months = 12

      namelist /dimensions/ &
         config_nvertlevels, &
         config_nsoillevels, &
         config_nfglevels, &
         config_nfgsoillevels, &
         config_months
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, dimensions, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_int(dminfo, config_nvertlevels)
         call mpas_dmpar_bcast_int(dminfo, config_nsoillevels)
         call mpas_dmpar_bcast_int(dminfo, config_nfglevels)
         call mpas_dmpar_bcast_int(dminfo, config_nfgsoillevels)
         call mpas_dmpar_bcast_int(dminfo, config_months)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record dimensions')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_nvertlevels = $i', intArgs=(/config_nvertlevels/))
            call mpas_log_write('        config_nsoillevels = $i', intArgs=(/config_nsoillevels/))
            call mpas_log_write('        config_nfglevels = $i', intArgs=(/config_nfglevels/))
            call mpas_log_write('        config_nfgsoillevels = $i', intArgs=(/config_nfgsoillevels/))
            call mpas_log_write('        config_months = $i', intArgs=(/config_months/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record dimensions.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_nvertlevels', config_nvertlevels)
      call mpas_pool_add_config(configPool, 'config_nsoillevels', config_nsoillevels)
      call mpas_pool_add_config(configPool, 'config_nfglevels', config_nfglevels)
      call mpas_pool_add_config(configPool, 'config_nfgsoillevels', config_nfgsoillevels)
      call mpas_pool_add_config(configPool, 'config_months', config_months)

   end subroutine init_atm_setup_nmlrec_dimensions


   subroutine init_atm_setup_nmlrec_data_sources(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      character (len=StrKIND) :: config_geog_data_path = '/glade/work/wrfhelp/WPS_GEOG/'
      character (len=StrKIND) :: config_met_prefix = 'CFSR'
      character (len=StrKIND) :: config_sfc_prefix = 'SST'
      integer :: config_fg_interval = 86400
      character (len=StrKIND) :: config_landuse_data = 'MODIFIED_IGBP_MODIS_NOAH'
      character (len=StrKIND) :: config_topo_data = 'GMTED2010'
      character (len=StrKIND) :: config_vegfrac_data = 'MODIS'
      character (len=StrKIND) :: config_albedo_data = 'MODIS'
      character (len=StrKIND) :: config_maxsnowalbedo_data = 'MODIS'
      integer :: config_supersample_factor = 3
      logical :: config_use_spechumd = .false.

      namelist /data_sources/ &
         config_geog_data_path, &
         config_met_prefix, &
         config_sfc_prefix, &
         config_fg_interval, &
         config_landuse_data, &
         config_topo_data, &
         config_vegfrac_data, &
         config_albedo_data, &
         config_maxsnowalbedo_data, &
         config_supersample_factor, &
         config_use_spechumd
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, data_sources, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_char(dminfo, config_geog_data_path)
         call mpas_dmpar_bcast_char(dminfo, config_met_prefix)
         call mpas_dmpar_bcast_char(dminfo, config_sfc_prefix)
         call mpas_dmpar_bcast_int(dminfo, config_fg_interval)
         call mpas_dmpar_bcast_char(dminfo, config_landuse_data)
         call mpas_dmpar_bcast_char(dminfo, config_topo_data)
         call mpas_dmpar_bcast_char(dminfo, config_vegfrac_data)
         call mpas_dmpar_bcast_char(dminfo, config_albedo_data)
         call mpas_dmpar_bcast_char(dminfo, config_maxsnowalbedo_data)
         call mpas_dmpar_bcast_int(dminfo, config_supersample_factor)
         call mpas_dmpar_bcast_logical(dminfo, config_use_spechumd)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record data_sources')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_geog_data_path = '//mpas_log_escape_dollars(config_geog_data_path))
            call mpas_log_write('        config_met_prefix = '//mpas_log_escape_dollars(config_met_prefix))
            call mpas_log_write('        config_sfc_prefix = '//mpas_log_escape_dollars(config_sfc_prefix))
            call mpas_log_write('        config_fg_interval = $i', intArgs=(/config_fg_interval/))
            call mpas_log_write('        config_landuse_data = '//mpas_log_escape_dollars(config_landuse_data))
            call mpas_log_write('        config_topo_data = '//mpas_log_escape_dollars(config_topo_data))
            call mpas_log_write('        config_vegfrac_data = '//mpas_log_escape_dollars(config_vegfrac_data))
            call mpas_log_write('        config_albedo_data = '//mpas_log_escape_dollars(config_albedo_data))
            call mpas_log_write('        config_maxsnowalbedo_data = '//mpas_log_escape_dollars(config_maxsnowalbedo_data))
            call mpas_log_write('        config_supersample_factor = $i', intArgs=(/config_supersample_factor/))
            call mpas_log_write('        config_use_spechumd = $l', logicArgs=(/config_use_spechumd/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record data_sources.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_geog_data_path', config_geog_data_path)
      call mpas_pool_add_config(configPool, 'config_met_prefix', config_met_prefix)
      call mpas_pool_add_config(configPool, 'config_sfc_prefix', config_sfc_prefix)
      call mpas_pool_add_config(configPool, 'config_fg_interval', config_fg_interval)
      call mpas_pool_add_config(configPool, 'config_landuse_data', config_landuse_data)
      call mpas_pool_add_config(configPool, 'config_topo_data', config_topo_data)
      call mpas_pool_add_config(configPool, 'config_vegfrac_data', config_vegfrac_data)
      call mpas_pool_add_config(configPool, 'config_albedo_data', config_albedo_data)
      call mpas_pool_add_config(configPool, 'config_maxsnowalbedo_data', config_maxsnowalbedo_data)
      call mpas_pool_add_config(configPool, 'config_supersample_factor', config_supersample_factor)
      call mpas_pool_add_config(configPool, 'config_use_spechumd', config_use_spechumd)

   end subroutine init_atm_setup_nmlrec_data_sources


   subroutine init_atm_setup_nmlrec_vertical_grid(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      real (kind=RKIND) :: config_ztop = 30000.000000
      integer :: config_nsmterrain = 1
      logical :: config_smooth_surfaces = .true.
      real (kind=RKIND) :: config_dzmin = 0.300000
      integer :: config_nsm = 30
      logical :: config_tc_vertical_grid = .true.
      character (len=StrKIND) :: config_specified_zeta_levels = ''
      logical :: config_blend_bdy_terrain = .false.

      namelist /vertical_grid/ &
         config_ztop, &
         config_nsmterrain, &
         config_smooth_surfaces, &
         config_dzmin, &
         config_nsm, &
         config_tc_vertical_grid, &
         config_specified_zeta_levels, &
         config_blend_bdy_terrain
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, vertical_grid, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_real(dminfo, config_ztop)
         call mpas_dmpar_bcast_int(dminfo, config_nsmterrain)
         call mpas_dmpar_bcast_logical(dminfo, config_smooth_surfaces)
         call mpas_dmpar_bcast_real(dminfo, config_dzmin)
         call mpas_dmpar_bcast_int(dminfo, config_nsm)
         call mpas_dmpar_bcast_logical(dminfo, config_tc_vertical_grid)
         call mpas_dmpar_bcast_char(dminfo, config_specified_zeta_levels)
         call mpas_dmpar_bcast_logical(dminfo, config_blend_bdy_terrain)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record vertical_grid')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_ztop = $r', realArgs=(/config_ztop/))
            call mpas_log_write('        config_nsmterrain = $i', intArgs=(/config_nsmterrain/))
            call mpas_log_write('        config_smooth_surfaces = $l', logicArgs=(/config_smooth_surfaces/))
            call mpas_log_write('        config_dzmin = $r', realArgs=(/config_dzmin/))
            call mpas_log_write('        config_nsm = $i', intArgs=(/config_nsm/))
            call mpas_log_write('        config_tc_vertical_grid = $l', logicArgs=(/config_tc_vertical_grid/))
            call mpas_log_write('        config_specified_zeta_levels = '//mpas_log_escape_dollars(config_specified_zeta_levels))
            call mpas_log_write('        config_blend_bdy_terrain = $l', logicArgs=(/config_blend_bdy_terrain/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record vertical_grid.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_ztop', config_ztop)
      call mpas_pool_add_config(configPool, 'config_nsmterrain', config_nsmterrain)
      call mpas_pool_add_config(configPool, 'config_smooth_surfaces', config_smooth_surfaces)
      call mpas_pool_add_config(configPool, 'config_dzmin', config_dzmin)
      call mpas_pool_add_config(configPool, 'config_nsm', config_nsm)
      call mpas_pool_add_config(configPool, 'config_tc_vertical_grid', config_tc_vertical_grid)
      call mpas_pool_add_config(configPool, 'config_specified_zeta_levels', config_specified_zeta_levels)
      call mpas_pool_add_config(configPool, 'config_blend_bdy_terrain', config_blend_bdy_terrain)

   end subroutine init_atm_setup_nmlrec_vertical_grid


   subroutine init_atm_setup_nmlrec_interpolation_control(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      character (len=StrKIND) :: config_extrap_airtemp = 'lapse-rate'

      namelist /interpolation_control/ &
         config_extrap_airtemp
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, interpolation_control, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_char(dminfo, config_extrap_airtemp)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record interpolation_control')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_extrap_airtemp = '//mpas_log_escape_dollars(config_extrap_airtemp))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record interpolation_control.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_extrap_airtemp', config_extrap_airtemp)

   end subroutine init_atm_setup_nmlrec_interpolation_control


   subroutine init_atm_setup_nmlrec_preproc_stages(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      logical :: config_static_interp = .true.
      logical :: config_native_gwd_static = .true.
      real (kind=RKIND) :: config_gwd_cell_scaling = 1.000000
      logical :: config_vertical_grid = .true.
      logical :: config_met_interp = .true.
      logical :: config_input_sst = .false.
      logical :: config_frac_seaice = .true.

      namelist /preproc_stages/ &
         config_static_interp, &
         config_native_gwd_static, &
         config_gwd_cell_scaling, &
         config_vertical_grid, &
         config_met_interp, &
         config_input_sst, &
         config_frac_seaice
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, preproc_stages, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_logical(dminfo, config_static_interp)
         call mpas_dmpar_bcast_logical(dminfo, config_native_gwd_static)
         call mpas_dmpar_bcast_real(dminfo, config_gwd_cell_scaling)
         call mpas_dmpar_bcast_logical(dminfo, config_vertical_grid)
         call mpas_dmpar_bcast_logical(dminfo, config_met_interp)
         call mpas_dmpar_bcast_logical(dminfo, config_input_sst)
         call mpas_dmpar_bcast_logical(dminfo, config_frac_seaice)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record preproc_stages')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_static_interp = $l', logicArgs=(/config_static_interp/))
            call mpas_log_write('        config_native_gwd_static = $l', logicArgs=(/config_native_gwd_static/))
            call mpas_log_write('        config_gwd_cell_scaling = $r', realArgs=(/config_gwd_cell_scaling/))
            call mpas_log_write('        config_vertical_grid = $l', logicArgs=(/config_vertical_grid/))
            call mpas_log_write('        config_met_interp = $l', logicArgs=(/config_met_interp/))
            call mpas_log_write('        config_input_sst = $l', logicArgs=(/config_input_sst/))
            call mpas_log_write('        config_frac_seaice = $l', logicArgs=(/config_frac_seaice/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record preproc_stages.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_static_interp', config_static_interp)
      call mpas_pool_add_config(configPool, 'config_native_gwd_static', config_native_gwd_static)
      call mpas_pool_add_config(configPool, 'config_gwd_cell_scaling', config_gwd_cell_scaling)
      call mpas_pool_add_config(configPool, 'config_vertical_grid', config_vertical_grid)
      call mpas_pool_add_config(configPool, 'config_met_interp', config_met_interp)
      call mpas_pool_add_config(configPool, 'config_input_sst', config_input_sst)
      call mpas_pool_add_config(configPool, 'config_frac_seaice', config_frac_seaice)

   end subroutine init_atm_setup_nmlrec_preproc_stages


   subroutine init_atm_setup_nmlrec_physics(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      real (kind=RKIND) :: config_tsk_seaice_threshold = 100.000000

      namelist /physics/ &
         config_tsk_seaice_threshold
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, physics, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_real(dminfo, config_tsk_seaice_threshold)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record physics')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_tsk_seaice_threshold = $r', realArgs=(/config_tsk_seaice_threshold/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record physics.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_tsk_seaice_threshold', config_tsk_seaice_threshold)

   end subroutine init_atm_setup_nmlrec_physics


   subroutine init_atm_setup_nmlrec_io(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      integer :: config_pio_num_iotasks = 0
      integer :: config_pio_stride = 1

      namelist /io/ &
         config_pio_num_iotasks, &
         config_pio_stride
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, io, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_int(dminfo, config_pio_num_iotasks)
         call mpas_dmpar_bcast_int(dminfo, config_pio_stride)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record io')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_pio_num_iotasks = $i', intArgs=(/config_pio_num_iotasks/))
            call mpas_log_write('        config_pio_stride = $i', intArgs=(/config_pio_stride/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record io.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_pio_num_iotasks', config_pio_num_iotasks)
      call mpas_pool_add_config(configPool, 'config_pio_stride', config_pio_stride)

   end subroutine init_atm_setup_nmlrec_io


   subroutine init_atm_setup_nmlrec_decomposition(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      character (len=StrKIND) :: config_block_decomp_file_prefix = 'x1.40962.graph.info.part.'
      integer :: config_number_of_blocks = 0
      logical :: config_explicit_proc_decomp = .false.
      character (len=StrKIND) :: config_proc_decomp_file_prefix = 'graph.info.part.'

      namelist /decomposition/ &
         config_block_decomp_file_prefix, &
         config_number_of_blocks, &
         config_explicit_proc_decomp, &
         config_proc_decomp_file_prefix
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, decomposition, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_char(dminfo, config_block_decomp_file_prefix)
         call mpas_dmpar_bcast_int(dminfo, config_number_of_blocks)
         call mpas_dmpar_bcast_logical(dminfo, config_explicit_proc_decomp)
         call mpas_dmpar_bcast_char(dminfo, config_proc_decomp_file_prefix)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record decomposition')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_block_decomp_file_prefix = '&
//''//mpas_log_escape_dollars(config_block_decomp_file_prefix))
            call mpas_log_write('        config_number_of_blocks = $i', intArgs=(/config_number_of_blocks/))
            call mpas_log_write('        config_explicit_proc_decomp = $l', logicArgs=(/config_explicit_proc_decomp/))
            call mpas_log_write('        config_proc_decomp_file_prefix = '&
//''//mpas_log_escape_dollars(config_proc_decomp_file_prefix))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record decomposition.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_block_decomp_file_prefix', config_block_decomp_file_prefix)
      call mpas_pool_add_config(configPool, 'config_number_of_blocks', config_number_of_blocks)
      call mpas_pool_add_config(configPool, 'config_explicit_proc_decomp', config_explicit_proc_decomp)
      call mpas_pool_add_config(configPool, 'config_proc_decomp_file_prefix', config_proc_decomp_file_prefix)

   end subroutine init_atm_setup_nmlrec_decomposition


