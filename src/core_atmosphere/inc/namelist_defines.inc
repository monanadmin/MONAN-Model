   subroutine atm_setup_nmlrec_nhyd_model(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      character (len=StrKIND) :: config_time_integration = 'SRK3'
      integer :: config_time_integration_order = 2
      real (kind=RKIND) :: config_dt = 720.000000
      character (len=StrKIND) :: config_calendar_type = 'gregorian'
      character (len=StrKIND) :: config_start_time = '2010-10-23_00:00:00'
      character (len=StrKIND) :: config_stop_time = 'none'
      character (len=StrKIND) :: config_run_duration = '5_00:00:00'
      logical :: config_split_dynamics_transport = .true.
      integer :: config_number_of_sub_steps = 2
      integer :: config_dynamics_split_steps = 3
      real (kind=RKIND) :: config_h_mom_eddy_visc2 = 0.000000
      real (kind=RKIND) :: config_h_mom_eddy_visc4 = 0.000000
      real (kind=RKIND) :: config_v_mom_eddy_visc2 = 0.000000
      real (kind=RKIND) :: config_h_theta_eddy_visc2 = 0.000000
      real (kind=RKIND) :: config_h_theta_eddy_visc4 = 0.000000
      real (kind=RKIND) :: config_v_theta_eddy_visc2 = 0.000000
      character (len=StrKIND) :: config_horiz_mixing = '2d_smagorinsky'
      real (kind=RKIND) :: config_len_disp = 0.000000
      real (kind=RKIND) :: config_visc4_2dsmag = 0.050000
      real (kind=RKIND) :: config_del4u_div_factor = 10.000000
      integer :: config_w_adv_order = 3
      integer :: config_theta_adv_order = 3
      integer :: config_scalar_adv_order = 3
      integer :: config_u_vadv_order = 3
      integer :: config_w_vadv_order = 3
      integer :: config_theta_vadv_order = 3
      integer :: config_scalar_vadv_order = 3
      logical :: config_scalar_advection = .true.
      logical :: config_positive_definite = .false.
      logical :: config_monotonic = .true.
      real (kind=RKIND) :: config_coef_3rd_order = 0.250000
      real (kind=RKIND) :: config_smagorinsky_coef = 0.125000
      logical :: config_mix_full = .true.
      real (kind=RKIND) :: config_epssm = 0.100000
      real (kind=RKIND) :: config_smdiv = 0.100000
      real (kind=RKIND) :: config_apvm_upwinding = 0.500000
      logical :: config_h_ScaleWithMesh = .true.
      integer :: config_num_halos = 2
      real (kind=RKIND) :: config_relax_zone_divdamp_coef = 6.000000

      namelist /nhyd_model/ &
         config_time_integration, &
         config_time_integration_order, &
         config_dt, &
         config_calendar_type, &
         config_start_time, &
         config_stop_time, &
         config_run_duration, &
         config_split_dynamics_transport, &
         config_number_of_sub_steps, &
         config_dynamics_split_steps, &
         config_h_mom_eddy_visc2, &
         config_h_mom_eddy_visc4, &
         config_v_mom_eddy_visc2, &
         config_h_theta_eddy_visc2, &
         config_h_theta_eddy_visc4, &
         config_v_theta_eddy_visc2, &
         config_horiz_mixing, &
         config_len_disp, &
         config_visc4_2dsmag, &
         config_del4u_div_factor, &
         config_w_adv_order, &
         config_theta_adv_order, &
         config_scalar_adv_order, &
         config_u_vadv_order, &
         config_w_vadv_order, &
         config_theta_vadv_order, &
         config_scalar_vadv_order, &
         config_scalar_advection, &
         config_positive_definite, &
         config_monotonic, &
         config_coef_3rd_order, &
         config_smagorinsky_coef, &
         config_mix_full, &
         config_epssm, &
         config_smdiv, &
         config_apvm_upwinding, &
         config_h_ScaleWithMesh, &
         config_num_halos, &
         config_relax_zone_divdamp_coef
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
         call mpas_dmpar_bcast_char(dminfo, config_time_integration)
         call mpas_dmpar_bcast_int(dminfo, config_time_integration_order)
         call mpas_dmpar_bcast_real(dminfo, config_dt)
         call mpas_dmpar_bcast_char(dminfo, config_calendar_type)
         call mpas_dmpar_bcast_char(dminfo, config_start_time)
         call mpas_dmpar_bcast_char(dminfo, config_stop_time)
         call mpas_dmpar_bcast_char(dminfo, config_run_duration)
         call mpas_dmpar_bcast_logical(dminfo, config_split_dynamics_transport)
         call mpas_dmpar_bcast_int(dminfo, config_number_of_sub_steps)
         call mpas_dmpar_bcast_int(dminfo, config_dynamics_split_steps)
         call mpas_dmpar_bcast_real(dminfo, config_h_mom_eddy_visc2)
         call mpas_dmpar_bcast_real(dminfo, config_h_mom_eddy_visc4)
         call mpas_dmpar_bcast_real(dminfo, config_v_mom_eddy_visc2)
         call mpas_dmpar_bcast_real(dminfo, config_h_theta_eddy_visc2)
         call mpas_dmpar_bcast_real(dminfo, config_h_theta_eddy_visc4)
         call mpas_dmpar_bcast_real(dminfo, config_v_theta_eddy_visc2)
         call mpas_dmpar_bcast_char(dminfo, config_horiz_mixing)
         call mpas_dmpar_bcast_real(dminfo, config_len_disp)
         call mpas_dmpar_bcast_real(dminfo, config_visc4_2dsmag)
         call mpas_dmpar_bcast_real(dminfo, config_del4u_div_factor)
         call mpas_dmpar_bcast_int(dminfo, config_w_adv_order)
         call mpas_dmpar_bcast_int(dminfo, config_theta_adv_order)
         call mpas_dmpar_bcast_int(dminfo, config_scalar_adv_order)
         call mpas_dmpar_bcast_int(dminfo, config_u_vadv_order)
         call mpas_dmpar_bcast_int(dminfo, config_w_vadv_order)
         call mpas_dmpar_bcast_int(dminfo, config_theta_vadv_order)
         call mpas_dmpar_bcast_int(dminfo, config_scalar_vadv_order)
         call mpas_dmpar_bcast_logical(dminfo, config_scalar_advection)
         call mpas_dmpar_bcast_logical(dminfo, config_positive_definite)
         call mpas_dmpar_bcast_logical(dminfo, config_monotonic)
         call mpas_dmpar_bcast_real(dminfo, config_coef_3rd_order)
         call mpas_dmpar_bcast_real(dminfo, config_smagorinsky_coef)
         call mpas_dmpar_bcast_logical(dminfo, config_mix_full)
         call mpas_dmpar_bcast_real(dminfo, config_epssm)
         call mpas_dmpar_bcast_real(dminfo, config_smdiv)
         call mpas_dmpar_bcast_real(dminfo, config_apvm_upwinding)
         call mpas_dmpar_bcast_logical(dminfo, config_h_ScaleWithMesh)
         call mpas_dmpar_bcast_int(dminfo, config_num_halos)
         call mpas_dmpar_bcast_real(dminfo, config_relax_zone_divdamp_coef)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record nhyd_model')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_time_integration = '//mpas_log_escape_dollars(config_time_integration))
            call mpas_log_write('        config_time_integration_order = $i', intArgs=(/config_time_integration_order/))
            call mpas_log_write('        config_dt = $r', realArgs=(/config_dt/))
            call mpas_log_write('        config_calendar_type = '//mpas_log_escape_dollars(config_calendar_type))
            call mpas_log_write('        config_start_time = '//mpas_log_escape_dollars(config_start_time))
            call mpas_log_write('        config_stop_time = '//mpas_log_escape_dollars(config_stop_time))
            call mpas_log_write('        config_run_duration = '//mpas_log_escape_dollars(config_run_duration))
            call mpas_log_write('        config_split_dynamics_transport = $l', logicArgs=(/config_split_dynamics_transport/))
            call mpas_log_write('        config_number_of_sub_steps = $i', intArgs=(/config_number_of_sub_steps/))
            call mpas_log_write('        config_dynamics_split_steps = $i', intArgs=(/config_dynamics_split_steps/))
            call mpas_log_write('        config_h_mom_eddy_visc2 = $r', realArgs=(/config_h_mom_eddy_visc2/))
            call mpas_log_write('        config_h_mom_eddy_visc4 = $r', realArgs=(/config_h_mom_eddy_visc4/))
            call mpas_log_write('        config_v_mom_eddy_visc2 = $r', realArgs=(/config_v_mom_eddy_visc2/))
            call mpas_log_write('        config_h_theta_eddy_visc2 = $r', realArgs=(/config_h_theta_eddy_visc2/))
            call mpas_log_write('        config_h_theta_eddy_visc4 = $r', realArgs=(/config_h_theta_eddy_visc4/))
            call mpas_log_write('        config_v_theta_eddy_visc2 = $r', realArgs=(/config_v_theta_eddy_visc2/))
            call mpas_log_write('        config_horiz_mixing = '//mpas_log_escape_dollars(config_horiz_mixing))
            call mpas_log_write('        config_len_disp = $r', realArgs=(/config_len_disp/))
            call mpas_log_write('        config_visc4_2dsmag = $r', realArgs=(/config_visc4_2dsmag/))
            call mpas_log_write('        config_del4u_div_factor = $r', realArgs=(/config_del4u_div_factor/))
            call mpas_log_write('        config_w_adv_order = $i', intArgs=(/config_w_adv_order/))
            call mpas_log_write('        config_theta_adv_order = $i', intArgs=(/config_theta_adv_order/))
            call mpas_log_write('        config_scalar_adv_order = $i', intArgs=(/config_scalar_adv_order/))
            call mpas_log_write('        config_u_vadv_order = $i', intArgs=(/config_u_vadv_order/))
            call mpas_log_write('        config_w_vadv_order = $i', intArgs=(/config_w_vadv_order/))
            call mpas_log_write('        config_theta_vadv_order = $i', intArgs=(/config_theta_vadv_order/))
            call mpas_log_write('        config_scalar_vadv_order = $i', intArgs=(/config_scalar_vadv_order/))
            call mpas_log_write('        config_scalar_advection = $l', logicArgs=(/config_scalar_advection/))
            call mpas_log_write('        config_positive_definite = $l', logicArgs=(/config_positive_definite/))
            call mpas_log_write('        config_monotonic = $l', logicArgs=(/config_monotonic/))
            call mpas_log_write('        config_coef_3rd_order = $r', realArgs=(/config_coef_3rd_order/))
            call mpas_log_write('        config_smagorinsky_coef = $r', realArgs=(/config_smagorinsky_coef/))
            call mpas_log_write('        config_mix_full = $l', logicArgs=(/config_mix_full/))
            call mpas_log_write('        config_epssm = $r', realArgs=(/config_epssm/))
            call mpas_log_write('        config_smdiv = $r', realArgs=(/config_smdiv/))
            call mpas_log_write('        config_apvm_upwinding = $r', realArgs=(/config_apvm_upwinding/))
            call mpas_log_write('        config_h_ScaleWithMesh = $l', logicArgs=(/config_h_ScaleWithMesh/))
            call mpas_log_write('        config_num_halos = $i', intArgs=(/config_num_halos/))
            call mpas_log_write('        config_relax_zone_divdamp_coef = $r', realArgs=(/config_relax_zone_divdamp_coef/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record nhyd_model.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_time_integration', config_time_integration)
      call mpas_pool_add_config(configPool, 'config_time_integration_order', config_time_integration_order)
      call mpas_pool_add_config(configPool, 'config_dt', config_dt)
      call mpas_pool_add_config(configPool, 'config_calendar_type', config_calendar_type)
      call mpas_pool_add_config(configPool, 'config_start_time', config_start_time)
      call mpas_pool_add_config(configPool, 'config_stop_time', config_stop_time)
      call mpas_pool_add_config(configPool, 'config_run_duration', config_run_duration)
      call mpas_pool_add_config(configPool, 'config_split_dynamics_transport', config_split_dynamics_transport)
      call mpas_pool_add_config(configPool, 'config_number_of_sub_steps', config_number_of_sub_steps)
      call mpas_pool_add_config(configPool, 'config_dynamics_split_steps', config_dynamics_split_steps)
      call mpas_pool_add_config(configPool, 'config_h_mom_eddy_visc2', config_h_mom_eddy_visc2)
      call mpas_pool_add_config(configPool, 'config_h_mom_eddy_visc4', config_h_mom_eddy_visc4)
      call mpas_pool_add_config(configPool, 'config_v_mom_eddy_visc2', config_v_mom_eddy_visc2)
      call mpas_pool_add_config(configPool, 'config_h_theta_eddy_visc2', config_h_theta_eddy_visc2)
      call mpas_pool_add_config(configPool, 'config_h_theta_eddy_visc4', config_h_theta_eddy_visc4)
      call mpas_pool_add_config(configPool, 'config_v_theta_eddy_visc2', config_v_theta_eddy_visc2)
      call mpas_pool_add_config(configPool, 'config_horiz_mixing', config_horiz_mixing)
      call mpas_pool_add_config(configPool, 'config_len_disp', config_len_disp)
      call mpas_pool_add_config(configPool, 'config_visc4_2dsmag', config_visc4_2dsmag)
      call mpas_pool_add_config(configPool, 'config_del4u_div_factor', config_del4u_div_factor)
      call mpas_pool_add_config(configPool, 'config_w_adv_order', config_w_adv_order)
      call mpas_pool_add_config(configPool, 'config_theta_adv_order', config_theta_adv_order)
      call mpas_pool_add_config(configPool, 'config_scalar_adv_order', config_scalar_adv_order)
      call mpas_pool_add_config(configPool, 'config_u_vadv_order', config_u_vadv_order)
      call mpas_pool_add_config(configPool, 'config_w_vadv_order', config_w_vadv_order)
      call mpas_pool_add_config(configPool, 'config_theta_vadv_order', config_theta_vadv_order)
      call mpas_pool_add_config(configPool, 'config_scalar_vadv_order', config_scalar_vadv_order)
      call mpas_pool_add_config(configPool, 'config_scalar_advection', config_scalar_advection)
      call mpas_pool_add_config(configPool, 'config_positive_definite', config_positive_definite)
      call mpas_pool_add_config(configPool, 'config_monotonic', config_monotonic)
      call mpas_pool_add_config(configPool, 'config_coef_3rd_order', config_coef_3rd_order)
      call mpas_pool_add_config(configPool, 'config_smagorinsky_coef', config_smagorinsky_coef)
      call mpas_pool_add_config(configPool, 'config_mix_full', config_mix_full)
      call mpas_pool_add_config(configPool, 'config_epssm', config_epssm)
      call mpas_pool_add_config(configPool, 'config_smdiv', config_smdiv)
      call mpas_pool_add_config(configPool, 'config_apvm_upwinding', config_apvm_upwinding)
      call mpas_pool_add_config(configPool, 'config_h_ScaleWithMesh', config_h_ScaleWithMesh)
      call mpas_pool_add_config(configPool, 'config_num_halos', config_num_halos)
      call mpas_pool_add_config(configPool, 'config_relax_zone_divdamp_coef', config_relax_zone_divdamp_coef)

   end subroutine atm_setup_nmlrec_nhyd_model


   subroutine atm_setup_nmlrec_damping(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      real (kind=RKIND) :: config_zd = 22000.000000
      real (kind=RKIND) :: config_xnutr = 0.200000
      real (kind=RKIND) :: config_mpas_cam_coef = 0.000000
      integer :: config_number_cam_damping_levels = 4
      logical :: config_rayleigh_damp_u = .false.
      real (kind=RKIND) :: config_rayleigh_damp_u_timescale_days = 5.000000
      integer :: config_number_rayleigh_damp_u_levels = 6

      namelist /damping/ &
         config_zd, &
         config_xnutr, &
         config_mpas_cam_coef, &
         config_number_cam_damping_levels, &
         config_rayleigh_damp_u, &
         config_rayleigh_damp_u_timescale_days, &
         config_number_rayleigh_damp_u_levels
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, damping, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_real(dminfo, config_zd)
         call mpas_dmpar_bcast_real(dminfo, config_xnutr)
         call mpas_dmpar_bcast_real(dminfo, config_mpas_cam_coef)
         call mpas_dmpar_bcast_int(dminfo, config_number_cam_damping_levels)
         call mpas_dmpar_bcast_logical(dminfo, config_rayleigh_damp_u)
         call mpas_dmpar_bcast_real(dminfo, config_rayleigh_damp_u_timescale_days)
         call mpas_dmpar_bcast_int(dminfo, config_number_rayleigh_damp_u_levels)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record damping')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_zd = $r', realArgs=(/config_zd/))
            call mpas_log_write('        config_xnutr = $r', realArgs=(/config_xnutr/))
            call mpas_log_write('        config_mpas_cam_coef = $r', realArgs=(/config_mpas_cam_coef/))
            call mpas_log_write('        config_number_cam_damping_levels = $i', intArgs=(/config_number_cam_damping_levels/))
            call mpas_log_write('        config_rayleigh_damp_u = $l', logicArgs=(/config_rayleigh_damp_u/))
            call mpas_log_write('        config_rayleigh_damp_u_timescale_days = $r', &
realArgs=(/config_rayleigh_damp_u_timescale_days/))
            call mpas_log_write('        config_number_rayleigh_damp_u_levels = $i', &
intArgs=(/config_number_rayleigh_damp_u_levels/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record damping.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_zd', config_zd)
      call mpas_pool_add_config(configPool, 'config_xnutr', config_xnutr)
      call mpas_pool_add_config(configPool, 'config_mpas_cam_coef', config_mpas_cam_coef)
      call mpas_pool_add_config(configPool, 'config_number_cam_damping_levels', config_number_cam_damping_levels)
      call mpas_pool_add_config(configPool, 'config_rayleigh_damp_u', config_rayleigh_damp_u)
      call mpas_pool_add_config(configPool, 'config_rayleigh_damp_u_timescale_days', config_rayleigh_damp_u_timescale_days)
      call mpas_pool_add_config(configPool, 'config_number_rayleigh_damp_u_levels', config_number_rayleigh_damp_u_levels)

   end subroutine atm_setup_nmlrec_damping


   subroutine atm_setup_nmlrec_limited_area(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      logical :: config_apply_lbcs = .false.

      namelist /limited_area/ &
         config_apply_lbcs
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, limited_area, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_logical(dminfo, config_apply_lbcs)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record limited_area')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_apply_lbcs = $l', logicArgs=(/config_apply_lbcs/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record limited_area.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_apply_lbcs', config_apply_lbcs)

   end subroutine atm_setup_nmlrec_limited_area


   subroutine atm_setup_nmlrec_io(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      character (len=StrKIND) :: config_restart_timestamp_name = 'restart_timestamp'
      integer :: config_pio_num_iotasks = 0
      integer :: config_pio_stride = 1

      namelist /io/ &
         config_restart_timestamp_name, &
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
         call mpas_dmpar_bcast_char(dminfo, config_restart_timestamp_name)
         call mpas_dmpar_bcast_int(dminfo, config_pio_num_iotasks)
         call mpas_dmpar_bcast_int(dminfo, config_pio_stride)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record io')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_restart_timestamp_name = '//mpas_log_escape_dollars(config_restart_timestamp_name))
            call mpas_log_write('        config_pio_num_iotasks = $i', intArgs=(/config_pio_num_iotasks/))
            call mpas_log_write('        config_pio_stride = $i', intArgs=(/config_pio_stride/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record io.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_restart_timestamp_name', config_restart_timestamp_name)
      call mpas_pool_add_config(configPool, 'config_pio_num_iotasks', config_pio_num_iotasks)
      call mpas_pool_add_config(configPool, 'config_pio_stride', config_pio_stride)

   end subroutine atm_setup_nmlrec_io


   subroutine atm_setup_nmlrec_decomposition(configPool, unitNumber, dminfo)
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

   end subroutine atm_setup_nmlrec_decomposition


   subroutine atm_setup_nmlrec_restart(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      logical :: config_do_restart = .false.
      logical :: config_do_DAcycling = .false.

      namelist /restart/ &
         config_do_restart, &
         config_do_DAcycling
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, restart, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_logical(dminfo, config_do_restart)
         call mpas_dmpar_bcast_logical(dminfo, config_do_DAcycling)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record restart')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_do_restart = $l', logicArgs=(/config_do_restart/))
            call mpas_log_write('        config_do_DAcycling = $l', logicArgs=(/config_do_DAcycling/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record restart.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_do_restart', config_do_restart)
      call mpas_pool_add_config(configPool, 'config_do_DAcycling', config_do_DAcycling)

   end subroutine atm_setup_nmlrec_restart


   subroutine atm_setup_nmlrec_printout(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      logical :: config_print_global_minmax_vel = .true.
      logical :: config_print_detailed_minmax_vel = .false.
      logical :: config_print_global_minmax_sca = .false.

      namelist /printout/ &
         config_print_global_minmax_vel, &
         config_print_detailed_minmax_vel, &
         config_print_global_minmax_sca
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, printout, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_logical(dminfo, config_print_global_minmax_vel)
         call mpas_dmpar_bcast_logical(dminfo, config_print_detailed_minmax_vel)
         call mpas_dmpar_bcast_logical(dminfo, config_print_global_minmax_sca)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record printout')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_print_global_minmax_vel = $l', logicArgs=(/config_print_global_minmax_vel/))
            call mpas_log_write('        config_print_detailed_minmax_vel = $l', logicArgs=(/config_print_detailed_minmax_vel/))
            call mpas_log_write('        config_print_global_minmax_sca = $l', logicArgs=(/config_print_global_minmax_sca/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record printout.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_print_global_minmax_vel', config_print_global_minmax_vel)
      call mpas_pool_add_config(configPool, 'config_print_detailed_minmax_vel', config_print_detailed_minmax_vel)
      call mpas_pool_add_config(configPool, 'config_print_global_minmax_sca', config_print_global_minmax_sca)

   end subroutine atm_setup_nmlrec_printout


   subroutine atm_setup_nmlrec_IAU(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      character (len=StrKIND) :: config_IAU_option = 'off'
      real (kind=RKIND) :: config_IAU_window_length_s = 21600.000000

      namelist /IAU/ &
         config_IAU_option, &
         config_IAU_window_length_s
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, IAU, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_char(dminfo, config_IAU_option)
         call mpas_dmpar_bcast_real(dminfo, config_IAU_window_length_s)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record IAU')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_IAU_option = '//mpas_log_escape_dollars(config_IAU_option))
            call mpas_log_write('        config_IAU_window_length_s = $r', realArgs=(/config_IAU_window_length_s/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record IAU.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_IAU_option', config_IAU_option)
      call mpas_pool_add_config(configPool, 'config_IAU_window_length_s', config_IAU_window_length_s)

   end subroutine atm_setup_nmlrec_IAU


   subroutine atm_setup_nmlrec_assimilation(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      logical :: config_jedi_da = .false.

      namelist /assimilation/ &
         config_jedi_da
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, assimilation, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_logical(dminfo, config_jedi_da)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record assimilation')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_jedi_da = $l', logicArgs=(/config_jedi_da/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record assimilation.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_jedi_da', config_jedi_da)

   end subroutine atm_setup_nmlrec_assimilation


   subroutine atm_setup_nmlrec_development(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      character (len=StrKIND) :: config_halo_exch_method = 'mpas_halo'

      namelist /development/ &
         config_halo_exch_method
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, development, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_char(dminfo, config_halo_exch_method)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record development')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_halo_exch_method = '//mpas_log_escape_dollars(config_halo_exch_method))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record development.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_halo_exch_method', config_halo_exch_method)

   end subroutine atm_setup_nmlrec_development


   subroutine atm_setup_nmlrec_physics(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      character (len=StrKIND) :: input_soil_data = 'STAS'
      integer :: input_soil_temperature_lag = 140
      integer :: num_soil_layers = 4
      integer :: months = 12
      integer :: noznlev = 59
      integer :: naerlev = 29
      integer :: camdim1 = 4
      logical :: config_frac_seaice = .true.
      logical :: config_sfc_albedo = .true.
      logical :: config_sfc_snowalbedo = .true.
      logical :: config_sst_update = .false.
      logical :: config_sstdiurn_update = .false.
      logical :: config_deepsoiltemp_update = .false.
      logical :: config_o3climatology = .true.
      logical :: config_microp_re = .false.
      logical :: config_ysu_pblmix = .false.
      integer :: config_n_microp = 1
      character (len=StrKIND) :: config_radtlw_interval = '00:30:00'
      character (len=StrKIND) :: config_radtsw_interval = '00:30:00'
      character (len=StrKIND) :: config_conv_interval = 'none'
      character (len=StrKIND) :: config_pbl_interval = 'none'
      character (len=StrKIND) :: config_camrad_abs_update = '06:00:00'
      character (len=StrKIND) :: config_greeness_update = '24:00:00'
      character (len=StrKIND) :: config_bucket_update = 'none'
      character (len=StrKIND) :: config_physics_suite = 'mesoscale_reference'
      character (len=StrKIND) :: config_microp_scheme = 'suite'
      character (len=StrKIND) :: config_convection_scheme = 'suite'
      character (len=StrKIND) :: config_lsm_scheme = 'suite'
      character (len=StrKIND) :: config_pbl_scheme = 'suite'
      character (len=StrKIND) :: config_gwdo_scheme = 'suite'
      character (len=StrKIND) :: config_radt_cld_scheme = 'suite'
      character (len=StrKIND) :: config_radt_lw_scheme = 'suite'
      character (len=StrKIND) :: config_radt_sw_scheme = 'suite'
      character (len=StrKIND) :: config_sfclayer_scheme = 'suite'
      integer :: config_gfconv_closure_deep = 0
      integer :: config_gfconv_closure_shallow = 8
      real (kind=RKIND) :: config_bucket_radt = 1000000000.000000
      real (kind=RKIND) :: config_bucket_rainc = 100.000000
      real (kind=RKIND) :: config_bucket_rainnc = 100.000000
      logical :: config_oml1d = .false.
      real (kind=RKIND) :: config_oml_hml0 = 30.000000
      real (kind=RKIND) :: config_oml_gamma = 0.140000
      real (kind=RKIND) :: config_oml_relaxation_time = 864000.000000

      namelist /physics/ &
         input_soil_data, &
         input_soil_temperature_lag, &
         num_soil_layers, &
         months, &
         noznlev, &
         naerlev, &
         camdim1, &
         config_frac_seaice, &
         config_sfc_albedo, &
         config_sfc_snowalbedo, &
         config_sst_update, &
         config_sstdiurn_update, &
         config_deepsoiltemp_update, &
         config_o3climatology, &
         config_microp_re, &
         config_ysu_pblmix, &
         config_n_microp, &
         config_radtlw_interval, &
         config_radtsw_interval, &
         config_conv_interval, &
         config_pbl_interval, &
         config_camrad_abs_update, &
         config_greeness_update, &
         config_bucket_update, &
         config_physics_suite, &
         config_microp_scheme, &
         config_convection_scheme, &
         config_lsm_scheme, &
         config_pbl_scheme, &
         config_gwdo_scheme, &
         config_radt_cld_scheme, &
         config_radt_lw_scheme, &
         config_radt_sw_scheme, &
         config_sfclayer_scheme, &
         config_gfconv_closure_deep, &
         config_gfconv_closure_shallow, &
         config_bucket_radt, &
         config_bucket_rainc, &
         config_bucket_rainnc, &
         config_oml1d, &
         config_oml_hml0, &
         config_oml_gamma, &
         config_oml_relaxation_time
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
         call mpas_dmpar_bcast_char(dminfo, input_soil_data)
         call mpas_dmpar_bcast_int(dminfo, input_soil_temperature_lag)
         call mpas_dmpar_bcast_int(dminfo, num_soil_layers)
         call mpas_dmpar_bcast_int(dminfo, months)
         call mpas_dmpar_bcast_int(dminfo, noznlev)
         call mpas_dmpar_bcast_int(dminfo, naerlev)
         call mpas_dmpar_bcast_int(dminfo, camdim1)
         call mpas_dmpar_bcast_logical(dminfo, config_frac_seaice)
         call mpas_dmpar_bcast_logical(dminfo, config_sfc_albedo)
         call mpas_dmpar_bcast_logical(dminfo, config_sfc_snowalbedo)
         call mpas_dmpar_bcast_logical(dminfo, config_sst_update)
         call mpas_dmpar_bcast_logical(dminfo, config_sstdiurn_update)
         call mpas_dmpar_bcast_logical(dminfo, config_deepsoiltemp_update)
         call mpas_dmpar_bcast_logical(dminfo, config_o3climatology)
         call mpas_dmpar_bcast_logical(dminfo, config_microp_re)
         call mpas_dmpar_bcast_logical(dminfo, config_ysu_pblmix)
         call mpas_dmpar_bcast_int(dminfo, config_n_microp)
         call mpas_dmpar_bcast_char(dminfo, config_radtlw_interval)
         call mpas_dmpar_bcast_char(dminfo, config_radtsw_interval)
         call mpas_dmpar_bcast_char(dminfo, config_conv_interval)
         call mpas_dmpar_bcast_char(dminfo, config_pbl_interval)
         call mpas_dmpar_bcast_char(dminfo, config_camrad_abs_update)
         call mpas_dmpar_bcast_char(dminfo, config_greeness_update)
         call mpas_dmpar_bcast_char(dminfo, config_bucket_update)
         call mpas_dmpar_bcast_char(dminfo, config_physics_suite)
         call mpas_dmpar_bcast_char(dminfo, config_microp_scheme)
         call mpas_dmpar_bcast_char(dminfo, config_convection_scheme)
         call mpas_dmpar_bcast_char(dminfo, config_lsm_scheme)
         call mpas_dmpar_bcast_char(dminfo, config_pbl_scheme)
         call mpas_dmpar_bcast_char(dminfo, config_gwdo_scheme)
         call mpas_dmpar_bcast_char(dminfo, config_radt_cld_scheme)
         call mpas_dmpar_bcast_char(dminfo, config_radt_lw_scheme)
         call mpas_dmpar_bcast_char(dminfo, config_radt_sw_scheme)
         call mpas_dmpar_bcast_char(dminfo, config_sfclayer_scheme)
         call mpas_dmpar_bcast_int(dminfo, config_gfconv_closure_deep)
         call mpas_dmpar_bcast_int(dminfo, config_gfconv_closure_shallow)
         call mpas_dmpar_bcast_real(dminfo, config_bucket_radt)
         call mpas_dmpar_bcast_real(dminfo, config_bucket_rainc)
         call mpas_dmpar_bcast_real(dminfo, config_bucket_rainnc)
         call mpas_dmpar_bcast_logical(dminfo, config_oml1d)
         call mpas_dmpar_bcast_real(dminfo, config_oml_hml0)
         call mpas_dmpar_bcast_real(dminfo, config_oml_gamma)
         call mpas_dmpar_bcast_real(dminfo, config_oml_relaxation_time)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record physics')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        input_soil_data = '//mpas_log_escape_dollars(input_soil_data))
            call mpas_log_write('        input_soil_temperature_lag = $i', intArgs=(/input_soil_temperature_lag/))
            call mpas_log_write('        num_soil_layers = $i', intArgs=(/num_soil_layers/))
            call mpas_log_write('        months = $i', intArgs=(/months/))
            call mpas_log_write('        noznlev = $i', intArgs=(/noznlev/))
            call mpas_log_write('        naerlev = $i', intArgs=(/naerlev/))
            call mpas_log_write('        camdim1 = $i', intArgs=(/camdim1/))
            call mpas_log_write('        config_frac_seaice = $l', logicArgs=(/config_frac_seaice/))
            call mpas_log_write('        config_sfc_albedo = $l', logicArgs=(/config_sfc_albedo/))
            call mpas_log_write('        config_sfc_snowalbedo = $l', logicArgs=(/config_sfc_snowalbedo/))
            call mpas_log_write('        config_sst_update = $l', logicArgs=(/config_sst_update/))
            call mpas_log_write('        config_sstdiurn_update = $l', logicArgs=(/config_sstdiurn_update/))
            call mpas_log_write('        config_deepsoiltemp_update = $l', logicArgs=(/config_deepsoiltemp_update/))
            call mpas_log_write('        config_o3climatology = $l', logicArgs=(/config_o3climatology/))
            call mpas_log_write('        config_microp_re = $l', logicArgs=(/config_microp_re/))
            call mpas_log_write('        config_ysu_pblmix = $l', logicArgs=(/config_ysu_pblmix/))
            call mpas_log_write('        config_n_microp = $i', intArgs=(/config_n_microp/))
            call mpas_log_write('        config_radtlw_interval = '//mpas_log_escape_dollars(config_radtlw_interval))
            call mpas_log_write('        config_radtsw_interval = '//mpas_log_escape_dollars(config_radtsw_interval))
            call mpas_log_write('        config_conv_interval = '//mpas_log_escape_dollars(config_conv_interval))
            call mpas_log_write('        config_pbl_interval = '//mpas_log_escape_dollars(config_pbl_interval))
            call mpas_log_write('        config_camrad_abs_update = '//mpas_log_escape_dollars(config_camrad_abs_update))
            call mpas_log_write('        config_greeness_update = '//mpas_log_escape_dollars(config_greeness_update))
            call mpas_log_write('        config_bucket_update = '//mpas_log_escape_dollars(config_bucket_update))
            call mpas_log_write('        config_physics_suite = '//mpas_log_escape_dollars(config_physics_suite))
            call mpas_log_write('        config_microp_scheme = '//mpas_log_escape_dollars(config_microp_scheme))
            call mpas_log_write('        config_convection_scheme = '//mpas_log_escape_dollars(config_convection_scheme))
            call mpas_log_write('        config_lsm_scheme = '//mpas_log_escape_dollars(config_lsm_scheme))
            call mpas_log_write('        config_pbl_scheme = '//mpas_log_escape_dollars(config_pbl_scheme))
            call mpas_log_write('        config_gwdo_scheme = '//mpas_log_escape_dollars(config_gwdo_scheme))
            call mpas_log_write('        config_radt_cld_scheme = '//mpas_log_escape_dollars(config_radt_cld_scheme))
            call mpas_log_write('        config_radt_lw_scheme = '//mpas_log_escape_dollars(config_radt_lw_scheme))
            call mpas_log_write('        config_radt_sw_scheme = '//mpas_log_escape_dollars(config_radt_sw_scheme))
            call mpas_log_write('        config_sfclayer_scheme = '//mpas_log_escape_dollars(config_sfclayer_scheme))
            call mpas_log_write('        config_gfconv_closure_deep = $i', intArgs=(/config_gfconv_closure_deep/))
            call mpas_log_write('        config_gfconv_closure_shallow = $i', intArgs=(/config_gfconv_closure_shallow/))
            call mpas_log_write('        config_bucket_radt = $r', realArgs=(/config_bucket_radt/))
            call mpas_log_write('        config_bucket_rainc = $r', realArgs=(/config_bucket_rainc/))
            call mpas_log_write('        config_bucket_rainnc = $r', realArgs=(/config_bucket_rainnc/))
            call mpas_log_write('        config_oml1d = $l', logicArgs=(/config_oml1d/))
            call mpas_log_write('        config_oml_hml0 = $r', realArgs=(/config_oml_hml0/))
            call mpas_log_write('        config_oml_gamma = $r', realArgs=(/config_oml_gamma/))
            call mpas_log_write('        config_oml_relaxation_time = $r', realArgs=(/config_oml_relaxation_time/))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record physics.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'input_soil_data', input_soil_data)
      call mpas_pool_add_config(configPool, 'input_soil_temperature_lag', input_soil_temperature_lag)
      call mpas_pool_add_config(configPool, 'num_soil_layers', num_soil_layers)
      call mpas_pool_add_config(configPool, 'months', months)
      call mpas_pool_add_config(configPool, 'noznlev', noznlev)
      call mpas_pool_add_config(configPool, 'naerlev', naerlev)
      call mpas_pool_add_config(configPool, 'camdim1', camdim1)
      call mpas_pool_add_config(configPool, 'config_frac_seaice', config_frac_seaice)
      call mpas_pool_add_config(configPool, 'config_sfc_albedo', config_sfc_albedo)
      call mpas_pool_add_config(configPool, 'config_sfc_snowalbedo', config_sfc_snowalbedo)
      call mpas_pool_add_config(configPool, 'config_sst_update', config_sst_update)
      call mpas_pool_add_config(configPool, 'config_sstdiurn_update', config_sstdiurn_update)
      call mpas_pool_add_config(configPool, 'config_deepsoiltemp_update', config_deepsoiltemp_update)
      call mpas_pool_add_config(configPool, 'config_o3climatology', config_o3climatology)
      call mpas_pool_add_config(configPool, 'config_microp_re', config_microp_re)
      call mpas_pool_add_config(configPool, 'config_ysu_pblmix', config_ysu_pblmix)
      call mpas_pool_add_config(configPool, 'config_n_microp', config_n_microp)
      call mpas_pool_add_config(configPool, 'config_radtlw_interval', config_radtlw_interval)
      call mpas_pool_add_config(configPool, 'config_radtsw_interval', config_radtsw_interval)
      call mpas_pool_add_config(configPool, 'config_conv_interval', config_conv_interval)
      call mpas_pool_add_config(configPool, 'config_pbl_interval', config_pbl_interval)
      call mpas_pool_add_config(configPool, 'config_camrad_abs_update', config_camrad_abs_update)
      call mpas_pool_add_config(configPool, 'config_greeness_update', config_greeness_update)
      call mpas_pool_add_config(configPool, 'config_bucket_update', config_bucket_update)
      call mpas_pool_add_config(configPool, 'config_physics_suite', config_physics_suite)
      call mpas_pool_add_config(configPool, 'config_microp_scheme', config_microp_scheme)
      call mpas_pool_add_config(configPool, 'config_convection_scheme', config_convection_scheme)
      call mpas_pool_add_config(configPool, 'config_lsm_scheme', config_lsm_scheme)
      call mpas_pool_add_config(configPool, 'config_pbl_scheme', config_pbl_scheme)
      call mpas_pool_add_config(configPool, 'config_gwdo_scheme', config_gwdo_scheme)
      call mpas_pool_add_config(configPool, 'config_radt_cld_scheme', config_radt_cld_scheme)
      call mpas_pool_add_config(configPool, 'config_radt_lw_scheme', config_radt_lw_scheme)
      call mpas_pool_add_config(configPool, 'config_radt_sw_scheme', config_radt_sw_scheme)
      call mpas_pool_add_config(configPool, 'config_sfclayer_scheme', config_sfclayer_scheme)
      call mpas_pool_add_config(configPool, 'config_gfconv_closure_deep', config_gfconv_closure_deep)
      call mpas_pool_add_config(configPool, 'config_gfconv_closure_shallow', config_gfconv_closure_shallow)
      call mpas_pool_add_config(configPool, 'config_bucket_radt', config_bucket_radt)
      call mpas_pool_add_config(configPool, 'config_bucket_rainc', config_bucket_rainc)
      call mpas_pool_add_config(configPool, 'config_bucket_rainnc', config_bucket_rainnc)
      call mpas_pool_add_config(configPool, 'config_oml1d', config_oml1d)
      call mpas_pool_add_config(configPool, 'config_oml_hml0', config_oml_hml0)
      call mpas_pool_add_config(configPool, 'config_oml_gamma', config_oml_gamma)
      call mpas_pool_add_config(configPool, 'config_oml_relaxation_time', config_oml_relaxation_time)

   end subroutine atm_setup_nmlrec_physics


   subroutine atm_setup_nmlrec_soundings(configPool, unitNumber, dminfo)
      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool
      integer :: ierr

      character (len=StrKIND) :: config_sounding_interval = 'none'

      namelist /soundings/ &
         config_sounding_interval
      if (dminfo % my_proc_id == IO_NODE) then
! Rewinding before each read leads to errors when the code is built with
! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.
#ifndef NAG_COMPILER
         rewind(unitNumber)
#endif
         read(unitNumber, soundings, iostat=ierr)
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr <= 0) then
         call mpas_dmpar_bcast_char(dminfo, config_sounding_interval)
         if (ierr < 0) then
            call mpas_log_write('*** Encountered an issue while attempting to read namelist record soundings')
            call mpas_log_write('    The following values will be used for variables in this record:')
            call mpas_log_write(' ')
            call mpas_log_write('        config_sounding_interval = '//mpas_log_escape_dollars(config_sounding_interval))
            call mpas_log_write(' ')
         end if
      else if (ierr > 0) then
         call mpas_log_write('Error while reading namelist record soundings.', MPAS_LOG_CRIT)
      end if

      call mpas_pool_add_config(configPool, 'config_sounding_interval', config_sounding_interval)

   end subroutine atm_setup_nmlrec_soundings


