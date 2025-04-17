! Copyright (c) 2016,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
module mpas_isobaric_diagnostics

    use mpas_dmpar
    use mpas_kind_types
    use mpas_derived_types
    use mpas_pool_routines
    use mpas_constants
    use mpas_log, only : mpas_log_write
    use VerticalInterpolation

    type (MPAS_pool_type), pointer :: mesh
    type (MPAS_pool_type), pointer :: state
    type (MPAS_pool_type), pointer :: diag
    type (MPAS_pool_type), pointer :: diag_physics
#ifdef DO_PHYSICS
    type (MPAS_pool_type), pointer :: tend_physics
#endif

    type (MPAS_clock_type), pointer :: clock

    public :: isobaric_diagnostics_setup, &
              isobaric_diagnostics_compute

    private

    logical :: need_mslp, need_bl_mynn, need_cu_grell_freitas, need_mp_thompson, &
               need_t_isobaric, need_z_isobaric, need_u_isobaric, need_meanT_500_300

    logical :: need_cldfrac_isobaric, need_temperature_isobaric, need_dewpoint_isobaric, &
               need_relhum_isobaric, need_w_isobaric,need_uzonal_isobaric, &
               need_umeridional_isobaric, need_vorticity_isobaric, need_height_isobaric, &
               need_spechum_isobaric, need_omega_isobaric, &
               need_rqvblten_isobaric, need_rthblten_isobaric, need_rqvcuten_isobaric, &
               need_rthcuten_isobaric, need_rthratenlw_isobaric, need_rthratensw_isobaric, &
               need_rucuten_isobaric, need_rvcuten_isobaric, need_rvblten_isobaric, &
               need_rublten_isobaric, need_kzq_isobaric, need_kzh_isobaric, &
               need_kzm_isobaric, need_ni_isobaric, need_nr_isobaric, &
               need_qr_isobaric, need_qs_isobaric, need_qv_isobaric, &
               need_qc_isobaric, need_qg_isobaric, need_qi_isobaric, &
               need_zgeo_isobaric


    contains


    !-----------------------------------------------------------------------
    !  routine isobaric_diagnostics_setup
    !
    !> \brief Set up the isobaric diagnostics module
    !> \author Michael Duda
    !> \date   21 October 2016
    !> \details
    !>  This routine sets up the isobaric diagnostics module, principally by
    !>  saving pointers to pools that are used in the computation of diagnostics.
    !
    !-----------------------------------------------------------------------
    subroutine isobaric_diagnostics_setup(configs, all_pools, simulation_clock)

        use mpas_derived_types, only : MPAS_pool_type, MPAS_clock_type
        use mpas_pool_routines, only : mpas_pool_get_subpool

        implicit none

        type (MPAS_pool_type), pointer :: configs
        type (MPAS_pool_type), pointer :: all_pools
        type (MPAS_clock_type), pointer :: simulation_clock

        character(len=StrKIND), pointer :: micropScheme, convectionScheme, pblScheme

        clock => simulation_clock

        call mpas_pool_get_subpool(all_pools, 'mesh', mesh)
        call mpas_pool_get_subpool(all_pools, 'state', state)
        call mpas_pool_get_subpool(all_pools, 'diag', diag)
        call mpas_pool_get_subpool(all_pools, 'diag_physics', diag_physics)
#ifdef DO_PHYSICS
        call mpas_pool_get_subpool(all_pools, 'tend_physics', tend_physics)
#endif

        call mpas_pool_get_config(configs, 'config_microp_scheme', micropScheme)
        call mpas_pool_get_config(configs, 'config_convection_scheme', convectionScheme)
        call mpas_pool_get_config(configs, 'config_pbl_scheme', pblScheme)

        if (.not. associated(micropScheme)) then
            call mpas_log_write('config_microp_scheme is not a namelist option...')
            return
        end if

        if (.not. associated(convectionScheme)) then
            call mpas_log_write('config_convection_scheme is not a namelist option...')
            return
        end if

        if (.not. associated(pblScheme)) then
            call mpas_log_write('config_pbl_scheme is not a namelist option...')
            return
        end if        

        if (trim(micropScheme) == "mp_thompson") then
           need_mp_thompson = .true.
        else
           need_mp_thompson = .false.
        endif

        if (trim(convectionScheme) == "cu_grell_freitas") then
           need_cu_grell_freitas = .true.
        else
           need_cu_grell_freitas = .false.
        endif

        if (trim(pblScheme) == "bl_mynn") then
           need_bl_mynn = .true.
        else
           need_bl_mynn = .false.
        endif

!        print *, "mp_thompson=", need_mp_thompson
!        print *, "cu_grell_freitas=", need_cu_grell_freitas
!        print *, "bl_mynn=", need_bl_mynn

    end subroutine isobaric_diagnostics_setup


    !-----------------------------------------------------------------------
    !  routine isobaric_diagnostics_compute
    !
    !> \brief Compute isobaric diagnostic before model output is written
    !> \author Michael Duda
    !> \date   21 October 2016
    !> \details
    !>  Compute isobaric diagnostic before model output is written. Code called
    !>  from here was previously in mpas_atm_interp_diagnostics.F.
    !
    !-----------------------------------------------------------------------
    subroutine isobaric_diagnostics_compute()

        use mpas_atm_diagnostics_utils, only : MPAS_field_will_be_written
 
        implicit none

        logical :: need_any_diags

        need_any_diags = .false.

        need_zgeo_isobaric = .false.
        need_cldfrac_isobaric = .false.
        need_temperature_isobaric = .false.
        need_dewpoint_isobaric = .false.
        need_relhum_isobaric = .false.
        need_w_isobaric = .false.
        need_uzonal_isobaric = .false.
        need_umeridional_isobaric = .false.
        need_vorticity_isobaric = .false.
        need_height_isobaric = .false.
        need_spechum_isobaric = .false.
        need_omega_isobaric = .false.
        need_rqvblten_isobaric = .false.
        need_rthblten_isobaric = .false.
        need_rqvcuten_isobaric = .false.
        need_rthcuten_isobaric = .false.
        need_rthratenlw_isobaric = .false.
        need_rthratensw_isobaric = .false.
        need_rucuten_isobaric = .false.
        need_rvcuten_isobaric = .false.
        need_rvblten_isobaric = .false.
        need_rublten_isobaric = .false.
        need_kzq_isobaric = .false.
        need_kzh_isobaric = .false.
        need_kzm_isobaric = .false.
        need_ni_isobaric = .false.
        need_nr_isobaric = .false.
        need_qr_isobaric = .false.
        need_qs_isobaric = .false.
        need_qv_isobaric = .false.
        need_qc_isobaric = .false.
        need_qg_isobaric = .false.
        need_qi_isobaric = .false.

        need_mslp = MPAS_field_will_be_written('mslp')
        need_any_diags = need_any_diags .or. need_mslp
        

        need_t_isobaric = MPAS_field_will_be_written('t_isobaric')
        need_any_diags = need_any_diags .or. need_t_isobaric
        need_z_isobaric = MPAS_field_will_be_written('z_isobaric')
        need_any_diags = need_any_diags .or. need_z_isobaric
        need_u_isobaric = MPAS_field_will_be_written('uvel_isobaric')
        need_any_diags = need_any_diags .or. need_u_isobaric
        need_meanT_500_300 = MPAS_field_will_be_written('meanT_500_300')
        need_any_diags = need_any_diags .or. need_meanT_500_300
        need_zgeo_isobaric = MPAS_field_will_be_written('zgeo_isobaric')
        need_any_diags = need_any_diags .or. need_zgeo_isobaric


        need_cldfrac_isobaric = MPAS_field_will_be_written('cldfrac_isobaric')
        need_any_diags = need_any_diags .or. need_cldfrac_isobaric
        need_temperature_isobaric = MPAS_field_will_be_written('temperature_isobaric')
        need_any_diags = need_any_diags .or. need_temperature_isobaric
        need_dewpoint_isobaric = MPAS_field_will_be_written('dewpoint_isobaric')
        need_any_diags = need_any_diags .or. need_dewpoint_isobaric
        need_relhum_isobaric = MPAS_field_will_be_written('relhum_isobaric')
        need_any_diags = need_any_diags .or. need_relhum_isobaric
        need_w_isobaric = MPAS_field_will_be_written('w_isobaric')
        need_any_diags = need_any_diags .or. need_w_isobaric
        need_uzonal_isobaric = MPAS_field_will_be_written('uzonal_isobaric')
        need_any_diags = need_any_diags .or. need_uzonal_isobaric
        need_umeridional_isobaric = MPAS_field_will_be_written('umeridional_isobaric')
        need_any_diags = need_any_diags .or. need_umeridional_isobaric
        need_vorticity_isobaric = MPAS_field_will_be_written('vorticity_isobaric')
        need_any_diags = need_any_diags .or. need_vorticity_isobaric
        need_height_isobaric = MPAS_field_will_be_written('height_isobaric')
        need_any_diags = need_any_diags .or. need_height_isobaric
        need_spechum_isobaric = MPAS_field_will_be_written('spechum_isobaric')
        need_any_diags = need_any_diags .or. need_spechum_isobaric
        need_omega_isobaric = MPAS_field_will_be_written('omega_isobaric')
        need_any_diags = need_any_diags .or. need_omega_isobaric
        need_rqvblten_isobaric = MPAS_field_will_be_written('rqvblten_isobaric')
        need_any_diags = need_any_diags .or. need_rqvblten_isobaric
        need_rthblten_isobaric = MPAS_field_will_be_written('rthblten_isobaric')
        need_any_diags = need_any_diags .or. need_rthblten_isobaric
        need_rqvcuten_isobaric = MPAS_field_will_be_written('rqvcuten_isobaric')
        need_any_diags = need_any_diags .or. need_rqvcuten_isobaric
        need_rthcuten_isobaric = MPAS_field_will_be_written('rthcuten_isobaric')
        need_any_diags = need_any_diags .or. need_rthcuten_isobaric
        need_rthratenlw_isobaric = MPAS_field_will_be_written('rthratenlw_isobaric')
        need_any_diags = need_any_diags .or. need_rthratenlw_isobaric
        need_rthratensw_isobaric = MPAS_field_will_be_written('rthratensw_isobaric')
        need_any_diags = need_any_diags .or. need_rthratensw_isobaric
        need_rucuten_isobaric = MPAS_field_will_be_written('rucuten_isobaric')
        need_any_diags = need_any_diags .or. need_rucuten_isobaric
        need_rvcuten_isobaric = MPAS_field_will_be_written('rvcuten_isobaric')
        need_any_diags = need_any_diags .or. need_rvcuten_isobaric
        need_rvblten_isobaric = MPAS_field_will_be_written('rvblten_isobaric')
        need_any_diags = need_any_diags .or. need_rvblten_isobaric
        need_rublten_isobaric = MPAS_field_will_be_written('rublten_isobaric')
        need_any_diags = need_any_diags .or. need_rublten_isobaric
        need_kzq_isobaric = MPAS_field_will_be_written('kzq_isobaric')
        need_any_diags = need_any_diags .or. need_kzq_isobaric
        need_kzh_isobaric = MPAS_field_will_be_written('kzh_isobaric')
        need_any_diags = need_any_diags .or. need_kzh_isobaric
        need_kzm_isobaric = MPAS_field_will_be_written('kzm_isobaric')
        need_any_diags = need_any_diags .or. need_kzm_isobaric
        need_ni_isobaric = MPAS_field_will_be_written('ni_isobaric')
        need_any_diags = need_any_diags .or. need_ni_isobaric
        need_nr_isobaric = MPAS_field_will_be_written('nr_isobaric')
        need_any_diags = need_any_diags .or. need_nr_isobaric
        need_qr_isobaric = MPAS_field_will_be_written('qr_isobaric')
        need_any_diags = need_any_diags .or. need_qr_isobaric
        need_qs_isobaric = MPAS_field_will_be_written('qs_isobaric')
        need_any_diags = need_any_diags .or. need_qs_isobaric
        need_qv_isobaric = MPAS_field_will_be_written('qv_isobaric')
        need_any_diags = need_any_diags .or. need_qv_isobaric
        need_qc_isobaric = MPAS_field_will_be_written('qc_isobaric')
        need_any_diags = need_any_diags .or. need_qc_isobaric
        need_qg_isobaric = MPAS_field_will_be_written('qg_isobaric')
        need_any_diags = need_any_diags .or. need_qg_isobaric
        need_qi_isobaric = MPAS_field_will_be_written('qi_isobaric')
        need_any_diags = need_any_diags .or. need_qi_isobaric

!        print *, "mp_thompson2=", need_mp_thompson
!        print *, "cu_grell_freitas2=", need_cu_grell_freitas
!        print *, "bl_mynn2=", need_bl_mynn

        if (need_any_diags) then
            call interp_diagnostics(mesh, state, 1, diag, diag_physics, tend_physics)
        end if
   
    end subroutine isobaric_diagnostics_compute


   !==================================================================================================
    subroutine interp_diagnostics(mesh, state, time_lev, diag, diag_physics, tend_physics)
   !==================================================================================================

       !input arguments:
        type (mpas_pool_type), intent(in) :: mesh
        type (mpas_pool_type), intent(in) :: state
        integer, intent(in) :: time_lev              ! which time level to use from state
       
       !inout arguments:
        type (mpas_pool_type), intent(inout) :: diag
        type (mpas_pool_type), intent(inout) :: diag_physics
        type (mpas_pool_type), intent(inout) :: tend_physics
       
       !local variables:

        integer :: nIsoLevelsU
        integer :: iCell,iVert,iVertD,k,kk
        integer, pointer :: nCells, nCellsSolve, nVertLevels, nVertices, vertexDegree, nIsoLevelsT, nIsoLevelsZ
        integer :: nVertLevelsP1
        integer, pointer :: index_qv, index_ni, index_nr, num_scalars, index_qr, index_qs, index_qc, index_qg, index_qi
        integer, dimension(:,:), pointer :: cellsOnVertex

        type (field2DReal), pointer:: pressure_p_field
       
        real (kind=RKIND), dimension(:), pointer :: areaTriangle, latCell
        real (kind=RKIND), dimension(:,:), pointer :: kiteAreasOnVertex

        real (kind=RKIND), dimension(:,:), pointer :: rho
        real (kind=RKIND), dimension(:,:), pointer :: exner, height
        real (kind=RKIND), dimension(:,:), pointer :: pressure_b, pressure_p 
        real (kind=RKIND), dimension(:,:), pointer :: relhum, theta_m, vorticity
        real (kind=RKIND), dimension(:,:), pointer :: umeridional, uzonal, vvel
        real (kind=RKIND), dimension(:,:), pointer :: rqvblten, rqvcuten, rthratenlw, rthratensw, rthcuten,rthblten
        real (kind=RKIND), dimension(:,:), pointer :: rucuten, rublten
        real (kind=RKIND), dimension(:,:), pointer :: rvcuten, rvblten
        real (kind=RKIND), dimension(:,:), pointer :: cldfrac
        real (kind=RKIND), dimension(:,:), pointer :: kzq, kzh, kzm
        real (kind=RKIND), dimension(:,:,:), pointer :: scalars
       
        real (kind=RKIND), dimension(:), pointer :: t_iso_levels
        real (kind=RKIND), dimension(:), pointer :: z_iso_levels
        real (kind=RKIND), dimension(:  ), pointer :: u_iso_levels
        real (kind=RKIND), dimension(:,:), pointer :: t_isobaric
        real (kind=RKIND), dimension(:,:), pointer :: z_isobaric
        real (kind=RKIND), dimension(:), pointer :: meanT_500_300
        real (kind=RKIND), dimension(:,:), pointer :: zgeo_isobaric

        real (kind=RKIND), dimension(:,:), pointer :: cldfrac_isobaric
        real (kind=RKIND), dimension(:,:), pointer :: temperature_isobaric
        real (kind=RKIND), dimension(:,:), pointer :: dewpoint_isobaric
        real (kind=RKIND), dimension(:,:), pointer :: relhum_isobaric
        real (kind=RKIND), dimension(:,:), pointer :: w_isobaric
        real (kind=RKIND), dimension(:,:), pointer :: uzonal_isobaric
        real (kind=RKIND), dimension(:,:), pointer :: umeridional_isobaric
        real (kind=RKIND), dimension(:,:), pointer :: vorticity_isobaric
        real (kind=RKIND), dimension(:,:), pointer :: height_isobaric
        real (kind=RKIND), dimension(:,:), pointer :: spechum_isobaric
        real (kind=RKIND), dimension(:,:), pointer :: omega_isobaric
        real (kind=RKIND), dimension(:,:), pointer :: rqvblten_isobaric
        real (kind=RKIND), dimension(:,:), pointer :: rthblten_isobaric
        real (kind=RKIND), dimension(:,:), pointer :: rqvcuten_isobaric
        real (kind=RKIND), dimension(:,:), pointer :: rthcuten_isobaric
        real (kind=RKIND), dimension(:,:), pointer :: rthratenlw_isobaric
        real (kind=RKIND), dimension(:,:), pointer :: rthratensw_isobaric
        real (kind=RKIND), dimension(:,:), pointer :: rucuten_isobaric
        real (kind=RKIND), dimension(:,:), pointer :: rvcuten_isobaric
        real (kind=RKIND), dimension(:,:), pointer :: rvblten_isobaric
        real (kind=RKIND), dimension(:,:), pointer :: rublten_isobaric
        real (kind=RKIND), dimension(:,:), pointer :: kzq_isobaric
        real (kind=RKIND), dimension(:,:), pointer :: kzh_isobaric
        real (kind=RKIND), dimension(:,:), pointer :: kzm_isobaric
        real (kind=RKIND), dimension(:,:), pointer :: ni_isobaric
        real (kind=RKIND), dimension(:,:), pointer :: nr_isobaric
        real (kind=RKIND), dimension(:,:), pointer :: qr_isobaric
        real (kind=RKIND), dimension(:,:), pointer :: qs_isobaric
        real (kind=RKIND), dimension(:,:), pointer :: qv_isobaric
        real (kind=RKIND), dimension(:,:), pointer :: qc_isobaric
        real (kind=RKIND), dimension(:,:), pointer :: qg_isobaric
        real (kind=RKIND), dimension(:,:), pointer :: qi_isobaric
         
       
        real (kind=RKIND) :: evp
       
       !--------------------
       
        real (kind=RKIND), dimension(:), pointer :: mslp
       
        real (kind=RKIND), dimension(:,:), allocatable :: pressure, pressureCp1, pressure2, pressure_v, temperature
        real (kind=RKIND), dimension(:,:), allocatable :: dewpoint, spechum
       
       !local interpolated fields:
        integer :: nIntP
        real (kind=RKIND) :: w1,w2,z0,z1,z2
        real (kind=RKIND), dimension(:,:), allocatable :: field_in,press_in
        real (kind=RKIND), dimension(:,:), allocatable :: field_interp,press_interp
        
        parameter (grav=9.80616) !LFR
        
       !--------------------------------------------------------------------------------------------------
       
       ! call mpas_log_write('')
       ! call mpas_log_write('--- enter subroutine interp_diagnostics:')
       
        call mpas_pool_get_dimension(mesh, 'nCells', nCells)
        call mpas_pool_get_dimension(mesh, 'nCellsSolve', nCellsSolve)
        call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)
        call mpas_pool_get_dimension(mesh, 'nVertices', nVertices)
        call mpas_pool_get_dimension(mesh, 'vertexDegree', vertexDegree)
        call mpas_pool_get_dimension(mesh, 'nIsoLevelsT', nIsoLevelsT)
        call mpas_pool_get_dimension(mesh, 'nIsoLevelsZ', nIsoLevelsZ)
        call mpas_pool_get_dimension(state, 'index_qv', index_qv)
        call mpas_pool_get_dimension(state, 'index_ni', index_ni)
        call mpas_pool_get_dimension(state, 'index_nr', index_nr)
        call mpas_pool_get_dimension(state, 'index_qr', index_qr)
        call mpas_pool_get_dimension(state, 'index_qs', index_qs)
        call mpas_pool_get_dimension(state, 'index_qc', index_qc)
        call mpas_pool_get_dimension(state, 'index_qg', index_qg)
        call mpas_pool_get_dimension(state, 'index_qi', index_qi)
        call mpas_pool_get_dimension(state, 'num_scalars', num_scalars)
       
        nVertLevelsP1 = nVertLevels + 1
       
        call mpas_pool_get_array(mesh, 'cellsOnVertex', cellsOnVertex)
        call mpas_pool_get_array(mesh, 'areaTriangle', areaTriangle)
        call mpas_pool_get_array(mesh, 'kiteAreasOnVertex', kiteAreasOnVertex)
       
        call mpas_pool_get_array(mesh, 'latCell', latCell)
        call mpas_pool_get_array(mesh, 'zgrid', height)
        call mpas_pool_get_array(state, 'w', vvel, time_lev)
        call mpas_pool_get_array(state, 'theta_m', theta_m, time_lev)
        call mpas_pool_get_array(state, 'scalars', scalars, time_lev)

        call mpas_pool_get_array(tend_physics, 'rqvcuten', rqvcuten)
        call mpas_pool_get_array(tend_physics, 'rthcuten', rthcuten)

        call mpas_pool_get_array(tend_physics, 'rqvblten', rqvblten)
        call mpas_pool_get_array(tend_physics, 'rthblten', rthblten)

        call mpas_pool_get_array(tend_physics, 'rthratenlw', rthratenlw)
        call mpas_pool_get_array(tend_physics, 'rthratensw', rthratensw)
        call mpas_pool_get_array(tend_physics, 'rucuten', rucuten)
        call mpas_pool_get_array(tend_physics, 'rvcuten', rvcuten)
        call mpas_pool_get_array(tend_physics, 'rvblten', rvblten)
        call mpas_pool_get_array(tend_physics, 'rublten', rublten)

        call mpas_pool_get_array(diag_physics, 'kzq', kzq)
        call mpas_pool_get_array(diag_physics, 'kzh', kzh)
        call mpas_pool_get_array(diag_physics, 'kzm', kzm)
        call mpas_pool_get_array(diag_physics,'cldfrac',cldfrac)


        call mpas_pool_get_field(diag, 'pressure_p', pressure_p_field)
        call mpas_dmpar_exch_halo_field(pressure_p_field)
       
        call mpas_pool_get_array(diag, 'exner', exner)
        call mpas_pool_get_array(diag, 'pressure_base', pressure_b)
        call mpas_pool_get_array(diag, 'pressure_p', pressure_p)
        call mpas_pool_get_array(diag, 'vorticity', vorticity)
        call mpas_pool_get_array(diag, 'uReconstructMeridional', umeridional)
        call mpas_pool_get_array(diag, 'uReconstructZonal', uzonal)
        call mpas_pool_get_array(diag, 'relhum', relhum)
        call mpas_pool_get_array(diag, 'rho', rho) !LFR
   
        call mpas_pool_get_array(diag, 'zgeo_isobaric', zgeo_isobaric)
    
        call mpas_pool_get_array(diag, 't_iso_levels', t_iso_levels)
        call mpas_pool_get_array(diag, 'u_iso_levels', u_iso_levels)
        call mpas_pool_get_array(diag, 'z_iso_levels', z_iso_levels)
        call mpas_pool_get_array(diag, 't_isobaric', t_isobaric)
        call mpas_pool_get_array(diag, 'z_isobaric', z_isobaric)
        call mpas_pool_get_array(diag, 'meanT_500_300', meanT_500_300)

        call mpas_pool_get_array(diag, 'cldfrac_isobaric', cldfrac_isobaric)
        call mpas_pool_get_array(diag, 'temperature_isobaric', temperature_isobaric)
        call mpas_pool_get_array(diag, 'dewpoint_isobaric', dewpoint_isobaric)
        call mpas_pool_get_array(diag, 'relhum_isobaric', relhum_isobaric)
        call mpas_pool_get_array(diag, 'w_isobaric', w_isobaric)
        call mpas_pool_get_array(diag, 'uzonal_isobaric', uzonal_isobaric)
        call mpas_pool_get_array(diag, 'umeridional_isobaric', umeridional_isobaric)
        call mpas_pool_get_array(diag, 'vorticity_isobaric', vorticity_isobaric)
        call mpas_pool_get_array(diag, 'height_isobaric', height_isobaric)
        call mpas_pool_get_array(diag, 'spechum_isobaric', spechum_isobaric)
        call mpas_pool_get_array(diag, 'omega_isobaric', omega_isobaric)
        call mpas_pool_get_array(diag, 'rqvblten_isobaric', rqvblten_isobaric)
        call mpas_pool_get_array(diag, 'rthblten_isobaric', rthblten_isobaric)
        call mpas_pool_get_array(diag, 'rqvcuten_isobaric', rqvcuten_isobaric)
        call mpas_pool_get_array(diag, 'rthcuten_isobaric', rthcuten_isobaric)
        call mpas_pool_get_array(diag, 'rthratenlw_isobaric', rthratenlw_isobaric)
        call mpas_pool_get_array(diag, 'rthratensw_isobaric', rthratensw_isobaric)
        call mpas_pool_get_array(diag, 'rucuten_isobaric', rucuten_isobaric)
        call mpas_pool_get_array(diag, 'rvcuten_isobaric', rvcuten_isobaric)
        call mpas_pool_get_array(diag, 'rvblten_isobaric', rvblten_isobaric)
        call mpas_pool_get_array(diag, 'rublten_isobaric', rublten_isobaric)
        call mpas_pool_get_array(diag, 'kzq_isobaric', kzq_isobaric)
        call mpas_pool_get_array(diag, 'kzh_isobaric', kzh_isobaric)
        call mpas_pool_get_array(diag, 'kzm_isobaric', kzm_isobaric)
        call mpas_pool_get_array(diag, 'ni_isobaric', ni_isobaric)
        call mpas_pool_get_array(diag, 'nr_isobaric', nr_isobaric)
        call mpas_pool_get_array(diag, 'qr_isobaric', qr_isobaric)
        call mpas_pool_get_array(diag, 'qs_isobaric', qs_isobaric)
        call mpas_pool_get_array(diag, 'qv_isobaric', qv_isobaric)
        call mpas_pool_get_array(diag, 'qc_isobaric', qc_isobaric)
        call mpas_pool_get_array(diag, 'qg_isobaric', qg_isobaric)
        call mpas_pool_get_array(diag, 'qi_isobaric', qi_isobaric)
 

       
        call mpas_pool_get_array(diag, 'mslp', mslp)
       
        if(.not.allocated(pressure)    ) allocate(pressure(nVertLevels,nCells)      )
        if(.not.allocated(pressureCp1) ) allocate(pressureCp1(nVertLevels,nCells+1) )
        if(.not.allocated(pressure2)   ) allocate(pressure2(nVertLevelsP1,nCells)   )
        if(.not.allocated(pressure_v)  ) allocate(pressure_v(nVertLevels,nVertices) )
        if(.not.allocated(temperature) ) allocate(temperature(nVertLevels,nCells)   )
        if(.not.allocated(spechum)  ) allocate(spechum(nVertLevels,nCells)    ) !LFR
        if(.not.allocated(dewpoint) ) allocate(dewpoint(nVertLevels,nCells)   ) 
       
        if (need_t_isobaric .or. need_zgeo_isobaric) then
                t_iso_levels(:) = (/   100000.0,   & !CR: surface 
                                       97500.0,    &
                                       95000.0,    &
                                       92500.0,    &
                                       90000.0,    &
                                       87500.0,    &
                                       85000.0,    &
                                       82500.0,    &
                                       70000.0,    &
                                       60000.0,    &
                                       50000.0,    &
                                       40000.0,    &
                                       30000.0,    &
                                       25000.0,    &
                                       20000.0,    &
                                       15000.0,    &
                                       10000.0,    &
                                       7000.0,     &
                                       5000.0,     &
                                       3000.0,     &
                                       2000.0,     &
                                       1500.0      /) !CR: top
              z_iso_levels = t_iso_levels
              u_iso_levels = t_iso_levels
        end if
       
        if (need_z_isobaric) then
           z_iso_levels = t_iso_levels
            nIsoLevelsZ = size(z_iso_levels)
        end if

        if (need_u_isobaric) then
            u_iso_levels = t_iso_levels
            nIsoLevelsU = size(u_iso_levels)
        end if
       
       !calculation of total pressure at cell centers (at mass points):
        do iCell = 1, nCells
        do k = 1, nVertLevels
           pressure(k,iCell)    = (pressure_p(k,iCell) + pressure_b(k,iCell)) / 100._RKIND
           pressureCp1(k,iCell) = pressure(k,iCell)
        enddo
        enddo
        do iCell = nCells+1, nCells+1
        do k = 1, nVertLevels
           pressureCp1(k,iCell)   = (pressure_p(k,iCell) + pressure_b(k,iCell)) / 100._RKIND
        enddo
        enddo
       
       !calculation of total pressure at cell centers (at vertical velocity points):
        k = nVertLevelsP1
        do iCell = 1, nCells
           z0 = height(k,iCell)
           z1 = 0.5*(height(k,iCell)+height(k-1,iCell)) 
           z2 = 0.5*(height(k-1,iCell)+height(k-2,iCell))
           w1 = (z0-z2)/(z1-z2)
           w2 = 1.-w1
           !use log of pressure to avoid occurrences of negative top-of-the-model pressure.
           pressure2(k,iCell) = exp(w1*log(pressure(k-1,iCell))+w2*log(pressure(k-2,iCell)))
        enddo
        do k = 2, nVertLevels
        do iCell = 1, nCells
           w1 = (height(k,iCell)-height(k-1,iCell)) / (height(k+1,iCell)-height(k-1,iCell))
           w2 = (height(k+1,iCell)-height(k,iCell)) / (height(k+1,iCell)-height(k-1,iCell))
           ! pressure2(k,iCell) = w1*pressure(k,iCell) + w2*pressure(k-1,iCell)
           !
           ! switch to use ln(pressure) for more accurate vertical interpolation, WCS 20230407
           pressure2(k,iCell) = exp(w1*log(pressure(k,iCell))+w2*log(pressure(k-1,iCell)))
        enddo
        enddo
        k = 1
        do iCell = 1, nCells
           z0 = height(k,iCell)
           z1 = 0.5*(height(k,iCell)+height(k+1,iCell)) 
           z2 = 0.5*(height(k+1,iCell)+height(k+2,iCell))
           w1 = (z0-z2)/(z1-z2)
           w2 = 1.-w1
           ! pressure2(k,iCell) = w1*pressure(k,iCell)+w2*pressure(k+1,iCell)
           !
           ! switch to use ln(pressure) for more accurate vertical interpolation, WCS 20230407
           pressure2(k,iCell) = exp(w1*log(pressure(k,iCell))+w2*log(pressure(k+1,iCell)))
        enddo
       
       !calculation of total pressure at cell vertices (at mass points):
        do iVert = 1, nVertices
           pressure_v(:,iVert) = 0._RKIND
       
           do k = 1, nVertLevels
           do iVertD = 1, vertexDegree
              pressure_v(k,iVert) = pressure_v(k,iVert) &
                      + kiteAreasOnVertex(iVertD,iVert)*pressureCp1(k,cellsOnVertex(iVertD,iVert))
           enddo
           pressure_v(k,iVert) = pressure_v(k,iVert) / areaTriangle(iVert)
           enddo
        enddo
       
        if (NEED_TEMPERATURE_ISOBARIC .or. NEED_RELHUM_ISOBARIC .or. NEED_SPECHUM_ISOBARIC &
            .or. NEED_DEWPOINT_ISOBARIC .or. need_mslp) then

           !calculation of temperature at cell centers:
            do iCell = 1,nCells
            do k = 1,nVertLevels
                temperature(k,iCell) = (theta_m(k,iCell)/(1._RKIND+rvord*scalars(index_qv,k,iCell)))*exner(k,iCell) 

                ! Vapor pressure (NB: pressure here is already in hPa)
                evp = pressure(k,iCell) * scalars(index_qv,k,iCell) / (scalars(index_qv,k,iCell) + 0.622_RKIND)
                evp = max(evp, 1.0e-8_RKIND)

                spechum(k,iCell) = 0.622*evp/(pressure2(k,iCell)-0.378*evp) !LFR

                ! Dewpoint temperature following Bolton (1980)
                dewpoint(k,iCell) = (243.5_RKIND * log(evp/6.112_RKIND)) / (17.67_RKIND - log(evp/6.112_RKIND))
                dewpoint(k,iCell) = dewpoint(k,iCell) + 273.15
            enddo
            enddo
        end if
       
       !interpolation to fixed pressure levels for fields located at cells centers and at mass points:
        !CR: nIntP = 22 this var should not be used any more! use nIsoLevelsT (came from registry xml)
        
        !CR: invertendo field_interp de (nCells,nIntP) para (nIsoLevelsT, nCells), to work in new vertical 
        !interpolation subroutine
        
        if(.not.allocated(field_interp)) allocate(field_interp(nIsoLevelsT, nCells) )
        if(.not.allocated(press_interp)) allocate(press_interp(nCells,nIsoLevelsT) ) !CR: not used!
        do iCell = 1, nCells
           press_interp(iCell, 1) = 15.0_RKIND
           press_interp(iCell, 2) = 20.0_RKIND
           press_interp(iCell, 3) = 30.0_RKIND
           press_interp(iCell, 4) = 50.0_RKIND
           press_interp(iCell, 5) = 70.0_RKIND
           press_interp(iCell, 6) = 100.0_RKIND
           press_interp(iCell, 7) = 150.0_RKIND
           press_interp(iCell, 8) = 200.0_RKIND
           press_interp(iCell, 9) = 250.0_RKIND
           press_interp(iCell,10) = 300.0_RKIND
           press_interp(iCell,11) = 400.0_RKIND
           press_interp(iCell,12) = 500.0_RKIND
           press_interp(iCell,13) = 600.0_RKIND
           press_interp(iCell,14) = 700.0_RKIND
           press_interp(iCell,15) = 825.0_RKIND
           press_interp(iCell,16) = 850.0_RKIND
           press_interp(iCell,17) = 875.0_RKIND
           press_interp(iCell,18) = 900.0_RKIND
           press_interp(iCell,19) = 925.0_RKIND
           press_interp(iCell,20) = 950.0_RKIND
           press_interp(iCell,21) = 975.0_RKIND
           press_interp(iCell,22) = 1000.0_RKIND
        enddo
       
       
        if(.not.allocated(field_in)) allocate(field_in(nCells,nVertLevels)) !CR: not used any more

        if (NEED_CLDFRAC_ISOBARIC) then
         !... cldfrac:
         
         !CR: 1o single
         !CR: before: call interp_tofixed_pressure(nCells,nVertLevels,nIntP,press_in,field_in,press_interp,field_interp)
            
         call VerticalInterpolation_Single ( nCells                                    , & !1
                                             nVertLevels                               , & !2
                                             nVertLevels                               , & !3 CR:duvida: ta certo?
                                             nIsoLevelsT                               , & !4 CR: ta certo?
                                             nIsoLevelsT                               , & !5 CR: ta certo?
                                             cldfrac     (1:nVertLevels,1:nCells)      , & !6
                                             field_interp(1:nIsoLevelsT,1:nCells)      , & !7
                                             pressure2   (1:nVertLevels,1:nCells)*100.0, & !8
                                             t_iso_levels(1:nIsoLevelsT)               )   !9
                                                
            cldfrac_isobaric(1:nIsoLevelsT,1:nCells) =  field_interp(1:nIsoLevelsT,1:nCells) 

            ! call mpas_log_write('--- end interpolate cldfrac:')
        end if
        !CR: new vertical interpolation done! 1


        if (NEED_TEMPERATURE_ISOBARIC) then
           !... temperature:
           !CR: 2o thermal2
           !CR: call interp_tofixed_pressure(nCells,nVertLevels,nIntP,press_in,field_in,press_interp,field_interp)

           call VerticalInterpolation_Thermal2 ( nCells                                    , & !1
                                                 nVertLevels                               , & !2
                                                 nIsoLevelsT                               , & !3
                                                 pressure2   (1:nVertLevels,1:nCells)*100.0, & !4
                                                 temperature (1:nVertLevels,1:nCells)      , & !5
                                                 field_interp(1:nIsoLevelsT,1:nCells)      , & !6
                                                 t_iso_levels(1:nIsoLevelsT)               )   !7

            temperature_isobaric(1:nIsoLevelsT,1:nCells) =  field_interp(1:nIsoLevelsT,1:nCells)

            ! call mpas_log_write('--- end interpolate temperature:')
        end if
        !CR: new vertical interpolation done! 2


        if (need_mp_thompson .eqv. .true. .or. need_bl_mynn .eqv. .true.) then
           if (NEED_NI_ISOBARIC) then
               !... ni:
               !CR: 3o single
               !CR: call interp_tofixed_pressure(nCells,nVertLevels,nIntP,press_in,field_in,press_interp,field_interp)
                
               call VerticalInterpolation_Single ( nCells                                        , & !1
                                                   nVertLevels                                   , & !2
                                                   nVertLevels                                   , & !3
                                                   nIsoLevelsT                                   , & !4
                                                   nIsoLevelsT                                   , & !5
                                                   scalars(index_ni,1:nVertLevels,1:nCells)      , & !6
                                                   field_interp    (1:nIsoLevelsT,1:nCells)      , & !7
                                                   pressure2       (1:nVertLevels,1:nCells)*100.0, & !8
                                                   t_iso_levels    (1:nIsoLevelsT)               )   !9
                

               ni_isobaric(1:nIsoLevelsT,1:nCells) =  field_interp(1:nIsoLevelsT,1:nCells)
 
               ! call mpas_log_write('--- end interpolate ni:')
            end if
        end if
        !CR: new vertical interpolation done! 3

        if (need_mp_thompson .eqv. .true.) then
            if (NEED_NR_ISOBARIC) then
               !... nr:
               !CR: 4o single
               !CR: call interp_tofixed_pressure(nCells,nVertLevels,nIntP,press_in,field_in,press_interp,field_interp)
                
               call VerticalInterpolation_Single ( nCells                                        , & !1
                                                   nVertLevels                                   , & !2
                                                   nVertLevels                                   , & !3
                                                   nIsoLevelsT                                   , & !4
                                                   nIsoLevelsT                                   , & !5
                                                   scalars(index_nr,1:nVertLevels,1:nCells)      , & !6
                                                   field_interp    (1:nIsoLevelsT,1:nCells)      , & !7
                                                   pressure2       (1:nVertLevels,1:nCells)*100.0, & !8
                                                   t_iso_levels    (1:nIsoLevelsT)               )   !9


               nr_isobaric(1:nIsoLevelsT,1:nCells) =  field_interp(1:nIsoLevelsT,1:nCells)
 
               ! call mpas_log_write('--- end interpolate nr:')
            end if
        end if
        !CR: new vertical interpolation done! 4

         if (NEED_QR_ISOBARIC) then
            !... qr:
            !CR: 5o single
            !CR: call interp_tofixed_pressure(nCells,nVertLevels,nIntP,press_in,field_in,press_interp,field_interp)

            call VerticalInterpolation_Single ( nCells                                        , & !1
                                                nVertLevels                                   , & !2
                                                nVertLevels                                   , & !3
                                                nIsoLevelsT                                   , & !4
                                                nIsoLevelsT                                   , & !5
                                                scalars(index_qr,1:nVertLevels,1:nCells)      , & !6
                                                field_interp    (1:nIsoLevelsT,1:nCells)      , & !7
                                                pressure2       (1:nVertLevels,1:nCells)*100.0, & !8
                                                t_iso_levels    (1:nIsoLevelsT)               )   !9
            

            qr_isobaric(1:nIsoLevelsT,1:nCells) =  field_interp(1:nIsoLevelsT,1:nCells)
 
            ! call mpas_log_write('--- end interpolate qr:')
        end if
        !CR: new vertical interpolation done! 5

         if (NEED_QS_ISOBARIC) then
            !... qs:
            !CR: 6o single
            !CR: call interp_tofixed_pressure(nCells,nVertLevels,nIntP,press_in,field_in,press_interp,field_interp)

            call VerticalInterpolation_Single ( nCells                                        , & !1
                                                nVertLevels                                   , & !2
                                                nVertLevels                                   , & !3
                                                nIsoLevelsT                                   , & !4
                                                nIsoLevelsT                                   , & !5
                                                scalars(index_qs,1:nVertLevels,1:nCells)      , & !6
                                                field_interp    (1:nIsoLevelsT,1:nCells)      , & !7
                                                pressure2       (1:nVertLevels,1:nCells)*100.0, & !8
                                                t_iso_levels    (1:nIsoLevelsT)               )   !9

            qs_isobaric(1:nIsoLevelsT,1:nCells) =  field_interp(1:nIsoLevelsT,1:nCells)
 
            ! call mpas_log_write('--- end interpolate qs:')
        end if
        !CR: new vertical interpolation done! 6

        if (NEED_QV_ISOBARIC) then
           !... qv:
            !CR: 7o single
            !CR: call interp_tofixed_pressure(nCells,nVertLevels,nIntP,press_in,field_in,press_interp,field_interp)

            call VerticalInterpolation_Single ( nCells                                        , & !1
                                                nVertLevels                                   , & !2
                                                nVertLevels                                   , & !3
                                                nIsoLevelsT                                   , & !4
                                                nIsoLevelsT                                   , & !5
                                                scalars(index_qv,1:nVertLevels,1:nCells)      , & !6
                                                field_interp    (1:nIsoLevelsT,1:nCells)      , & !7
                                                pressure2       (1:nVertLevels,1:nCells)*100.0, & !8
                                                t_iso_levels    (1:nIsoLevelsT)               )   !9

            qv_isobaric(1:nIsoLevelsT,1:nCells) =  field_interp(1:nIsoLevelsT,1:nCells)
 
            ! call mpas_log_write('--- end interpolate qv:')
        end if
        !CR: new vertical interpolation done! 7

        if (NEED_QC_ISOBARIC) then
           !... qc:
            !CR: 8o single
            !CR: call interp_tofixed_pressure(nCells,nVertLevels,nIntP,press_in,field_in,press_interp,field_interp)

            call VerticalInterpolation_Single ( nCells                                        , & !1
                                                nVertLevels                                   , & !2
                                                nVertLevels                                   , & !3
                                                nIsoLevelsT                                   , & !4
                                                nIsoLevelsT                                   , & !5
                                                scalars(index_qc,1:nVertLevels,1:nCells)      , & !6
                                                field_interp    (1:nIsoLevelsT,1:nCells)      , & !7
                                                pressure2       (1:nVertLevels,1:nCells)*100.0, & !8
                                                t_iso_levels    (1:nIsoLevelsT)               )   !9

            qc_isobaric(1:nIsoLevelsT,1:nCells) =  field_interp(1:nIsoLevelsT,1:nCells)

            ! call mpas_log_write('--- end interpolate qc:')
        end if
        !CR: new vertical interpolation done! 8

        if (NEED_QG_ISOBARIC) then
           !... qg:
            !CR: 9o single
            !CR: call interp_tofixed_pressure(nCells,nVertLevels,nIntP,press_in,field_in,press_interp,field_interp)

            call VerticalInterpolation_Single ( nCells                                        , & !1
                                                nVertLevels                                   , & !2
                                                nVertLevels                                   , & !3
                                                nIsoLevelsT                                   , & !4
                                                nIsoLevelsT                                   , & !5
                                                scalars(index_qg,1:nVertLevels,1:nCells)      , & !6
                                                field_interp    (1:nIsoLevelsT,1:nCells)      , & !7
                                                pressure2       (1:nVertLevels,1:nCells)*100.0, & !8
                                                t_iso_levels    (1:nIsoLevelsT)               )   !9

            qg_isobaric(1:nIsoLevelsT,1:nCells) =  field_interp(1:nIsoLevelsT,1:nCells)
 
            ! call mpas_log_write('--- end interpolate qg:')
        end if
        !CR: new vertical interpolation done! 9

        if (NEED_QI_ISOBARIC) then
           !... qi:
            !CR: 10o single
            !CR: call interp_tofixed_pressure(nCells,nVertLevels,nIntP,press_in,field_in,press_interp,field_interp)

            call VerticalInterpolation_Single ( nCells                                        , & !1
                                                nVertLevels                                   , & !2
                                                nVertLevels                                   , & !3
                                                nIsoLevelsT                                   , & !4
                                                nIsoLevelsT                                   , & !5
                                                scalars(index_qi,1:nVertLevels,1:nCells)      , & !6
                                                field_interp    (1:nIsoLevelsT,1:nCells)      , & !7
                                                pressure2       (1:nVertLevels,1:nCells)*100.0, & !8
                                                t_iso_levels    (1:nIsoLevelsT)               )   !9

            qi_isobaric(1:nIsoLevelsT,1:nCells) =  field_interp(1:nIsoLevelsT,1:nCells)

            ! call mpas_log_write('--- end interpolate qi:')
        end if
        !CR: new vertical interpolation done! 10

        if (NEED_RQVBLTEN_ISOBARIC) then
           !... rqvblten:
            !CR: 11o single
            !CR: call interp_tofixed_pressure(nCells,nVertLevels,nIntP,press_in,field_in,press_interp,field_interp)

            call VerticalInterpolation_Single ( nCells                                    , & !1
                                                nVertLevels                               , & !2
                                                nVertLevels                               , & !3
                                                nIsoLevelsT                               , & !4
                                                nIsoLevelsT                               , & !5
                                                rqvblten    (1:nVertLevels,1:nCells)      , & !6
                                                field_interp(1:nIsoLevelsT,1:nCells)      , & !7
                                                pressure2   (1:nVertLevels,1:nCells)*100.0, & !8
                                                t_iso_levels(1:nIsoLevelsT)               )   !9

            rqvblten_isobaric(1:nIsoLevelsT,1:nCells) =  field_interp(1:nIsoLevelsT,1:nCells)
 
            ! call mpas_log_write('--- end interpolate rqvblten:')
        end if
        !CR: new vertical interpolation done! 11
        

        if (NEED_RTHBLTEN_ISOBARIC) then
           !... rthblten:
           !EGK: updating interpolation method
           ! do iCell = 1, nCells
           ! do k = 1, nVertLevels
           !    kk = nVertLevels+1-k
           !    field_in(iCell,kk) = rthblten(k,iCell)
           ! enddo
           ! enddo
           ! call interp_tofixed_pressure(nCells,nVertLevels,nIntP,press_in,field_in,press_interp,field_interp)

            call VerticalInterpolation_Single ( nCells                                    , & !1
                                                nVertLevels                               , & !2
                                                nVertLevels                               , & !3
                                                nIsoLevelsT                               , & !4
                                                nIsoLevelsT                               , & !5
                                                rthblten    (1:nVertLevels,1:nCells)      , & !6
                                                field_interp(1:nIsoLevelsT,1:nCells)      , & !7
                                                pressure2   (1:nVertLevels,1:nCells)*100.0, & !8
                                                t_iso_levels(1:nIsoLevelsT)               )   !9

            rthblten_isobaric(1:nIsoLevelsT,1:nCells) =  field_interp(1:nIsoLevelsT,1:nCells)

            ! call mpas_log_write('--- end interpolate rthblten:')
        end if
        !EGK: new vertical interpolation done!


        if (NEED_RQVCUTEN_ISOBARIC) then
           !... rqvcuten:
            !CR: 12o single
            !CR: call interp_tofixed_pressure(nCells,nVertLevels,nIntP,press_in,field_in,press_interp,field_interp)

            call VerticalInterpolation_Single ( nCells                                    , & !1
                                                nVertLevels                               , & !2
                                                nVertLevels                               , & !3
                                                nIsoLevelsT                               , & !4
                                                nIsoLevelsT                               , & !5
                                                rqvcuten    (1:nVertLevels,1:nCells)      , & !6
                                                field_interp(1:nIsoLevelsT,1:nCells)      , & !7
                                                pressure2   (1:nVertLevels,1:nCells)*100.0, & !8
                                                t_iso_levels(1:nIsoLevelsT)               )   !9

            rqvcuten_isobaric(1:nIsoLevelsT,1:nCells) =  field_interp(1:nIsoLevelsT,1:nCells)

            ! call mpas_log_write('--- end interpolate rqvcuten:')
        end if
        !CR: new vertical interpolation done! 12


        if (NEED_RTHCUTEN_ISOBARIC) then
           !... rthcuten:
           !EGK: updating interpolation method
           ! do iCell = 1, nCells
           ! do k = 1, nVertLevels
           !    kk = nVertLevels+1-k
           !   field_in(iCell,kk) = rthcuten(k,iCell)
           ! enddo
           ! enddo
           ! call interp_tofixed_pressure(nCells,nVertLevels,nIntP,press_in,field_in,press_interp,field_interp)

            call VerticalInterpolation_Single ( nCells                                    , & !1
                                                nVertLevels                               , & !2
                                                nVertLevels                               , & !3
                                                nIsoLevelsT                               , & !4
                                                nIsoLevelsT                               , & !5
                                                rthcuten    (1:nVertLevels,1:nCells)      , & !6
                                                field_interp(1:nIsoLevelsT,1:nCells)      , & !7
                                                pressure2   (1:nVertLevels,1:nCells)*100.0, & !8
                                                t_iso_levels(1:nIsoLevelsT)               )   !9

            rthcuten_isobaric(1:nIsoLevelsT,1:nCells) =  field_interp(1:nIsoLevelsT,1:nCells)

            ! call mpas_log_write('--- end interpolate rthcuten:')
        end if
        !EGK: new vertical interpolation done!


        if (NEED_RTHRATENLW_ISOBARIC) then
           !... rthratenlw:
            !CR: 13o single
            !CR: call interp_tofixed_pressure(nCells,nVertLevels,nIntP,press_in,field_in,press_interp,field_interp)

            call VerticalInterpolation_Single ( nCells                                    , & !1
                                                nVertLevels                               , & !2
                                                nVertLevels                               , & !3
                                                nIsoLevelsT                               , & !4
                                                nIsoLevelsT                               , & !5
                                                rthratenlw  (1:nVertLevels,1:nCells)      , & !6
                                                field_interp(1:nIsoLevelsT,1:nCells)      , & !7
                                                pressure2   (1:nVertLevels,1:nCells)*100.0, & !8
                                                t_iso_levels(1:nIsoLevelsT)               )   !9

            rthratenlw_isobaric(1:nIsoLevelsT,1:nCells) =  field_interp(1:nIsoLevelsT,1:nCells)
 
            ! call mpas_log_write('--- end interpolate rthratenlw:')
        end if
        !CR: new vertical interpolation done! 13

        if (NEED_RTHRATENSW_ISOBARIC) then
           !... rthratenlw:
            !CR: 14o single
            !CR: call interp_tofixed_pressure(nCells,nVertLevels,nIntP,press_in,field_in,press_interp,field_interp)

            call VerticalInterpolation_Single ( nCells                                    , & !1
                                                nVertLevels                               , & !2
                                                nVertLevels                               , & !3
                                                nIsoLevelsT                               , & !4
                                                nIsoLevelsT                               , & !5
                                                rthratensw  (1:nVertLevels,1:nCells)      , & !6
                                                field_interp(1:nIsoLevelsT,1:nCells)      , & !7
                                                pressure2   (1:nVertLevels,1:nCells)*100.0, & !8
                                                t_iso_levels(1:nIsoLevelsT)               ) !9

            rthratensw_isobaric(1:nIsoLevelsT,1:nCells) =  field_interp(1:nIsoLevelsT,1:nCells)
 
            ! call mpas_log_write('--- end interpolate rthratenlw:')
        end if
        !CR: new vertical interpolation done! 14



         if (NEED_RUCUTEN_ISOBARIC) then
            !... rucuten:
            !CR: 15o single
            !CR: call interp_tofixed_pressure(nCells,nVertLevels,nIntP,press_in,field_in,press_interp,field_interp)

            call VerticalInterpolation_Single ( nCells                                    , & !1
                                                nVertLevels                               , & !2
                                                nVertLevels                               , & !3
                                                nIsoLevelsT                               , & !4
                                                nIsoLevelsT                               , & !5
                                                rucuten     (1:nVertLevels,1:nCells)      , & !6
                                                field_interp(1:nIsoLevelsT,1:nCells)      , & !7
                                                pressure2   (1:nVertLevels,1:nCells)*100.0, & !8
                                                t_iso_levels(1:nIsoLevelsT)               ) !9

            rucuten_isobaric(1:nIsoLevelsT,1:nCells) =  field_interp(1:nIsoLevelsT,1:nCells)
 
            ! call mpas_log_write('--- end interpolate rucuten:')
        end if
        !CR: new vertical interpolation done! 15

         if (NEED_RVCUTEN_ISOBARIC) then
            !... rvcuten:
            !CR: 16o single
            !CR: call interp_tofixed_pressure(nCells,nVertLevels,nIntP,press_in,field_in,press_interp,field_interp)

            call VerticalInterpolation_Single ( nCells                                    , & !1
                                                nVertLevels                               , & !2
                                                nVertLevels                               , & !3
                                                nIsoLevelsT                               , & !4
                                                nIsoLevelsT                               , & !5
                                                rvcuten     (1:nVertLevels,1:nCells)      , & !6
                                                field_interp(1:nIsoLevelsT,1:nCells)      , & !7
                                                pressure2   (1:nVertLevels,1:nCells)*100.0, & !8  
                                                t_iso_levels(1:nIsoLevelsT)               ) !9

            rvcuten_isobaric(1:nIsoLevelsT,1:nCells) =  field_interp(1:nIsoLevelsT,1:nCells)
 
          ! call mpas_log_write('--- end interpolate rvcuten:')
        end if
        !CR: new vertical interpolation done! 16

        if (NEED_RUBLTEN_ISOBARIC) then
            !... rublten:
            !CR: 17o single
            !CR: call interp_tofixed_pressure(nCells,nVertLevels,nIntP,press_in,field_in,press_interp,field_interp)


            call VerticalInterpolation_Single ( nCells                                    , & !1
                                                nVertLevels                               , & !2
                                                nVertLevels                               , & !3
                                                nIsoLevelsT                               , & !4
                                                nIsoLevelsT                               , & !5
                                                rublten     (1:nVertLevels,1:nCells)      , & !6
                                                field_interp(1:nIsoLevelsT,1:nCells)      , & !7
                                                pressure2   (1:nVertLevels,1:nCells)*100.0, & !8
                                                t_iso_levels(1:nIsoLevelsT)               )   !9

            rublten_isobaric(1:nIsoLevelsT,1:nCells) =  field_interp(1:nIsoLevelsT,1:nCells)
 
          ! call mpas_log_write('--- end interpolate rublten:')
        end if
        !CR: new vertical interpolation done! 17

     
        if (NEED_RVBLTEN_ISOBARIC) then
            !... rvblten:
            !CR: 18o single
            !CR: call interp_tofixed_pressure(nCells,nVertLevels,nIntP,press_in,field_in,press_interp,field_interp)

            call VerticalInterpolation_Single ( nCells                                    , & !1
                                                nVertLevels                               , & !2
                                                nVertLevels                               , & !3
                                                nIsoLevelsT                               , & !4
                                                nIsoLevelsT                               , & !5
                                                rvblten     (1:nVertLevels,1:nCells)      , & !6
                                                field_interp(1:nIsoLevelsT,1:nCells)      , & !7
                                                pressure2   (1:nVertLevels,1:nCells)*100.0, & !8
                                                t_iso_levels(1:nIsoLevelsT)               )   !9

            rvblten_isobaric(1:nIsoLevelsT,1:nCells) =  field_interp(1:nIsoLevelsT,1:nCells)
 
          ! call mpas_log_write('--- end interpolate rvblten:')
        end if
        !CR: new vertical interpolation done! 18

        if (NEED_RELHUM_ISOBARIC) then
           !... relative humidity:
            !CR: 19o Humid
            !CR: call interp_tofixed_pressure(nCells,nVertLevels,nIntP,press_in,field_in,press_interp,field_interp)
            
            call VerticalInterpolation_Humid ( nCells                                    , & !1
                                               nVertLevels                               , & !2
                                               nIsoLevelsT                               , & !3
                                               pressure2   (1:nVertLevels,1:nCells)*100.0, & !4
                                               relhum      (1:nVertLevels,1:nCells)      , & !5
                                               field_interp(1:nIsoLevelsT,1:nCells)      , & !6
                                               t_iso_levels(1:nIsoLevelsT)               ) !7
                                             

            relhum_isobaric(1:nIsoLevelsT,1:nCells) =  field_interp(1:nIsoLevelsT,1:nCells)             
            
           ! call mpas_log_write('--- end interpolate relative humidity:')
        end if
        !CR: new vertical interpolation done! 19
       
        if (NEED_KZQ_ISOBARIC) then
           !... kzq:
            !CR: 20o single
            !CR: call interp_tofixed_pressure(nCells,nVertLevels,nIntP,press_in,field_in,press_interp,field_interp)

            call VerticalInterpolation_Single ( nCells                                    , & !1
                                                nVertLevels                               , & !2
                                                nVertLevels                               , & !3
                                                nIsoLevelsT                               , & !4
                                                nIsoLevelsT                               , & !5
                                                kzq         (1:nVertLevels,1:nCells)      , & !6
                                                field_interp(1:nIsoLevelsT,1:nCells)      , & !7
                                                pressure2   (1:nVertLevels,1:nCells)*100.0, & !8
                                                t_iso_levels(1:nIsoLevelsT)               )   !9

            kzq_isobaric(1:nIsoLevelsT,1:nCells) =  field_interp(1:nIsoLevelsT,1:nCells)
 
            ! call mpas_log_write('--- end interpolate kzq:')
        end if
        !CR: new vertical interpolation done! 20

        if (NEED_KZH_ISOBARIC) then
           !... kzh:
            !CR: 21o single
            !CR: call interp_tofixed_pressure(nCells,nVertLevels,nIntP,press_in,field_in,press_interp,field_interp)

            call VerticalInterpolation_Single ( nCells                                    , & !1
                                                nVertLevels                               , & !2
                                                nVertLevels                               , & !3
                                                nIsoLevelsT                               , & !4
                                                nIsoLevelsT                               , & !5
                                                kzh         (1:nVertLevels,1:nCells)      , & !6
                                                field_interp(1:nIsoLevelsT,1:nCells)      , & !7
                                                pressure2   (1:nVertLevels,1:nCells)*100.0, & !8
                                                t_iso_levels(1:nIsoLevelsT)               )   !9

            kzh_isobaric(1:nIsoLevelsT,1:nCells) =  field_interp(1:nIsoLevelsT,1:nCells)
 
            ! call mpas_log_write('--- end interpolate kzh:')
        end if
        !CR: new vertical interpolation done! 21

        if (NEED_KZM_ISOBARIC) then
           !... kzm:
            !CR: 22o single
            !CR: call interp_tofixed_pressure(nCells,nVertLevels,nIntP,press_in,field_in,press_interp,field_interp)

            call VerticalInterpolation_Single ( nCells                                    , & !1
                                                nVertLevels                               , & !2
                                                nVertLevels                               , & !3
                                                nIsoLevelsT                               , & !4
                                                nIsoLevelsT                               , & !5
                                                kzm         (1:nVertLevels,1:nCells)      , & !6
                                                field_interp(1:nIsoLevelsT,1:nCells)      , & !7
                                                pressure2   (1:nVertLevels,1:nCells)*100.0, & !8
                                                t_iso_levels(1:nIsoLevelsT)               )   !9

            kzm_isobaric(1:nIsoLevelsT,1:nCells) =  field_interp(1:nIsoLevelsT,1:nCells)
 
            ! call mpas_log_write('--- end interpolate kzm:')
        end if
        !CR: new vertical interpolation done! 22


        if (NEED_DEWPOINT_ISOBARIC) then
           !... dewpoint
            !CR: 23o Thermal2
            !CR: call interp_tofixed_pressure(nCells,nVertLevels,nIntP,press_in,field_in,press_interp,field_interp)
            !CR: dewpoint was in 8 levels = now its 22, like the others!
            
            call VerticalInterpolation_Thermal2 ( nCells                                    , & !1
                                                  nVertLevels                               , & !2
                                                  nIsoLevelsT                               , & !3
                                                  pressure2   (1:nVertLevels,1:nCells)*100.0, & !4
                                                  dewpoint    (1:nVertLevels,1:nCells)      , & !5
                                                  field_interp(1:nIsoLevelsT,1:nCells)      , & !6
                                                  t_iso_levels(1:nIsoLevelsT)               )   !7

            dewpoint_isobaric(1:nIsoLevelsT,1:nCells) =  field_interp(1:nIsoLevelsT,1:nCells)
            
           ! call mpas_log_write('--- end interpolate relative humidity:')
        end if
        !CR: new vertical interpolation done! 23
        
       
        if (NEED_UZONAL_ISOBARIC) then
           !... u zonal wind:
           !CR: 24o Wind 
            !CR: oldone: call interp_tofixed_pressure(1.nCells, 2.nVertLevels, 3.nIntP, 4.press_in, 5.field_in, 6.press_interp, 7.field_interp)

            call VerticalInterpolation_Wind ( nCells                                    , & !1. 
                                              nVertLevels                               , & !2. 
                                              nIsoLevelsT                               , & !3. 
                                              pressure2    (1:nVertLevels,1:nCells)*100.0, & !4. 
                                              t_iso_levels(1:nIsoLevelsT)               , & !5. 
                                              uzonal      (1:nVertLevels, 1:nCells)     , & !6. 
                                              field_interp(1:nIsoLevelsT,1:nCells)      )   !7. 
            

            uzonal_isobaric(1:nIsoLevelsT,1:nCells) =  field_interp(1:nIsoLevelsT,1:nCells)


                       ! call mpas_log_write('--- end interpolate zonal wind:')
        end if
        !CR: new vertical interpolation done! 24
       
        if (NEED_UMERIDIONAL_ISOBARIC) then
           !... u meridional wind:
           !CR: 25o Wind
            !CR: call interp_tofixed_pressure(nCells,nVertLevels,nIntP,press_in,field_in,press_interp,field_interp)

            call VerticalInterpolation_Wind (nCells                                    , & !1. 
                                             nVertLevels                               , & !2. 
                                             nIsoLevelsT                               , & !3. 
                                             pressure2    (1:nVertLevels,1:nCells)*100.0, & !4. 
                                             t_iso_levels(1:nIsoLevelsT)               , & !5. 
                                             umeridional (1:nVertLevels,1:nCells)      , & !6. 
                                             field_interp(1:nIsoLevelsT,1:nCells)      )   !7. 


            umeridional_isobaric(1:nIsoLevelsT,1:nCells) =  field_interp(1:nIsoLevelsT,1:nCells)

           ! call mpas_log_write('--- end interpolate meridional wind:')
        end if
        !CR: new vertical interpolation done! 25
       
        if (NEED_SPECHUM_ISOBARIC) then
           !... spechum:
           !CR: 26o Humid
           !CR: call interp_tofixed_pressure(nCells,nVertLevels,nIntP,press_in,field_in,press_interp,field_interp)

            call VerticalInterpolation_Humid(nCells                                    , & !1
                                             nVertLevels                               , & !2
                                             nIsoLevelsT                               , & !3
                                             pressure2    (1:nVertLevels,1:nCells)*100.0, & !4
                                             spechum     (1:nVertLevels,1:nCells)      , & !5
                                             field_interp(1:nIsoLevelsT,1:nCells)      , & !6
                                             t_iso_levels(1:nIsoLevelsT)               )   !7

            spechum_isobaric(1:nIsoLevelsT,1:nCells) =  field_interp(1:nIsoLevelsT,1:nCells)
 
           ! call mpas_log_write('--- end interpolate spechum:')
        end if
        !CR: new vertical interpolation done! 26



        if(allocated(field_in)) deallocate(field_in)
        if(allocated(press_in)) deallocate(press_in)
       
        if (NEED_W_ISOBARIC .or. NEED_HEIGHT_ISOBARIC .or. NEED_OMEGA_ISOBARIC &
            .or. NEED_ZGEO_ISOBARIC) then
           !interpolation to fixed pressure levels for fields located at cells centers and at vertical
           !velocity points:
            if(.not.allocated(press_in)) allocate(press_in(nCells,nVertLevelsP1))
       
           !CR: 27o Single
            !CR: call interp_tofixed_pressure(nCells,nVertLevelsP1,nIntP,press_in,field_in,press_interp,field_interp)

            call VerticalInterpolation_Single ( nCells                                    , & !1
                                                nVertLevels                               , & !2
                                                nVertLevels                               , & !3
                                                nIsoLevelsT                               , & !4
                                                nIsoLevelsT                               , & !5
                                                height      (1:nVertLevels,1:nCells)      , & !6
                                                field_interp(1:nIsoLevelsT,1:nCells)      , & !7
                                                pressure2    (1:nVertLevels,1:nCells)*100.0, & !8
                                                t_iso_levels(1:nIsoLevelsT)               )   !9

            height_isobaric(1:nIsoLevelsT,1:nCells) =  field_interp(1:nIsoLevelsT,1:nCells)

           ! call mpas_log_write('--- end interpolate height:')
           !CR: new vertical interpolation done! 27
        
           !... vertical velocity
            !CR: 28o Wind
            !CR: call interp_tofixed_pressure(nCells,nVertLevelsP1,nIntP,press_in,field_in,press_interp,field_interp)

            call VerticalInterpolation_Wind (nCells                                    , & !1. 
                                             nVertLevels                               , & !2. 
                                             nIsoLevelsT                               , & !3. 
                                             pressure2    (1:nVertLevels,1:nCells)*100.0, & !4. 
                                             t_iso_levels(1:nIsoLevelsT)               , & !5. 
                                             vvel        (1:nVertLevels,1:nCells)      , & !6. 
                                             field_interp(1:nIsoLevelsT,1:nCells)      )   !7.

            w_isobaric(1:nIsoLevelsT,1:nCells) =  field_interp(1:nIsoLevelsT,1:nCells)

           ! call mpas_log_write('--- end interpolate vertical velocity:')
           !CR: new vertical interpolation done! 28



           !... zgeo:

            if(allocated(field_in)) deallocate(field_in)
            !CR:if(.not.allocated(field_in)) allocate(field_in(nCells,nVertLevelsP1))

            !CR: 29o zgeo
            !CR: call interp_tofixed_pressure(nCells,nVertLevelsP1,nIntP,press_in,field_in,press_interp,field_interp)

            call VerticalInterpolation_ZGEO(nVertLevels                               , & !1
                                            nCells                                    , & !2
                                            temperature (1:nVertLevels,1:nCells)      , & !3
                                            pressure2    (1:nVertLevels,1:nCells)*100.0, & !4 
                                            height      (1:1,1:nCells)                , & !5 
                                            nIsoLevelsT                               , & !6
                                            field_interp(1:nIsoLevelsT,1:nCells)      , & !7
                                            t_iso_levels(1:nIsoLevelsT)               )   !8
            
            zgeo_isobaric(1:nIsoLevelsT,1:nCells) =  field_interp(1:nIsoLevelsT,1:nCells)


           ! call mpas_log_write('--- end interpolate zgeo:')
           !CR: new vertical interpolation done! 29


            !CR:? ...omega:

            if(allocated(field_in)) deallocate(field_in)
            !CR:if(.not.allocated(field_in)) allocate(field_in(nCells,nVertLevels))

            !CR: 30o Wind
            !CR: call interp_tofixed_pressure(nCells,nVertLevels,nIntP,press_in,field_in,press_interp,field_interp)

            call VerticalInterpolation_Wind (nCells                                    , & !1. 
                                             nVertLevels                               , & !2. 
                                             nIsoLevelsT                               , & !3. 
                                             pressure2    (1:nVertLevels,1:nCells)*100.0, & !4. 
                                             t_iso_levels(1:nIsoLevelsT)               , & !5. 
                                             rho         (1:nVertLevels, 1:nCells)     , & !6. 
                                             field_interp(1:nIsoLevelsT,1:nCells)      )   !7.


            omega_isobaric(1:nIsoLevelsT,1:nCells) =  w_isobaric(1:nIsoLevelsT,1:nCells) * field_interp(1:nIsoLevelsT,1:nCells) *grav *(-1)

           ! call mpas_log_write('--- end interpolat omega:')
           !CR: new vertical interpolation done! 30
            
            
            

            if(allocated(field_in)) deallocate(field_in)
            if(allocated(press_in)) deallocate(press_in)
        end if

        if(allocated(field_interp)) deallocate(field_interp)
        if(allocated(press_interp)) deallocate(press_interp)
       
        if (NEED_VORTICITY_ISOBARIC) then
           !interpolation to fixed pressure levels for fields located at cell vertices and at mass points:
            nIntP = 8
            if(.not.allocated(field_interp)) allocate(field_interp(nVertices,nIntP) )
            if(.not.allocated(press_interp)) allocate(press_interp(nVertices,nIntP) )
            do iVert = 1, nVertices
               press_interp(iVert,1) = 50.0_RKIND
               press_interp(iVert,2) = 100.0_RKIND
               press_interp(iVert,3) = 200.0_RKIND
               press_interp(iVert,4) = 250.0_RKIND
               press_interp(iVert,5) = 500.0_RKIND
               press_interp(iVert,6) = 700.0_RKIND
               press_interp(iVert,7) = 850.0_RKIND
               press_interp(iVert,8) = 925.0_RKIND
            enddo
       
            if(.not.allocated(press_in)) allocate(press_in(nVertices,nVertLevels))
            do iVert = 1, nVertices
            do k = 1, nVertLevels
               kk = nVertLevels+1-k
               press_in(iVert,kk) = pressure_v(k,iVert)
            enddo
            enddo

            if(.not.allocated(field_in)) allocate(field_in(nVertices,nVertLevels)) 
            !... relative vorticity:
            do iVert = 1, nVertices
            do k = 1, nVertLevels
               kk = nVertLevels+1-k
               field_in(iVert,kk) = vorticity(k,iVert)
            enddo
            enddo
            call interp_tofixed_pressure(nVertices,nVertLevels,nIntP,press_in,field_in,press_interp,field_interp)
            !CR: passar vorticity de 8n para 22n padrao, mas antes eh necessario calcular vorticity 
            !CR: no centro da celula, pois ela eh calculada nas arestas. Continua original MPAS ainda. (TO-DO)
               
            !vorticity_50hPa(1:nVertices) = field_interp(1:nVertices,1)
            !vorticity_100hPa(1:nVertices) = field_interp(1:nVertices,2)
            !vorticity_200hPa(1:nVertices) = field_interp(1:nVertices,3)
            !vorticity_250hPa(1:nVertices) = field_interp(1:nVertices,4)
            !vorticity_500hPa(1:nVertices) = field_interp(1:nVertices,5)
            !vorticity_700hPa(1:nVertices) = field_interp(1:nVertices,6)
            !vorticity_850hPa(1:nVertices) = field_interp(1:nVertices,7)
            !vorticity_925hPa(1:nVertices) = field_interp(1:nVertices,8)
           ! call mpas_log_write('--- end interpolate relative vorticity:')

            if(allocated(field_interp)) deallocate(field_interp)
            if(allocated(press_interp)) deallocate(press_interp)
       
            if(allocated(field_in    )) deallocate(field_in)
            if(allocated(press_in    )) deallocate(press_in)
        end if

        if(allocated(pressureCp1) ) deallocate(pressureCp1 )
        if(allocated(pressure_v)  ) deallocate(pressure_v  )
       
        if (need_mslp) then
            !... compute SLP (requires temp, height, pressure, qvapor)
             call compute_slp(nCells, nVertLevels, num_scalars, temperature, height, pressure, index_qv, scalars, mslp)
             mslp(:) = mslp(:) * 100.0   ! Convert from hPa to Pa
            !... alternative way
            !do iCell = 1, nCells
            !   mslp(iCell) = diag % surface_pressure % array(iCell) + 11.38*height(1,iCell)
            !   mslp(iCell) = mslp(iCell)/100.
            !enddo
        end if
    
    
        !!!!!!!!!!! Additional temperature levels for vortex tracking !!!!!!!!!!!
        if (need_t_isobaric .or. need_meanT_500_300) then
     
            allocate(field_in(nCells, nVertLevels))
            allocate(press_in(nCells, nVertLevels))
            !PK allocate(field_interp(nCells, nIsoLevelsT))
            !PK allocate(press_interp(nCells, nIsoLevelsT))
     
            !PK do k=1,nIsoLevelsT
            !PK   press_interp(:,k) = t_iso_levels(k)
            !PK end do
     
            ! Additional temperature levels for vortex tracking
            do iCell=1,nCells
            do k=1,nVertLevels
               kk = nVertLevels+1-k
               field_in(iCell,kk) = temperature(k,iCell)
            end do
            end do
     
            do iCell=1,nCells
            do k=1,nVertLevels
               kk = nVertLevels+1-k
               press_in(iCell,kk) = pressure(k,iCell) * 100.0
            end do
            end do
     
            if (need_t_isobaric) then
               !PK call interp_tofixed_pressure(nCells, nVertLevels, nIsoLevelsT, press_in, field_in, press_interp, field_interp)
               CALL VerticalInterpolation_Thermal2(nCells                                       ,&
                                                   nVertLevels                                  ,&
                                                   nIsoLevelsT                                  ,&
                                                   pressure2     (1:nVertLevels, 1:nCells)*100.0,&
                                                   temperature   (1:nVertLevels, 1:nCells)      ,&
                                                   t_isobaric    (1:nIsoLevelsT, 1:nCells)      ,&
                                                   t_iso_levels  (1:nIsoLevelsT)                 )
        
               !PK do k=1,nIsoLevelsT
               !PK    t_isobaric(k,1:nCells) = field_interp(k,1:nCells)
               !PK end do
            end if
     
     
            !!!!!!!!!!! Calculate mean temperature in 500 hPa - 300 hPa layer !!!!!!!!!!!
     
            if (need_meanT_500_300) then
                call compute_layer_mean(meanT_500_300, 50000.0_RKIND, 30000.0_RKIND, field_in, press_in)
            end if
     
     
            deallocate(field_in)
            !PK deallocate(field_interp)
            deallocate(press_in)
            !PK deallocate(press_interp)
        end if
     
     
        !!!!!!!!!!! Additional height levels for vortex tracking !!!!!!!!!!!
        if (need_z_isobaric) then
            !PK allocate(field_in(nCells, nVertLevelsP1))
            !PK allocate(press_in(nCells, nVertLevelsP1))

            !PK  allocate(field_interp(nCells, nIsoLevelsZ))
            !PK allocate(press_interp(nCells, nIsoLevelsZ))
     
            !do k=1,nIsoLevelsZ
            !   press_interp(:,k) = z_iso_levels(k)
            !end do
     
            !PK do iCell=1,nCells
            !PK do k=1,nVertLevelsP1
            !PK    kk = nVertLevelsP1+1-k
            !PK    field_in(iCell,kk) = height(k,iCell)
            !PK end do
            !PK end do
     
            !PK do iCell=1,nCells
            !PK do k=1,nVertLevelsP1
            !PK    kk = nVertLevelsP1+1-k
            !PK    press_in(iCell,kk) = pressure2(k,iCell) * 100.0
            !PK end do
            !PK end do
     
            !PK call interp_tofixed_pressure(nCells, nVertLevelsP1, nIsoLevelsZ, press_in, field_in, press_interp, field_interp)
              CALL VerticalInterpolation_Wind(nCells                                               ,&
                                              nVertLevels                                          ,&
                                              nIsoLevelsT                                          ,&
                                              pressure2      (1:nVertLevels        ,1:nCells)*100.0,&
                                              t_iso_levels   (1:nIsoLevelsT)                       ,&
                                              height         (1:nVertLevels        ,1:nCells)      ,&
                                              z_isobaric     (1:nIsoLevelsT        ,1:nCells)       )
      
            !PK do k=1,nIsoLevelsZ
            !PK   z_isobaric(k,1:nCells) = field_interp(1:nCells,k)
            !PK end do
     
            !PKdeallocate(field_in)
            !PKdeallocate(field_interp)
            !PKdeallocate(press_in)
            !PKdeallocate(press_interp)
        end if
    
        if(allocated(temperature) ) deallocate(temperature )
        if(allocated(pressure2)   ) deallocate(pressure2   )
        if(allocated(pressure)    ) deallocate(pressure    )
        if(allocated(dewpoint)    ) deallocate(dewpoint )
        if(allocated(spechum)     ) deallocate(spechum  )

    end subroutine interp_diagnostics


!***********************************************************************


   !==================================================================================================
    subroutine interp_tofixed_pressure(ncol,nlev_in,nlev_out,pres_in,field_in,pres_out,field_out)
   !==================================================================================================
   
   !input arguments:
    integer,intent(in):: ncol,nlev_in,nlev_out
   
    real(kind=RKIND),intent(in),dimension(ncol,nlev_in) :: pres_in,field_in
    real(kind=RKIND),intent(in),dimension(ncol,nlev_out):: pres_out
   
   !output arguments:
    real(kind=RKIND),intent(out),dimension(ncol,nlev_out):: field_out
   
   !local variables:
   ! integer:: i1,i2,icol,k,kk
    integer:: icol,k,kk
    integer:: kkstart,kount
    integer,dimension(ncol):: kupper
   
    real(kind=RKIND):: dpl,dpu
   
   !--------------------------------------------------------------------------------------------------
   
   !call mpas_log_write('')
   !call mpas_log_write('--- enter subroutine interp_tofixed_pressure:')
   !call mpas_log_write('... ncol     = $i',intArgs=(/ncol/))
   !call mpas_log_write('... nlev_in  = $i',intArgs=(/nlev_in/))
   !call mpas_log_write('... nlev_out = $i',intArgs=(/nlev_out/))
   !i1=1 ; i2=ncol
   !do k = 1, nlev_in
   !   call mpas_log_write('$i $r $r $r $r', intArgs=(/k/), realArgs=(/pres_in(i1,k),field_in(i1,k),pres_in(i2,k),field_in(i2,k)/))
   !enddo
   !call mpas_log_write('')
   
    do icol = 1, ncol
       kupper(icol) = 1
    enddo
   
    do k = 1, nlev_out 
   
       kkstart = nlev_in
       do icol = 1, ncol
          kkstart = min0(kkstart,kupper(icol))
       enddo
       kount = 0
   
       do kk = kkstart, nlev_in-1
          do icol = 1, ncol
             if(pres_out(icol,k).gt.pres_in(icol,kk).and.pres_out(icol,k).le.pres_in(icol,kk+1)) then
                kupper(icol) = kk
                kount = kount + 1
   !            call mpas_log_write('$i $i $r $r $r', intArgs=(/k,kupper(icol)/), realArgs=(/pres_out(icol,k),pres_in(icol,kk),pres_in(icol,kk+1)/))
             endif
          enddo
   
          if(kount.eq.ncol) then
             do icol = 1, ncol
                dpu = pres_out(icol,k) - pres_in(icol,kupper(icol))
                dpl = pres_in(icol,kupper(icol)+1) - pres_out(icol,k)
                field_out(icol,k) = (field_in(icol,kupper(icol))*dpl &
                                  + field_in(icol,kupper(icol)+1)*dpu)/(dpl + dpu)
             end do
             goto 35
           end if
       enddo
   
       do icol = 1, ncol
          if(pres_out(icol,k) .lt. pres_in(icol,1)) then
             field_out(icol,k) = field_in(icol,1)*pres_out(icol,k)/pres_in(icol,1)
          elseif(pres_out(icol,k) .gt. pres_in(icol,nlev_in)) then
             field_out(icol,k) = field_in(icol,nlev_in)
          else
             dpu = pres_out(icol,k) - pres_in(icol,kupper(icol))
             dpl = pres_in(icol,kupper(icol)+1) - pres_out(icol,k)
             field_out(icol,k) = (field_in(icol,kupper(icol))*dpl &
                               + field_in(icol,kupper(icol)+1)*dpu)/(dpl + dpu)
          endif
       enddo
   
    35 continue
   !   call mpas_log_write('$i $r $r $r $r $r $r', intArgs=(/kupper(i1)/), &
   !                       realArgs=(/pres_out(i1,k),pres_in(i1,kupper(i1)),pres_in(i1,kupper(i1)+1),field_out(i1,k),field_in(i1,kupper(i1)),field_in(i1,kupper(i1)+1)/))
   !   call mpas_log_write('$i $r $r $r $r $r $r', intArgs=(/kupper(i2)/), &
   !                       realArgs=(/pres_out(i2,k),pres_in(i2,kupper(i2)),pres_in(i2,kupper(i2)+1),field_out(i2,k),field_in(i2,kupper(i2)),field_in(i2,kupper(i2)+1)/))
   
    enddo
   
    end subroutine interp_tofixed_pressure
   

    subroutine compute_slp(ncol,nlev_in,nscalars,t,height,p,index_qv,scalars,slp)
   
       implicit none
   
      !input arguments:
       integer, intent(in) :: ncol, nlev_in, nscalars
      
      !p: in mb
      !t: in K
      !scalars: in kg/kg
      !height: in m
       real(kind=RKIND), intent(in), dimension(nlev_in,ncol) :: p,t
       real(kind=RKIND), intent(in), dimension(nlev_in+1,ncol) :: height
       integer, intent(in) :: index_qv
       real(kind=RKIND), intent(in), dimension(nscalars,nlev_in,ncol) :: scalars
      
      !output arguments:
       real(kind=RKIND), intent(out), dimension(ncol) :: slp
      
      !local variables:
       integer :: icol, k, kcount
       integer :: klo, khi
      
       real(kind=RKIND) :: gamma, rr, grav
       parameter (rr=287.0, grav=9.80616, gamma=0.0065)
      
       real(kind=RKIND) :: tc, pconst
       parameter (tc=273.16+17.5, pconst=100.)
      
       logical mm5_test
       parameter (mm5_test=.true.)
      
       integer, dimension(:), allocatable :: level
       real(kind=RKIND), dimension(:), allocatable :: t_surf, t_msl
       real(kind=RKIND) :: plo , phi , tlo, thi , zlo , zhi
       real(kind=RKIND) :: p_at_pconst , t_at_pconst , z_at_pconst, z_half_lowest
      
       logical :: l1, l2, l3, found
      
      ! Find least zeta level that is PCONST Pa above the surface.  We later use this
      ! level to extrapolate a surface pressure and temperature, which is supposed
      ! to reduce the effect of the diurnal heating cycle in the pressure field.
      
       if (.not.allocated(level))  allocate(level(ncol))
       if (.not.allocated(t_surf)) allocate(t_surf(ncol))
       if (.not.allocated(t_msl))  allocate(t_msl(ncol))
      
       do icol = 1 , ncol
          level(icol) = -1
      
          k = 1
          found = .false.
          do while ( (.not. found) .and. (k.le.nlev_in))
                if ( p(k,icol) .lt. p(1,icol)-pconst ) then
                   level(icol) = k
                   found = .true.
                end if
                k = k+1
          end do
      
          if ( level(icol) .eq. -1 ) then
             call mpas_log_write('Troubles finding level $r above ground.', realArgs=(/pconst/))
             call mpas_log_write('Problems first occur at ($i)', intArgs=(/icol/))
             call mpas_log_write('Surface pressure = $r hPa.', realArgs=(/p(1,icol)/))
             call mpas_log_write('*** MSLP field will not be computed')
             slp(:) = 0.0
             return
          end if
      
       end do
      
      ! Get temperature PCONST hPa above surface.  Use this to extrapolate
      ! the temperature at the surface and down to sea level.
      
       do icol = 1 , ncol
      
          klo = max ( level(icol) - 1 , 1      )
          khi = min ( klo + 1        , nlev_in - 1 )
      
          if ( klo .eq. khi ) then
             call mpas_log_write('Trapping levels are weird.')
             call mpas_log_write('icol = $i', intArgs=(/icol/))
             call mpas_log_write('klo = $i, khi = $i: and they should not be equal.', intArgs=(/klo,khi/))
             call mpas_log_write('*** MSLP field will not be computed')
             slp(:) = 0.0
             return
          end if
      
          plo = p(klo,icol)
          phi = p(khi,icol)
          tlo = t(klo,icol) * (1. + 0.608 * scalars(index_qv,klo,icol))
          thi = t(khi,icol) * (1. + 0.608 * scalars(index_qv,khi,icol))
          zlo = 0.5*(height(klo,icol)+height(klo+1,icol))
          zhi = 0.5*(height(khi,icol)+height(khi+1,icol))
      
          p_at_pconst = p(1,icol) - pconst
          t_at_pconst = thi-(thi-tlo)*log(p_at_pconst/phi)*log(plo/phi)
          z_at_pconst = zhi-(zhi-zlo)*log(p_at_pconst/phi)*log(plo/phi)
      
          t_surf(icol) = t_at_pconst*(p(1,icol)/p_at_pconst)**(gamma*rr/grav)
          t_msl(icol) = t_at_pconst+gamma*z_at_pconst
      
       end do
      
      ! If we follow a traditional computation, there is a correction to the sea level
      ! temperature if both the surface and sea level temnperatures are *too* hot.
      
       if ( mm5_test ) then
          kcount = 0
          do icol = 1 , ncol
                l1 = t_msl(icol) .lt. tc
                l2 = t_surf(icol) .le. tc
                l3 = .not. l1
                if ( l2 .and. l3 ) then
                   t_msl(icol) = tc
                else
                   t_msl(icol) = tc - 0.005*(t_surf(icol)-tc)**2
                   kcount = kcount+1
                end if
          end do
      !   call mpas_log_write('These number of points had t_msl adjusted $i', intArgs=(/kcount/))
       end if
      
       do icol = 1 , ncol
          z_half_lowest=0.5*(height(1,icol)+height(2,icol))
          slp(icol) = p(1,icol) * exp((2.*grav*z_half_lowest)/ &
                                    (rr*(t_msl(icol)+t_surf(icol))))
       end do
      
       if (allocated(level))  deallocate(level)
       if (allocated(t_surf)) deallocate(t_surf)
       if (allocated(t_msl))  deallocate(t_msl)
   
    end subroutine compute_slp


   !***********************************************************************
   !
   !  routine compute_layer_mean
   !
   !> \brief   Computes the mean of a field in the specified layer.
   !> \author  Michael Duda
   !> \date    3 July 2014
   !> \details
   !>  Given a 3d pressure field, press_in(nCells,nVertLevels), with pressure 
   !>  increasing with vertical index, and a 3d field, 
   !>  field_in(nCells,nVertLevels) with levels in the same order, this routine
   !>  will compute the mean of the field for each column between pressures
   !>  p1 and p2.
   !
   !----------------------------------------------------------------------- 
    subroutine compute_layer_mean(layerMean, p1, p2, field_in, press_in)
   
       implicit none
   
       real(kind=RKIND), dimension(:), intent(out) :: layerMean
       real(kind=RKIND), intent(in) :: p1, p2
       real(kind=RKIND), dimension(:,:), intent(in) :: field_in
       real(kind=RKIND), dimension(:,:), intent(in) :: press_in
   
       integer :: nCells, nVertLevels
       integer :: iCell, k
       integer :: k_bot, k_top
       real(kind=RKIND) :: p_bot, p_top
       real(kind=RKIND) :: wtop_p, wtop_m
       real(kind=RKIND) :: wbot_p, wbot_m
       real(kind=RKIND) :: wtotal, w
       real(kind=RKIND) :: temp
   
   
       !
       ! Get dimensions of input arrays
       !
       nCells = size(field_in, 1) 
       nVertLevels = size(field_in, 2) 
   
   
       !
       ! Check that pressure is increasing with index
       !
       if (press_in(1,1) > press_in(1,nVertLevels)) then
           call mpas_log_write('Error in compute_layer_mean: pressure should increase with index', messageType=MPAS_LOG_ERR)
           layerMean(:) = 0.0
           return
       end if
   
       
       !
       ! Set the pressure at the top and bottom of the layer
       !
       if (p1 < p2) then
          p_top = p1
          p_bot = p2
       else
          p_top = p2
          p_bot = p1
       end if
   
   
       !
       ! For each column, compute the mean value of the field between p_bot and
       ! p_top, with the field weighted by delta-p in each layer
       !
       do iCell=1,nCells
          k_bot = -1
          k_top = -1
   
          ! Search for trapping levels: k_top is the index just above (or equal to)
          ! p_top, and k_bot is the index just above (or equal to) p_bot.
          do k=1,nVertLevels-1
             if (press_in(iCell,k) <= p_top .and. press_in(iCell,k+1) > p_top) then
                k_top = k
                wtop_p = (p_top - press_in(iCell,k)) / (press_in(iCell,k+1) - press_in(iCell,k))
                wtop_m = (press_in(iCell,k+1) - p_top) / (press_in(iCell,k+1) - press_in(iCell,k))
             end if
             if (press_in(iCell,k) <= p_bot .and. press_in(iCell,k+1) > p_bot) then
                k_bot = k
                wbot_m = (p_bot - press_in(iCell,k)) / (press_in(iCell,k+1) - press_in(iCell,k))
                wbot_p = (press_in(iCell,k+1) - p_bot) / (press_in(iCell,k+1) - press_in(iCell,k))
             end if
          end do
   
          if (k_top == -1 .or. k_bot == -1) then      ! Layer intersects top or bottom boundary
   
             layerMean(iCell) = 0.0  
   
          else if (k_top == k_bot) then               ! Layer lies entirely within a single model layer
   
             layerMean(iCell) = wtop_m * field_in(iCell,k_top) + wtop_p * field_in(iCell,k_top+1)
             layerMean(iCell) = layerMean(iCell) + wbot_m * field_in(iCell,k_bot) + wbot_p * field_in(iCell,k_bot+1)
             layerMean(iCell) = 0.5 * layerMean(iCell)
   
          else
   
             ! First layer: from p_top down to press_in(iCell,k_top+1)
             wtotal = press_in(iCell,k_top+1) - p_top  
             temp = wtop_m * field_in(iCell,k_top) + wtop_p * field_in(iCell,k_top+1)
             layerMean(iCell) = wtotal * 0.5 * (field_in(iCell,k_top+1) + temp)
   
             ! Middle layers
             do k=k_top+1,k_bot-1
                w = press_in(iCell,k+1) - press_in(iCell,k)
                wtotal = wtotal + w
                layerMean(iCell) = layerMean(iCell) + w * 0.5 * (field_in(iCell,k) + field_in(iCell,k+1))
             end do
   
             ! Last layer: from press_in(iCell,k_bot) down to p_bot
             w = p_bot - press_in(iCell,k_bot)
             wtotal = wtotal + w
             temp = wbot_m * field_in(iCell,k_bot) + wbot_p * field_in(iCell,k_bot+1)
             layerMean(iCell) = layerMean(iCell) + w * 0.5 * (field_in(iCell,k_bot) + temp)
   
             layerMean(iCell) = layerMean(iCell) / wtotal
          end if
   
       end do
   
    end subroutine compute_layer_mean
   
end module mpas_isobaric_diagnostics
