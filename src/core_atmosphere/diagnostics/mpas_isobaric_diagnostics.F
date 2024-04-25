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

    logical :: need_mslp, &
               need_cldfrac_15, need_cldfrac_20, need_cldfrac_30, need_cldfrac_50, need_cldfrac_70, need_cldfrac_100, need_cldfrac_150,need_cldfrac_200, need_cldfrac_250, need_cldfrac_300, need_cldfrac_400, need_cldfrac_500, need_cldfrac_600,need_cldfrac_700, need_cldfrac_825, need_cldfrac_850, need_cldfrac_875, need_cldfrac_900, need_cldfrac_925,need_cldfrac_950, need_cldfrac_975, need_cldfrac_1000, &
               need_relhum_15, need_relhum_20, need_relhum_30, need_relhum_50, need_relhum_70, need_relhum_100, need_relhum_150, need_relhum_200, need_relhum_250, need_relhum_300, need_relhum_400, need_relhum_500, need_relhum_600, need_relhum_700, need_relhum_825, need_relhum_850, need_relhum_875, need_relhum_900, need_relhum_925, need_relhum_950, need_relhum_975, need_relhum_1000, &
               need_dewpoint_50, need_dewpoint_100, need_dewpoint_200, need_dewpoint_250, need_dewpoint_500, need_dewpoint_700, need_dewpoint_850, need_dewpoint_925, &
               need_temp_15, need_temp_20, need_temp_30, need_temp_50, need_temp_70, need_temp_100, need_temp_150, need_temp_200, need_temp_250, need_temp_300, need_temp_400, need_temp_500, need_temp_600, need_temp_700, need_temp_825, need_temp_850, need_temp_875, need_temp_900, need_temp_925, need_temp_950, need_temp_975, need_temp_1000, &
               need_height_15, need_height_20, need_height_30, need_height_50, need_height_70, need_height_100, need_height_150, need_height_200, need_height_250, need_height_300, need_height_400, need_height_500, need_height_600, need_height_700, need_height_825, need_height_850, need_height_875, need_height_900, need_height_925, need_height_950, need_height_975, need_height_1000, &
               need_uzonal_15, need_uzonal_20, need_uzonal_30, need_uzonal_50, need_uzonal_70, need_uzonal_100, need_uzonal_150, need_uzonal_200, need_uzonal_250, need_uzonal_300, need_uzonal_400, need_uzonal_500, need_uzonal_600, need_uzonal_700, need_uzonal_825, need_uzonal_850, need_uzonal_875, need_uzonal_900, need_uzonal_925, need_uzonal_950, need_uzonal_975, need_uzonal_1000, &
               need_umeridional_15, need_umeridional_20, need_umeridional_30, need_umeridional_50, need_umeridional_70, need_umeridional_100, need_umeridional_150,need_umeridional_200,need_umeridional_250, need_umeridional_300, need_umeridional_400, need_umeridional_500, need_umeridional_600, need_umeridional_700, need_umeridional_825, need_umeridional_850, need_umeridional_875, need_umeridional_900, need_umeridional_925, need_umeridional_950, need_umeridional_975, need_umeridional_1000, &
               need_spechum_15, need_spechum_20, need_spechum_30, need_spechum_50, need_spechum_70, need_spechum_100, need_spechum_150, need_spechum_200, need_spechum_250, need_spechum_300, need_spechum_400, need_spechum_500, need_spechum_600, need_spechum_700, need_spechum_825, need_spechum_850, need_spechum_875, need_spechum_900, need_spechum_925, need_spechum_950, need_spechum_975, need_spechum_1000, &               
               need_w_15, need_w_20, need_w_30, need_w_50, need_w_70, need_w_100, need_w_150, need_w_200, need_w_250, need_w_300, need_w_400, need_w_500, need_w_600, need_w_700, need_w_825, need_w_850, need_w_875, need_w_900, need_w_925, need_w_950, need_w_975, need_w_1000, &
               need_geoph_15, need_geoph_20, need_geoph_30, need_geoph_50, need_geoph_70, need_geoph_100, need_geoph_150, need_geoph_200, need_geoph_250, need_geoph_300, need_geoph_400, need_geoph_500, need_geoph_600, need_geoph_700, need_geoph_825, need_geoph_850, need_geoph_875, need_geoph_900, need_geoph_925, need_geoph_950, need_geoph_975, need_geoph_1000, &
               need_rqvblten_15, need_rqvblten_20, need_rqvblten_30, need_rqvblten_50, need_rqvblten_70, need_rqvblten_100, need_rqvblten_150, need_rqvblten_200, need_rqvblten_250, need_rqvblten_300, need_rqvblten_400, need_rqvblten_500, need_rqvblten_600, need_rqvblten_700, need_rqvblten_825, need_rqvblten_850, need_rqvblten_875, need_rqvblten_900, need_rqvblten_925, need_rqvblten_950, need_rqvblten_975, need_rqvblten_1000, &
               need_rqvcuten_15, need_rqvcuten_20, need_rqvcuten_30, need_rqvcuten_50, need_rqvcuten_70, need_rqvcuten_100, need_rqvcuten_150, need_rqvcuten_200, need_rqvcuten_250, need_rqvcuten_300, need_rqvcuten_400, need_rqvcuten_500, need_rqvcuten_600, need_rqvcuten_700, need_rqvcuten_825, need_rqvcuten_850, need_rqvcuten_875, need_rqvcuten_900, need_rqvcuten_925, need_rqvcuten_950, need_rqvcuten_975, need_rqvcuten_1000, &
               need_rthratenlw_15, need_rthratenlw_20, need_rthratenlw_30, need_rthratenlw_50, need_rthratenlw_70, need_rthratenlw_100, need_rthratenlw_150, need_rthratenlw_200, need_rthratenlw_250, need_rthratenlw_300, need_rthratenlw_400, need_rthratenlw_500, need_rthratenlw_600, need_rthratenlw_700, need_rthratenlw_825, need_rthratenlw_850, need_rthratenlw_875, need_rthratenlw_900, need_rthratenlw_925, need_rthratenlw_950, need_rthratenlw_975, need_rthratenlw_1000, &
               need_rucuten_15, need_rucuten_20, need_rucuten_30, need_rucuten_50, need_rucuten_70, need_rucuten_100, need_rucuten_150, need_rucuten_200, need_rucuten_250, need_rucuten_300, need_rucuten_400, need_rucuten_500, need_rucuten_600, need_rucuten_700, need_rucuten_825, need_rucuten_850, need_rucuten_875, need_rucuten_900, need_rucuten_925, need_rucuten_950, need_rucuten_975, need_rucuten_1000, &
               need_kzq_15, need_kzq_20, need_kzq_30, need_kzq_50, need_kzq_70, need_kzq_100, need_kzq_150, need_kzq_200, need_kzq_250, need_kzq_300, need_kzq_400, need_kzq_500, need_kzq_600, need_kzq_700, need_kzq_825, need_kzq_850, need_kzq_875, need_kzq_900, need_kzq_925, need_kzq_950, need_kzq_975, need_kzq_1000, &
               need_kzh_15, need_kzh_20, need_kzh_30, need_kzh_50, need_kzh_70, need_kzh_100, need_kzh_150, need_kzh_200, need_kzh_250, need_kzh_300, need_kzh_400, need_kzh_500, need_kzh_600, need_kzh_700, need_kzh_825, need_kzh_850, need_kzh_875, need_kzh_900, need_kzh_925, need_kzh_950, need_kzh_975, need_kzh_1000, &
               need_kzm_15, need_kzm_20, need_kzm_30, need_kzm_50, need_kzm_70, need_kzm_100, need_kzm_150, need_kzm_200, need_kzm_250, need_kzm_300, need_kzm_400, need_kzm_500, need_kzm_600, need_kzm_700, need_kzm_825, need_kzm_850, need_kzm_875, need_kzm_900, need_kzm_925, need_kzm_950, need_kzm_975, need_kzm_1000, &
               need_ni_15, need_ni_20, need_ni_30, need_ni_50, need_ni_70, need_ni_100, need_ni_150, need_ni_200, need_ni_250, need_ni_300, need_ni_400, need_ni_500, need_ni_600, need_ni_700, need_ni_825, need_ni_850, need_ni_875, need_ni_900, need_ni_925, need_ni_950, need_ni_975, need_ni_1000, &
               need_nr_15, need_nr_20, need_nr_30, need_nr_50, need_nr_70, need_nr_100, need_nr_150, need_nr_200, need_nr_250, need_nr_300, need_nr_400, need_nr_500, need_nr_600, need_nr_700, need_nr_825, need_nr_850, need_nr_875, need_nr_900, need_nr_925, need_nr_950, need_nr_975, need_nr_1000, &
               need_qr_15, need_qr_20, need_qr_30, need_qr_50, need_qr_70, need_qr_100, need_qr_150, need_qr_200, need_qr_250, need_qr_300, need_qr_400, need_qr_500, need_qr_600, need_qr_700, need_qr_825, need_qr_850, need_qr_875, need_qr_900, need_qr_925, need_qr_950, need_qr_975, need_qr_1000, &
               need_qs_15, need_qs_20, need_qs_30, need_qs_50, need_qs_70, need_qs_100, need_qs_150, need_qs_200, need_qs_250, need_qs_300, need_qs_400, need_qs_500, need_qs_600, need_qs_700, need_qs_825, need_qs_850, need_qs_875, need_qs_900, need_qs_925, need_qs_950, need_qs_975, need_qs_1000, &
               need_qv_15, need_qv_20, need_qv_30, need_qv_50, need_qv_70, need_qv_100, need_qv_150, need_qv_200, need_qv_250, need_qv_300, need_qv_400, need_qv_500, need_qv_600, need_qv_700, need_qv_825, need_qv_850, need_qv_875, need_qv_900, need_qv_925, need_qv_950, need_qv_975, need_qv_1000, &
               need_qc_15, need_qc_20, need_qc_30, need_qc_50, need_qc_70, need_qc_100, need_qc_150, need_qc_200, need_qc_250, need_qc_300, need_qc_400, need_qc_500, need_qc_600, need_qc_700, need_qc_825, need_qc_850, need_qc_875, need_qc_900, need_qc_925, need_qc_950, need_qc_975, need_qc_1000, &
               need_qg_15, need_qg_20, need_qg_30, need_qg_50, need_qg_70, need_qg_100, need_qg_150, need_qg_200, need_qg_250, need_qg_300, need_qg_400, need_qg_500, need_qg_600, need_qg_700, need_qg_825, need_qg_850, need_qg_875, need_qg_900, need_qg_925, need_qg_950, need_qg_975, need_qg_1000, &
               need_qi_15, need_qi_20, need_qi_30, need_qi_50, need_qi_70, need_qi_100, need_qi_150, need_qi_200, need_qi_250, need_qi_300, need_qi_400, need_qi_500, need_qi_600, need_qi_700, need_qi_825, need_qi_850, need_qi_875, need_qi_900, need_qi_925, need_qi_950, need_qi_975, need_qi_1000, &
               need_vorticity_50, need_vorticity_100, need_vorticity_200, need_vorticity_250, need_vorticity_500, need_vorticity_700, need_vorticity_850, need_vorticity_925, &
               need_t_isobaric, need_z_isobaric, need_meanT_500_300
    logical :: need_temp, need_relhum, need_dewpoint, need_w, need_uzonal, need_umeridional, need_vorticity, & 
               need_height, need_spechum, need_geoph, need_rqvblten, need_rqvcuten, need_rthratenlw, need_rucuten, &
               need_kzq, need_kzh, need_kzm, need_cldfrac, & 
               need_ni, need_nr, need_qr, need_qs, need_qv, need_qc, need_qg, need_qi, &
               need_mp_thompson, need_cu_grell_freitas, need_bl_mynn


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

        print *, "mp_thompson=", need_mp_thompson
        print *, "cu_grell_freitas=", need_cu_grell_freitas
        print *, "bl_mynn=", need_bl_mynn

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

        need_cldfrac = .false.
        need_temp = .false.
        need_dewpoint = .false.
        need_relhum = .false.
        need_w = .false.
        need_uzonal = .false.
        need_umeridional = .false.
        need_vorticity = .false.
        need_height = .false.
        need_spechum = .false.
        need_geoph = .false.
        need_rqvblten = .false.
        need_rqvcuten = .false.
        need_rthratenlw = .false.
        need_rucuten = .false.
        need_kzq = .false.
        need_kzh = .false.
        need_kzm = .false.
        need_ni = .false.
        need_nr = .false.
        need_qr = .false.
        need_qs = .false.
        need_qv = .false.
        need_qc = .false.
        need_qg = .false.
        need_qi = .false.

        need_mslp = MPAS_field_will_be_written('mslp')
        need_any_diags = need_any_diags .or. need_mslp
        
        need_relhum_15 = MPAS_field_will_be_written('relhum_15hPa')
        need_relhum = need_relhum .or. need_relhum_15
        need_any_diags = need_any_diags .or. need_relhum_15
        need_relhum_20 = MPAS_field_will_be_written('relhum_20hPa')
        need_relhum = need_relhum .or. need_relhum_20
        need_any_diags = need_any_diags .or. need_relhum_20
        need_relhum_30 = MPAS_field_will_be_written('relhum_30hPa')
        need_relhum = need_relhum .or. need_relhum_30
        need_any_diags = need_any_diags .or. need_relhum_30
        need_relhum_50 = MPAS_field_will_be_written('relhum_50hPa')
        need_relhum = need_relhum .or. need_relhum_50
        need_any_diags = need_any_diags .or. need_relhum_50
        need_relhum_70 = MPAS_field_will_be_written('relhum_70hPa')
        need_relhum = need_relhum .or. need_relhum_70
        need_any_diags = need_any_diags .or. need_relhum_70
        need_relhum_100 = MPAS_field_will_be_written('relhum_100hPa')
        need_relhum = need_relhum .or. need_relhum_100
        need_any_diags = need_any_diags .or. need_relhum_100
        need_relhum_150 = MPAS_field_will_be_written('relhum_150hPa')
        need_relhum = need_relhum .or. need_relhum_150
        need_any_diags = need_any_diags .or. need_relhum_150
        need_relhum_200 = MPAS_field_will_be_written('relhum_200hPa')
        need_relhum = need_relhum .or. need_relhum_200
        need_any_diags = need_any_diags .or. need_relhum_200
        need_relhum_250 = MPAS_field_will_be_written('relhum_250hPa')
        need_relhum = need_relhum .or. need_relhum_250
        need_any_diags = need_any_diags .or. need_relhum_250
        need_relhum_300 = MPAS_field_will_be_written('relhum_300hPa')
        need_relhum = need_relhum .or. need_relhum_300
        need_any_diags = need_any_diags .or. need_relhum_300
        need_relhum_400 = MPAS_field_will_be_written('relhum_400hPa')
        need_relhum = need_relhum .or. need_relhum_400
        need_any_diags = need_any_diags .or. need_relhum_400
        need_relhum_500 = MPAS_field_will_be_written('relhum_500hPa')
        need_relhum = need_relhum .or. need_relhum_500
        need_any_diags = need_any_diags .or. need_relhum_500
        need_relhum_600 = MPAS_field_will_be_written('relhum_600hPa')
        need_relhum = need_relhum .or. need_relhum_600
        need_any_diags = need_any_diags .or. need_relhum_600
        need_relhum_700 = MPAS_field_will_be_written('relhum_700hPa')
        need_relhum = need_relhum .or. need_relhum_700
        need_any_diags = need_any_diags .or. need_relhum_700
        need_relhum_825 = MPAS_field_will_be_written('relhum_825hPa')
        need_relhum = need_relhum .or. need_relhum_825
        need_any_diags = need_any_diags .or. need_relhum_825
        need_relhum_850 = MPAS_field_will_be_written('relhum_850hPa')
        need_relhum = need_relhum .or. need_relhum_850
        need_any_diags = need_any_diags .or. need_relhum_850
        need_relhum_875 = MPAS_field_will_be_written('relhum_875hPa')
        need_relhum = need_relhum .or. need_relhum_875
        need_any_diags = need_any_diags .or. need_relhum_875
        need_relhum_900 = MPAS_field_will_be_written('relhum_900hPa')
        need_relhum = need_relhum .or. need_relhum_900
        need_any_diags = need_any_diags .or. need_relhum_900
        need_relhum_925 = MPAS_field_will_be_written('relhum_925hPa')
        need_relhum = need_relhum .or. need_relhum_925
        need_any_diags = need_any_diags .or. need_relhum_925
        need_relhum_950 = MPAS_field_will_be_written('relhum_950hPa')
        need_relhum = need_relhum .or. need_relhum_950
        need_any_diags = need_any_diags .or. need_relhum_950
        need_relhum_975 = MPAS_field_will_be_written('relhum_975hPa')
        need_relhum = need_relhum .or. need_relhum_975
        need_any_diags = need_any_diags .or. need_relhum_975
        need_relhum_1000 = MPAS_field_will_be_written('relhum_1000hPa')
        need_relhum = need_relhum .or. need_relhum_1000
        need_any_diags = need_any_diags .or. need_relhum_1000
        need_dewpoint_50 = MPAS_field_will_be_written('dewpoint_50hPa')
        need_dewpoint = need_dewpoint .or. need_dewpoint_50
        need_any_diags = need_any_diags .or. need_dewpoint_50
        need_dewpoint_100 = MPAS_field_will_be_written('dewpoint_100hPa')
        need_dewpoint = need_dewpoint .or. need_dewpoint_100
        need_any_diags = need_any_diags .or. need_dewpoint_100
        need_dewpoint_200 = MPAS_field_will_be_written('dewpoint_200hPa')
        need_dewpoint = need_dewpoint .or. need_dewpoint_200
        need_any_diags = need_any_diags .or. need_dewpoint_200
        need_dewpoint_250 = MPAS_field_will_be_written('dewpoint_250hPa')
        need_dewpoint = need_dewpoint .or. need_dewpoint_250
        need_any_diags = need_any_diags .or. need_dewpoint_250
        need_dewpoint_500 = MPAS_field_will_be_written('dewpoint_500hPa')
        need_dewpoint = need_dewpoint .or. need_dewpoint_500
        need_any_diags = need_any_diags .or. need_dewpoint_500
        need_dewpoint_700 = MPAS_field_will_be_written('dewpoint_700hPa')
        need_dewpoint = need_dewpoint .or. need_dewpoint_700
        need_any_diags = need_any_diags .or. need_dewpoint_700
        need_dewpoint_850 = MPAS_field_will_be_written('dewpoint_850hPa')
        need_dewpoint = need_dewpoint .or. need_dewpoint_850
        need_any_diags = need_any_diags .or. need_dewpoint_850
        need_dewpoint_925 = MPAS_field_will_be_written('dewpoint_925hPa')
        need_dewpoint = need_dewpoint .or. need_dewpoint_925
        need_any_diags = need_any_diags .or. need_dewpoint_925

        need_cldfrac_15 = MPAS_field_will_be_written('cldfrac_15hPa')
        need_cldfrac = need_cldfrac .or. need_cldfrac_15
        need_any_diags = need_any_diags .or. need_cldfrac_15
        need_cldfrac_20 = MPAS_field_will_be_written('cldfrac_20hPa')
        need_cldfrac = need_cldfrac .or. need_cldfrac_20
        need_any_diags = need_any_diags .or. need_cldfrac_20
        need_cldfrac_30 = MPAS_field_will_be_written('cldfrac_30hPa')
        need_cldfrac = need_cldfrac .or. need_cldfrac_30
        need_any_diags = need_any_diags .or. need_cldfrac_30
        need_cldfrac_50 = MPAS_field_will_be_written('cldfrac_50hPa')
        need_cldfrac = need_cldfrac .or. need_cldfrac_50
        need_any_diags = need_any_diags .or. need_cldfrac_50
        need_cldfrac_70 = MPAS_field_will_be_written('cldfrac_70hPa')
        need_cldfrac = need_cldfrac .or. need_cldfrac_70
        need_any_diags = need_any_diags .or. need_cldfrac_70
        need_cldfrac_100 = MPAS_field_will_be_written('cldfrac_100hPa')
        need_cldfrac = need_cldfrac .or. need_cldfrac_100
        need_any_diags = need_any_diags .or. need_cldfrac_100
        need_cldfrac_150 = MPAS_field_will_be_written('cldfrac_150hPa')
        need_cldfrac = need_cldfrac .or. need_cldfrac_150
        need_any_diags = need_any_diags .or. need_cldfrac_150
        need_cldfrac_200 = MPAS_field_will_be_written('cldfrac_200hPa')
        need_cldfrac = need_cldfrac .or. need_cldfrac_200
        need_any_diags = need_any_diags .or. need_cldfrac_200
        need_cldfrac_250 = MPAS_field_will_be_written('cldfrac_250hPa')
        need_cldfrac = need_cldfrac .or. need_cldfrac_250
        need_any_diags = need_any_diags .or. need_cldfrac_250
        need_cldfrac_300 = MPAS_field_will_be_written('cldfrac_300hPa')
        need_cldfrac = need_cldfrac .or. need_cldfrac_300
        need_any_diags = need_any_diags .or. need_cldfrac_300
        need_cldfrac_400 = MPAS_field_will_be_written('cldfrac_400hPa')
        need_cldfrac = need_cldfrac .or. need_cldfrac_400
        need_any_diags = need_any_diags .or. need_cldfrac_400
        need_cldfrac_500 = MPAS_field_will_be_written('cldfrac_500hPa')
        need_cldfrac = need_cldfrac .or. need_cldfrac_500
        need_any_diags = need_any_diags .or. need_cldfrac_500
        need_cldfrac_600 = MPAS_field_will_be_written('cldfrac_600hPa')
        need_cldfrac = need_cldfrac .or. need_cldfrac_600
        need_any_diags = need_any_diags .or. need_cldfrac_600
        need_cldfrac_700 = MPAS_field_will_be_written('cldfrac_700hPa')
        need_cldfrac = need_cldfrac .or. need_cldfrac_700
        need_any_diags = need_any_diags .or. need_cldfrac_700
        need_cldfrac_825 = MPAS_field_will_be_written('cldfrac_825hPa')
        need_cldfrac = need_cldfrac .or. need_cldfrac_825
        need_any_diags = need_any_diags .or. need_cldfrac_825
        need_cldfrac_850 = MPAS_field_will_be_written('cldfrac_850hPa')
        need_cldfrac = need_cldfrac .or. need_cldfrac_850
        need_any_diags = need_any_diags .or. need_cldfrac_850
        need_cldfrac_875 = MPAS_field_will_be_written('cldfrac_875hPa')
        need_cldfrac = need_cldfrac .or. need_cldfrac_875
        need_any_diags = need_any_diags .or. need_cldfrac_875
        need_cldfrac_900 = MPAS_field_will_be_written('cldfrac_900hPa')
        need_cldfrac = need_cldfrac .or. need_cldfrac_900
        need_any_diags = need_any_diags .or. need_cldfrac_900
        need_cldfrac_925 = MPAS_field_will_be_written('cldfrac_925hPa')
        need_cldfrac = need_cldfrac .or. need_cldfrac_925
        need_any_diags = need_any_diags .or. need_cldfrac_925
        need_cldfrac_950 = MPAS_field_will_be_written('cldfrac_950hPa')
        need_cldfrac = need_cldfrac .or. need_cldfrac_950
        need_any_diags = need_any_diags .or. need_cldfrac_950
        need_cldfrac_975 = MPAS_field_will_be_written('cldfrac_975hPa')
        need_cldfrac = need_cldfrac .or. need_cldfrac_975
        need_any_diags = need_any_diags .or. need_cldfrac_975
        need_cldfrac_1000 = MPAS_field_will_be_written('cldfrac_1000hPa')
        need_cldfrac = need_cldfrac .or. need_cldfrac_1000
        need_any_diags = need_any_diags .or. need_cldfrac_1000


        need_temp_15 = MPAS_field_will_be_written('temperature_15hPa')
        need_temp = need_temp .or. need_temp_15
        need_any_diags = need_any_diags .or. need_temp_15
        need_temp_20 = MPAS_field_will_be_written('temperature_20hPa')
        need_temp = need_temp .or. need_temp_20
        need_any_diags = need_any_diags .or. need_temp_20
        need_temp_30 = MPAS_field_will_be_written('temperature_30hPa')
        need_temp = need_temp .or. need_temp_30
        need_any_diags = need_any_diags .or. need_temp_30
        need_temp_50 = MPAS_field_will_be_written('temperature_50hPa')
        need_temp = need_temp .or. need_temp_50
        need_any_diags = need_any_diags .or. need_temp_50
        need_temp_70 = MPAS_field_will_be_written('temperature_70hPa')
        need_temp = need_temp .or. need_temp_70
        need_any_diags = need_any_diags .or. need_temp_70
        need_temp_100 = MPAS_field_will_be_written('temperature_100hPa')
        need_temp = need_temp .or. need_temp_100
        need_any_diags = need_any_diags .or. need_temp_100
        need_temp_150 = MPAS_field_will_be_written('temperature_150hPa')
        need_temp = need_temp .or. need_temp_150
        need_any_diags = need_any_diags .or. need_temp_150
        need_temp_200 = MPAS_field_will_be_written('temperature_200hPa')
        need_temp = need_temp .or. need_temp_200
        need_any_diags = need_any_diags .or. need_temp_200
        need_temp_250 = MPAS_field_will_be_written('temperature_250hPa')
        need_temp = need_temp .or. need_temp_250
        need_any_diags = need_any_diags .or. need_temp_250
        need_temp_300 = MPAS_field_will_be_written('temperature_300hPa')
        need_temp = need_temp .or. need_temp_300
        need_any_diags = need_any_diags .or. need_temp_300
        need_temp_400 = MPAS_field_will_be_written('temperature_400hPa')
        need_temp = need_temp .or. need_temp_400
        need_any_diags = need_any_diags .or. need_temp_400
        need_temp_500 = MPAS_field_will_be_written('temperature_500hPa')
        need_temp = need_temp .or. need_temp_500
        need_any_diags = need_any_diags .or. need_temp_500
        need_temp_600 = MPAS_field_will_be_written('temperature_600hPa')
        need_temp = need_temp .or. need_temp_600
        need_any_diags = need_any_diags .or. need_temp_600
        need_temp_700 = MPAS_field_will_be_written('temperature_700hPa')
        need_temp = need_temp .or. need_temp_700
        need_any_diags = need_any_diags .or. need_temp_700
        need_temp_825 = MPAS_field_will_be_written('temperature_825hPa')
        need_temp = need_temp .or. need_temp_825
        need_any_diags = need_any_diags .or. need_temp_825
        need_temp_850 = MPAS_field_will_be_written('temperature_850hPa')
        need_temp = need_temp .or. need_temp_850
        need_any_diags = need_any_diags .or. need_temp_850
        need_temp_875 = MPAS_field_will_be_written('temperature_875hPa')
        need_temp = need_temp .or. need_temp_875
        need_any_diags = need_any_diags .or. need_temp_875
        need_temp_900 = MPAS_field_will_be_written('temperature_900hPa')
        need_temp = need_temp .or. need_temp_900
        need_any_diags = need_any_diags .or. need_temp_900
        need_temp_925 = MPAS_field_will_be_written('temperature_925hPa')
        need_temp = need_temp .or. need_temp_925
        need_any_diags = need_any_diags .or. need_temp_925
        need_temp_950 = MPAS_field_will_be_written('temperature_950hPa')
        need_temp = need_temp .or. need_temp_950
        need_any_diags = need_any_diags .or. need_temp_950
        need_temp_975 = MPAS_field_will_be_written('temperature_975hPa')
        need_temp = need_temp .or. need_temp_975
        need_any_diags = need_any_diags .or. need_temp_975
        need_temp_1000 = MPAS_field_will_be_written('temperature_1000hPa')
        need_temp = need_temp .or. need_temp_1000
        need_any_diags = need_any_diags .or. need_temp_1000

        need_height_15 = MPAS_field_will_be_written('height_15hPa')
        need_height = need_height .or. need_height_15
        need_any_diags = need_any_diags .or. need_height_15

        need_height_20 = MPAS_field_will_be_written('height_20hPa')
        need_height = need_height .or. need_height_20
        need_any_diags = need_any_diags .or. need_height_20

        need_height_30 = MPAS_field_will_be_written('height_30hPa')
        need_height = need_height .or. need_height_30
        need_any_diags = need_any_diags .or. need_height_30
        need_height_50 = MPAS_field_will_be_written('height_50hPa')
        need_height = need_height .or. need_height_50
        need_any_diags = need_any_diags .or. need_height_50
        need_height_70 = MPAS_field_will_be_written('height_70hPa')
        need_height = need_height .or. need_height_70
        need_any_diags = need_any_diags .or. need_height_70
        need_height_100 = MPAS_field_will_be_written('height_100hPa')
        need_height = need_height .or. need_height_100
        need_any_diags = need_any_diags .or. need_height_100
        need_height_150 = MPAS_field_will_be_written('height_150hPa')
        need_height = need_height .or. need_height_150
        need_any_diags = need_any_diags .or. need_height_150
        need_height_200 = MPAS_field_will_be_written('height_200hPa')
        need_height = need_height .or. need_height_200
        need_any_diags = need_any_diags .or. need_height_200
        need_height_250 = MPAS_field_will_be_written('height_250hPa')
        need_height = need_height .or. need_height_250
        need_any_diags = need_any_diags .or. need_height_250
        need_height_300 = MPAS_field_will_be_written('height_300hPa')
        need_height = need_height .or. need_height_300
        need_any_diags = need_any_diags .or. need_height_300
        need_height_400 = MPAS_field_will_be_written('height_400hPa')
        need_height = need_height .or. need_height_400
        need_any_diags = need_any_diags .or. need_height_400
        need_height_500 = MPAS_field_will_be_written('height_500hPa')
        need_height = need_height .or. need_height_500
        need_any_diags = need_any_diags .or. need_height_500
        need_height_600 = MPAS_field_will_be_written('height_600hPa')
        need_height = need_height .or. need_height_600
        need_any_diags = need_any_diags .or. need_height_600
        need_height_700 = MPAS_field_will_be_written('height_700hPa')
        need_height = need_height .or. need_height_700
        need_any_diags = need_any_diags .or. need_height_700
        need_height_825 = MPAS_field_will_be_written('height_825hPa')
        need_height = need_height .or. need_height_825
        need_any_diags = need_any_diags .or. need_height_825
        need_height_850 = MPAS_field_will_be_written('height_850hPa')
        need_height = need_height .or. need_height_850
        need_any_diags = need_any_diags .or. need_height_850
        need_height_875 = MPAS_field_will_be_written('height_875hPa')
        need_height = need_height .or. need_height_875
        need_any_diags = need_any_diags .or. need_height_875
        need_height_900 = MPAS_field_will_be_written('height_900hPa')
        need_height = need_height .or. need_height_900
        need_any_diags = need_any_diags .or. need_height_900
        need_height_925 = MPAS_field_will_be_written('height_925hPa')
        need_height = need_height .or. need_height_925
        need_any_diags = need_any_diags .or. need_height_925
        need_height_950 = MPAS_field_will_be_written('height_950hPa')
        need_height = need_height .or. need_height_950
        need_any_diags = need_any_diags .or. need_height_950
        need_height_975 = MPAS_field_will_be_written('height_975hPa')
        need_height = need_height .or. need_height_975
        need_any_diags = need_any_diags .or. need_height_975
        need_height_1000 = MPAS_field_will_be_written('height_1000hPa')
        need_height = need_height .or. need_height_1000
        need_any_diags = need_any_diags .or. need_height_1000

        need_uzonal_15 = MPAS_field_will_be_written('uzonal_15hPa')
        need_uzonal = need_uzonal .or. need_uzonal_15
        need_any_diags = need_any_diags .or. need_uzonal_15
        need_uzonal_20 = MPAS_field_will_be_written('uzonal_20hPa')
        need_uzonal = need_uzonal .or. need_uzonal_20
        need_any_diags = need_any_diags .or. need_uzonal_20
        need_uzonal_30 = MPAS_field_will_be_written('uzonal_30hPa')
        need_uzonal = need_uzonal .or. need_uzonal_30
        need_any_diags = need_any_diags .or. need_uzonal_30
        need_uzonal_50 = MPAS_field_will_be_written('uzonal_50hPa')
        need_uzonal = need_uzonal .or. need_uzonal_50
        need_any_diags = need_any_diags .or. need_uzonal_50
        need_uzonal_70 = MPAS_field_will_be_written('uzonal_70hPa')
        need_uzonal = need_uzonal .or. need_uzonal_70
        need_any_diags = need_any_diags .or. need_uzonal_70
        need_uzonal_100 = MPAS_field_will_be_written('uzonal_100hPa')
        need_uzonal = need_uzonal .or. need_uzonal_100
        need_any_diags = need_any_diags .or. need_uzonal_100
        need_uzonal_150 = MPAS_field_will_be_written('uzonal_150hPa')
        need_uzonal = need_uzonal .or. need_uzonal_150
        need_any_diags = need_any_diags .or. need_uzonal_150
        need_uzonal_200 = MPAS_field_will_be_written('uzonal_200hPa')
        need_uzonal = need_uzonal .or. need_uzonal_200
        need_any_diags = need_any_diags .or. need_uzonal_200
        need_uzonal_250 = MPAS_field_will_be_written('uzonal_250hPa')
        need_uzonal = need_uzonal .or. need_uzonal_250
        need_any_diags = need_any_diags .or. need_uzonal_250
        need_uzonal_300 = MPAS_field_will_be_written('uzonal_300hPa')
        need_uzonal = need_uzonal .or. need_uzonal_300
        need_any_diags = need_any_diags .or. need_uzonal_300
        need_uzonal_400 = MPAS_field_will_be_written('uzonal_400hPa')
        need_uzonal = need_uzonal .or. need_uzonal_400
        need_any_diags = need_any_diags .or. need_uzonal_400
        need_uzonal_500 = MPAS_field_will_be_written('uzonal_500hPa')
        need_uzonal = need_uzonal .or. need_uzonal_500
        need_any_diags = need_any_diags .or. need_uzonal_500
        need_uzonal_600 = MPAS_field_will_be_written('uzonal_600hPa')
        need_uzonal = need_uzonal .or. need_uzonal_600
        need_any_diags = need_any_diags .or. need_uzonal_600
        need_uzonal_700 = MPAS_field_will_be_written('uzonal_700hPa')
        need_uzonal = need_uzonal .or. need_uzonal_700
        need_any_diags = need_any_diags .or. need_uzonal_700
        need_uzonal_825 = MPAS_field_will_be_written('uzonal_825hPa')
        need_uzonal = need_uzonal .or. need_uzonal_825
        need_any_diags = need_any_diags .or. need_uzonal_825
        need_uzonal_850 = MPAS_field_will_be_written('uzonal_850hPa')
        need_uzonal = need_uzonal .or. need_uzonal_850
        need_any_diags = need_any_diags .or. need_uzonal_850
        need_uzonal_875 = MPAS_field_will_be_written('uzonal_875hPa')
        need_uzonal = need_uzonal .or. need_uzonal_875
        need_any_diags = need_any_diags .or. need_uzonal_875
        need_uzonal_900 = MPAS_field_will_be_written('uzonal_900hPa')
        need_uzonal = need_uzonal .or. need_uzonal_900
        need_any_diags = need_any_diags .or. need_uzonal_900
        need_uzonal_925 = MPAS_field_will_be_written('uzonal_925hPa')
        need_uzonal = need_uzonal .or. need_uzonal_925
        need_any_diags = need_any_diags .or. need_uzonal_925
        need_uzonal_950 = MPAS_field_will_be_written('uzonal_950hPa')
        need_uzonal = need_uzonal .or. need_uzonal_950
        need_any_diags = need_any_diags .or. need_uzonal_950
        need_uzonal_975 = MPAS_field_will_be_written('uzonal_975hPa')
        need_uzonal = need_uzonal .or. need_uzonal_975
        need_any_diags = need_any_diags .or. need_uzonal_975
        need_uzonal_1000 = MPAS_field_will_be_written('uzonal_1000hPa')
        need_uzonal = need_uzonal .or. need_uzonal_1000
        need_any_diags = need_any_diags .or. need_uzonal_1000

        need_umeridional_15 = MPAS_field_will_be_written('umeridional_15hPa')
        need_umeridional = need_umeridional .or. need_umeridional_15
        need_any_diags = need_any_diags .or. need_umeridional_15
        need_umeridional_20 = MPAS_field_will_be_written('umeridional_20hPa')
        need_umeridional = need_umeridional .or. need_umeridional_20
        need_any_diags = need_any_diags .or. need_umeridional_20
        need_umeridional_30 = MPAS_field_will_be_written('umeridional_30hPa')
        need_umeridional = need_umeridional .or. need_umeridional_30
        need_any_diags = need_any_diags .or. need_umeridional_30
        
        need_umeridional_50 = MPAS_field_will_be_written('umeridional_50hPa')
        need_umeridional = need_umeridional .or. need_umeridional_50
        need_any_diags = need_any_diags .or. need_umeridional_50
        need_umeridional_70 = MPAS_field_will_be_written('umeridional_70hPa')
        need_umeridional = need_umeridional .or. need_umeridional_70
        need_any_diags = need_any_diags .or. need_umeridional_70
        need_umeridional_100 = MPAS_field_will_be_written('umeridional_100hPa')
        need_umeridional = need_umeridional .or. need_umeridional_100
        need_any_diags = need_any_diags .or. need_umeridional_100
        need_umeridional_150 = MPAS_field_will_be_written('umeridional_150hPa')
        need_umeridional = need_umeridional .or. need_umeridional_150
        need_any_diags = need_any_diags .or. need_umeridional_150
        need_umeridional_200 = MPAS_field_will_be_written('umeridional_200hPa')
        need_umeridional = need_umeridional .or. need_umeridional_200
        need_any_diags = need_any_diags .or. need_umeridional_200
        need_umeridional_250 = MPAS_field_will_be_written('umeridional_250hPa')
        need_umeridional = need_umeridional .or. need_umeridional_250
        need_any_diags = need_any_diags .or. need_umeridional_250
        need_umeridional_300 = MPAS_field_will_be_written('umeridional_300hPa')
        need_umeridional = need_umeridional .or. need_umeridional_300
        need_any_diags = need_any_diags .or. need_umeridional_300
        need_umeridional_400 = MPAS_field_will_be_written('umeridional_400hPa')
        need_umeridional = need_umeridional .or. need_umeridional_400
        need_any_diags = need_any_diags .or. need_umeridional_400
        need_umeridional_500 = MPAS_field_will_be_written('umeridional_500hPa')
        need_umeridional = need_umeridional .or. need_umeridional_500
        need_any_diags = need_any_diags .or. need_umeridional_500
        need_umeridional_600 = MPAS_field_will_be_written('umeridional_600hPa')
        need_umeridional = need_umeridional .or. need_umeridional_600
        need_any_diags = need_any_diags .or. need_umeridional_600
        need_umeridional_700 = MPAS_field_will_be_written('umeridional_700hPa')
        need_umeridional = need_umeridional .or. need_umeridional_700
        need_any_diags = need_any_diags .or. need_umeridional_700
        need_umeridional_825 = MPAS_field_will_be_written('umeridional_825hPa')
        need_umeridional = need_umeridional .or. need_umeridional_825
        need_any_diags = need_any_diags .or. need_umeridional_825
        need_umeridional_850 = MPAS_field_will_be_written('umeridional_850hPa')
        need_umeridional = need_umeridional .or. need_umeridional_850
        need_any_diags = need_any_diags .or. need_umeridional_850
        need_umeridional_875 = MPAS_field_will_be_written('umeridional_875hPa')
        need_umeridional = need_umeridional .or. need_umeridional_875
        need_any_diags = need_any_diags .or. need_umeridional_875
        need_umeridional_900 = MPAS_field_will_be_written('umeridional_900hPa')
        need_umeridional = need_umeridional .or. need_umeridional_900
        need_any_diags = need_any_diags .or. need_umeridional_900
        need_umeridional_925 = MPAS_field_will_be_written('umeridional_925hPa')
        need_umeridional = need_umeridional .or. need_umeridional_925
        need_any_diags = need_any_diags .or. need_umeridional_925
        need_umeridional_950 = MPAS_field_will_be_written('umeridional_950hPa')
        need_umeridional = need_umeridional .or. need_umeridional_950
        need_any_diags = need_any_diags .or. need_umeridional_950
        need_umeridional_975 = MPAS_field_will_be_written('umeridional_975hPa')
        need_umeridional = need_umeridional .or. need_umeridional_975
        need_any_diags = need_any_diags .or. need_umeridional_975
        need_umeridional_1000 = MPAS_field_will_be_written('umeridional_1000hPa')
        need_umeridional = need_umeridional .or. need_umeridional_1000
        need_any_diags = need_any_diags .or. need_umeridional_1000
        need_spechum_15 = MPAS_field_will_be_written('spechum_15hPa')
        need_spechum = need_spechum .or. need_spechum_15
        need_any_diags = need_any_diags .or. need_spechum_15
        need_spechum_20 = MPAS_field_will_be_written('spechum_20hPa')
        need_spechum = need_spechum .or. need_spechum_20
        need_any_diags = need_any_diags .or. need_spechum_20
        need_spechum_30 = MPAS_field_will_be_written('spechum_30hPa')
        need_spechum = need_spechum .or. need_spechum_30
        need_any_diags = need_any_diags .or. need_spechum_30
        need_spechum_50 = MPAS_field_will_be_written('spechum_50hPa')
        need_spechum = need_spechum .or. need_spechum_50
        need_any_diags = need_any_diags .or. need_spechum_50
        need_spechum_70 = MPAS_field_will_be_written('spechum_70hPa')
        need_spechum = need_spechum .or. need_spechum_70
        need_any_diags = need_any_diags .or. need_spechum_70
        need_spechum_100 = MPAS_field_will_be_written('spechum_100hPa')
        need_spechum = need_spechum .or. need_spechum_100
        need_any_diags = need_any_diags .or. need_spechum_100
        need_spechum_150 = MPAS_field_will_be_written('spechum_150hPa')
        need_spechum = need_spechum .or. need_spechum_150
        need_any_diags = need_any_diags .or. need_spechum_150
        need_spechum_200 = MPAS_field_will_be_written('spechum_200hPa')
        need_spechum = need_spechum .or. need_spechum_200
        need_any_diags = need_any_diags .or. need_spechum_200
        need_spechum_250 = MPAS_field_will_be_written('spechum_250hPa')
        need_spechum = need_spechum .or. need_spechum_250
        need_any_diags = need_any_diags .or. need_spechum_250
        need_spechum_300 = MPAS_field_will_be_written('spechum_300hPa')
        need_spechum = need_spechum .or. need_spechum_300
        need_any_diags = need_any_diags .or. need_spechum_300
        need_spechum_400 = MPAS_field_will_be_written('spechum_400hPa')
        need_spechum = need_spechum .or. need_spechum_400
        need_any_diags = need_any_diags .or. need_spechum_400
        need_spechum_500 = MPAS_field_will_be_written('spechum_500hPa')
        need_spechum = need_spechum .or. need_spechum_500
        need_any_diags = need_any_diags .or. need_spechum_500
        need_spechum_600 = MPAS_field_will_be_written('spechum_600hPa')
        need_spechum = need_spechum .or. need_spechum_600
        need_any_diags = need_any_diags .or. need_spechum_600
        need_spechum_700 = MPAS_field_will_be_written('spechum_700hPa')
        need_spechum = need_spechum .or. need_spechum_700
        need_any_diags = need_any_diags .or. need_spechum_700
        need_spechum_825 = MPAS_field_will_be_written('spechum_825hPa')
        need_spechum = need_spechum .or. need_spechum_825
        need_any_diags = need_any_diags .or. need_spechum_825
        need_spechum_850 = MPAS_field_will_be_written('spechum_850hPa')
        need_spechum = need_spechum .or. need_spechum_850
        need_any_diags = need_any_diags .or. need_spechum_850
        need_spechum_875 = MPAS_field_will_be_written('spechum_875hPa')
        need_spechum = need_spechum .or. need_spechum_875
        need_any_diags = need_any_diags .or. need_spechum_875
        need_spechum_900 = MPAS_field_will_be_written('spechum_900hPa')
        need_spechum = need_spechum .or. need_spechum_900
        need_any_diags = need_any_diags .or. need_spechum_900
        need_spechum_925 = MPAS_field_will_be_written('spechum_925hPa')
        need_spechum = need_spechum .or. need_spechum_925
        need_any_diags = need_any_diags .or. need_spechum_925
        need_spechum_950 = MPAS_field_will_be_written('spechum_950hPa')
        need_spechum = need_spechum .or. need_spechum_950
        need_any_diags = need_any_diags .or. need_spechum_950
        need_spechum_975 = MPAS_field_will_be_written('spechum_975hPa')
        need_spechum = need_spechum .or. need_spechum_975
        need_any_diags = need_any_diags .or. need_spechum_975
        need_spechum_1000 = MPAS_field_will_be_written('spechum_1000hPa')
        need_spechum = need_spechum .or. need_spechum_1000
        need_any_diags = need_any_diags .or. need_spechum_1000

        need_w_15 = MPAS_field_will_be_written('w_15hPa')
        need_w = need_w .or. need_w_15
        need_any_diags = need_any_diags .or. need_w_15
        need_w_20 = MPAS_field_will_be_written('w_20hPa')
        need_w = need_w .or. need_w_20
        need_any_diags = need_any_diags .or. need_w_20
        need_w_30 = MPAS_field_will_be_written('w_30hPa')
        need_w = need_w .or. need_w_30
        need_any_diags = need_any_diags .or. need_w_30
        need_w_50 = MPAS_field_will_be_written('w_50hPa')
        need_w = need_w .or. need_w_50
        need_any_diags = need_any_diags .or. need_w_50
        need_w_70 = MPAS_field_will_be_written('w_70hPa')
        need_w = need_w .or. need_w_70
        need_any_diags = need_any_diags .or. need_w_70
        need_w_100 = MPAS_field_will_be_written('w_100hPa')
        need_w = need_w .or. need_w_100
        need_any_diags = need_any_diags .or. need_w_100
        need_w_150 = MPAS_field_will_be_written('w_150hPa')
        need_w = need_w .or. need_w_150
        need_any_diags = need_any_diags .or. need_w_150
        need_w_200 = MPAS_field_will_be_written('w_200hPa')
        need_w = need_w .or. need_w_200
        need_any_diags = need_any_diags .or. need_w_200
        need_w_250 = MPAS_field_will_be_written('w_250hPa')
        need_w = need_w .or. need_w_250
        need_any_diags = need_any_diags .or. need_w_250
        need_w_300 = MPAS_field_will_be_written('w_300hPa')
        need_w = need_w .or. need_w_300
        need_any_diags = need_any_diags .or. need_w_300
        need_w_400 = MPAS_field_will_be_written('w_400hPa')
        need_w = need_w .or. need_w_400
        need_any_diags = need_any_diags .or. need_w_400
        need_w_500 = MPAS_field_will_be_written('w_500hPa')
        need_w = need_w .or. need_w_500
        need_any_diags = need_any_diags .or. need_w_500
        need_w_600 = MPAS_field_will_be_written('w_600hPa')
        need_w = need_w .or. need_w_600
        need_any_diags = need_any_diags .or. need_w_600
        need_w_700 = MPAS_field_will_be_written('w_700hPa')
        need_w = need_w .or. need_w_700
        need_any_diags = need_any_diags .or. need_w_700
        need_w_825 = MPAS_field_will_be_written('w_825hPa')
        need_w = need_w .or. need_w_825
        need_any_diags = need_any_diags .or. need_w_825
        need_w_850 = MPAS_field_will_be_written('w_850hPa')
        need_w = need_w .or. need_w_850
        need_any_diags = need_any_diags .or. need_w_850
        need_w_875 = MPAS_field_will_be_written('w_875hPa')
        need_w = need_w .or. need_w_875
        need_any_diags = need_any_diags .or. need_w_875
        need_w_900 = MPAS_field_will_be_written('w_900hPa')
        need_w = need_w .or. need_w_900
        need_any_diags = need_any_diags .or. need_w_900
        need_w_925 = MPAS_field_will_be_written('w_925hPa')
        need_w = need_w .or. need_w_925
        need_any_diags = need_any_diags .or. need_w_925
        need_w_950 = MPAS_field_will_be_written('w_950hPa')
        need_w = need_w .or. need_w_950
        need_any_diags = need_any_diags .or. need_w_950
        need_w_975 = MPAS_field_will_be_written('w_975hPa')
        need_w = need_w .or. need_w_975
        need_any_diags = need_any_diags .or. need_w_975
        need_w_1000 = MPAS_field_will_be_written('w_1000hPa')
        need_w = need_w .or. need_w_1000
        need_any_diags = need_any_diags .or. need_w_1000

        need_geoph_15 = MPAS_field_will_be_written('geoph_15hPa')
        need_geoph = need_geoph .or. need_geoph_15
        need_any_diags = need_any_diags .or. need_geoph_15
        need_geoph_20 = MPAS_field_will_be_written('geoph_20hPa')
        need_geoph = need_geoph .or. need_geoph_20
        need_any_diags = need_any_diags .or. need_geoph_20
        need_geoph_30 = MPAS_field_will_be_written('geoph_30hPa')
        need_geoph = need_geoph .or. need_geoph_30
        need_any_diags = need_any_diags .or. need_geoph_30
        need_geoph_50 = MPAS_field_will_be_written('geoph_50hPa')
        need_geoph = need_geoph .or. need_geoph_50
        need_any_diags = need_any_diags .or. need_geoph_50
        need_geoph_70 = MPAS_field_will_be_written('geoph_70hPa')
        need_geoph = need_geoph .or. need_geoph_70
        need_any_diags = need_any_diags .or. need_geoph_70
        need_geoph_100 = MPAS_field_will_be_written('geoph_100hPa')
        need_geoph = need_geoph .or. need_geoph_100
        need_any_diags = need_any_diags .or. need_geoph_100
        need_geoph_150 = MPAS_field_will_be_written('geoph_150hPa')
        need_geoph = need_geoph .or. need_geoph_150
        need_any_diags = need_any_diags .or. need_geoph_150
        need_geoph_200 = MPAS_field_will_be_written('geoph_200hPa')
        need_geoph = need_geoph .or. need_geoph_200
        need_any_diags = need_any_diags .or. need_geoph_200
        need_geoph_250 = MPAS_field_will_be_written('geoph_250hPa')
        need_geoph = need_geoph .or. need_geoph_250
        need_any_diags = need_any_diags .or. need_geoph_250
        need_geoph_300 = MPAS_field_will_be_written('geoph_300hPa')
        need_geoph = need_geoph .or. need_geoph_300
        need_any_diags = need_any_diags .or. need_geoph_300
        need_geoph_400 = MPAS_field_will_be_written('geoph_400hPa')
        need_geoph = need_geoph .or. need_geoph_400
        need_any_diags = need_any_diags .or. need_geoph_400
        need_geoph_500 = MPAS_field_will_be_written('geoph_500hPa')
        need_geoph = need_geoph .or. need_geoph_500
        need_any_diags = need_any_diags .or. need_geoph_500
        need_geoph_600 = MPAS_field_will_be_written('geoph_600hPa')
        need_geoph = need_geoph .or. need_geoph_600
        need_any_diags = need_any_diags .or. need_geoph_600
        need_geoph_700 = MPAS_field_will_be_written('geoph_700hPa')
        need_geoph = need_geoph .or. need_geoph_700
        need_any_diags = need_any_diags .or. need_geoph_700
        need_geoph_825 = MPAS_field_will_be_written('geoph_825hPa')
        need_geoph = need_geoph .or. need_geoph_825
        need_any_diags = need_any_diags .or. need_geoph_825
        need_geoph_850 = MPAS_field_will_be_written('geoph_850hPa')
        need_geoph = need_geoph .or. need_geoph_850
        need_any_diags = need_any_diags .or. need_geoph_850
        need_geoph_875 = MPAS_field_will_be_written('geoph_875hPa')
        need_geoph = need_geoph .or. need_geoph_875
        need_any_diags = need_any_diags .or. need_geoph_875
        need_geoph_900 = MPAS_field_will_be_written('geoph_900hPa')
        need_geoph = need_geoph .or. need_geoph_900
        need_any_diags = need_any_diags .or. need_geoph_900
        need_geoph_925 = MPAS_field_will_be_written('geoph_925hPa')
        need_geoph = need_geoph .or. need_geoph_925
        need_any_diags = need_any_diags .or. need_geoph_925
        need_geoph_950 = MPAS_field_will_be_written('geoph_950hPa')
        need_geoph = need_geoph .or. need_geoph_950
        need_any_diags = need_any_diags .or. need_geoph_950
        need_geoph_975 = MPAS_field_will_be_written('geoph_975hPa')
        need_geoph = need_geoph .or. need_geoph_975
        need_any_diags = need_any_diags .or. need_geoph_975
        need_geoph_1000 = MPAS_field_will_be_written('geoph_1000hPa')
        need_geoph = need_geoph .or. need_geoph_1000
        need_any_diags = need_any_diags .or. need_geoph_1000

        need_rqvblten_15 = MPAS_field_will_be_written('rqvblten_15hPa')
        need_rqvblten = need_rqvblten .or. need_rqvblten_15
        need_any_diags = need_any_diags .or. need_rqvblten_15
        need_rqvblten_20 = MPAS_field_will_be_written('rqvblten_20hPa')
        need_rqvblten = need_rqvblten .or. need_rqvblten_20
        need_any_diags = need_any_diags .or. need_rqvblten_20
        need_rqvblten_30 = MPAS_field_will_be_written('rqvblten_30hPa')
        need_rqvblten = need_rqvblten .or. need_rqvblten_30
        need_any_diags = need_any_diags .or. need_rqvblten_30
        need_rqvblten_50 = MPAS_field_will_be_written('rqvblten_50hPa')
        need_rqvblten = need_rqvblten .or. need_rqvblten_50
        need_any_diags = need_any_diags .or. need_rqvblten_50
        need_rqvblten_70 = MPAS_field_will_be_written('rqvblten_70hPa')
        need_rqvblten = need_rqvblten .or. need_rqvblten_70
        need_any_diags = need_any_diags .or. need_rqvblten_70
        need_rqvblten_100 = MPAS_field_will_be_written('rqvblten_100hPa')
        need_rqvblten = need_rqvblten .or. need_rqvblten_100
        need_any_diags = need_any_diags .or. need_rqvblten_100
        need_rqvblten_150 = MPAS_field_will_be_written('rqvblten_150hPa')
        need_rqvblten = need_rqvblten .or. need_rqvblten_150
        need_any_diags = need_any_diags .or. need_rqvblten_150
        need_rqvblten_200 = MPAS_field_will_be_written('rqvblten_200hPa')
        need_rqvblten = need_rqvblten .or. need_rqvblten_200
        need_any_diags = need_any_diags .or. need_rqvblten_200
        need_rqvblten_250 = MPAS_field_will_be_written('rqvblten_250hPa')
        need_rqvblten = need_rqvblten .or. need_rqvblten_250
        need_any_diags = need_any_diags .or. need_rqvblten_250
        need_rqvblten_300 = MPAS_field_will_be_written('rqvblten_300hPa')
        need_rqvblten = need_rqvblten .or. need_rqvblten_300
        need_any_diags = need_any_diags .or. need_rqvblten_300
        need_rqvblten_400 = MPAS_field_will_be_written('rqvblten_400hPa')
        need_rqvblten = need_rqvblten .or. need_rqvblten_400
        need_any_diags = need_any_diags .or. need_rqvblten_400
        need_rqvblten_500 = MPAS_field_will_be_written('rqvblten_500hPa')
        need_rqvblten = need_rqvblten .or. need_rqvblten_500
        need_any_diags = need_any_diags .or. need_rqvblten_500
        need_rqvblten_600 = MPAS_field_will_be_written('rqvblten_600hPa')
        need_rqvblten = need_rqvblten .or. need_rqvblten_600
        need_any_diags = need_any_diags .or. need_rqvblten_600
        need_rqvblten_700 = MPAS_field_will_be_written('rqvblten_700hPa')
        need_rqvblten = need_rqvblten .or. need_rqvblten_700
        need_any_diags = need_any_diags .or. need_rqvblten_700
        need_rqvblten_825 = MPAS_field_will_be_written('rqvblten_825hPa')
        need_rqvblten = need_rqvblten .or. need_rqvblten_825
        need_any_diags = need_any_diags .or. need_rqvblten_825
        need_rqvblten_850 = MPAS_field_will_be_written('rqvblten_850hPa')
        need_rqvblten = need_rqvblten .or. need_rqvblten_850
        need_any_diags = need_any_diags .or. need_rqvblten_850
        need_rqvblten_875 = MPAS_field_will_be_written('rqvblten_875hPa')
        need_rqvblten = need_rqvblten .or. need_rqvblten_875
        need_any_diags = need_any_diags .or. need_rqvblten_875
        need_rqvblten_900 = MPAS_field_will_be_written('rqvblten_900hPa')
        need_rqvblten = need_rqvblten .or. need_rqvblten_900
        need_any_diags = need_any_diags .or. need_rqvblten_900
        need_rqvblten_925 = MPAS_field_will_be_written('rqvblten_925hPa')
        need_rqvblten = need_rqvblten .or. need_rqvblten_925
        need_any_diags = need_any_diags .or. need_rqvblten_925
        need_rqvblten_950 = MPAS_field_will_be_written('rqvblten_950hPa')
        need_rqvblten = need_rqvblten .or. need_rqvblten_950
        need_any_diags = need_any_diags .or. need_rqvblten_950
        need_rqvblten_975 = MPAS_field_will_be_written('rqvblten_975hPa')
        need_rqvblten = need_rqvblten .or. need_rqvblten_975
        need_any_diags = need_any_diags .or. need_rqvblten_975
        need_rqvblten_1000 = MPAS_field_will_be_written('rqvblten_1000hPa')
        need_rqvblten = need_rqvblten .or. need_rqvblten_1000
        need_any_diags = need_any_diags .or. need_rqvblten_1000

        need_rqvcuten_15 = MPAS_field_will_be_written('rqvcuten_15hPa')
        need_rqvcuten = need_rqvcuten .or. need_rqvcuten_15
        need_any_diags = need_any_diags .or. need_rqvcuten_15
        need_rqvcuten_20 = MPAS_field_will_be_written('rqvcuten_20hPa')
        need_rqvcuten = need_rqvcuten .or. need_rqvcuten_20
        need_any_diags = need_any_diags .or. need_rqvcuten_20
        need_rqvcuten_30 = MPAS_field_will_be_written('rqvcuten_30hPa')
        need_rqvcuten = need_rqvcuten .or. need_rqvcuten_30
        need_any_diags = need_any_diags .or. need_rqvcuten_30
        need_rqvcuten_50 = MPAS_field_will_be_written('rqvcuten_50hPa')
        need_rqvcuten = need_rqvcuten .or. need_rqvcuten_50
        need_any_diags = need_any_diags .or. need_rqvcuten_50
        need_rqvcuten_70 = MPAS_field_will_be_written('rqvcuten_70hPa')
        need_rqvcuten = need_rqvcuten .or. need_rqvcuten_70
        need_any_diags = need_any_diags .or. need_rqvcuten_70
        need_rqvcuten_100 = MPAS_field_will_be_written('rqvcuten_100hPa')
        need_rqvcuten = need_rqvcuten .or. need_rqvcuten_100
        need_any_diags = need_any_diags .or. need_rqvcuten_100
        need_rqvcuten_150 = MPAS_field_will_be_written('rqvcuten_150hPa')
        need_rqvcuten = need_rqvcuten .or. need_rqvcuten_150
        need_any_diags = need_any_diags .or. need_rqvcuten_150
        need_rqvcuten_200 = MPAS_field_will_be_written('rqvcuten_200hPa')
        need_rqvcuten = need_rqvcuten .or. need_rqvcuten_200
        need_any_diags = need_any_diags .or. need_rqvcuten_200
        need_rqvcuten_250 = MPAS_field_will_be_written('rqvcuten_250hPa')
        need_rqvcuten = need_rqvcuten .or. need_rqvcuten_250
        need_any_diags = need_any_diags .or. need_rqvcuten_250
        need_rqvcuten_300 = MPAS_field_will_be_written('rqvcuten_300hPa')
        need_rqvcuten = need_rqvcuten .or. need_rqvcuten_300
        need_any_diags = need_any_diags .or. need_rqvcuten_300
        need_rqvcuten_400 = MPAS_field_will_be_written('rqvcuten_400hPa')
        need_rqvcuten = need_rqvcuten .or. need_rqvcuten_400
        need_any_diags = need_any_diags .or. need_rqvcuten_400
        need_rqvcuten_500 = MPAS_field_will_be_written('rqvcuten_500hPa')
        need_rqvcuten = need_rqvcuten .or. need_rqvcuten_500
        need_any_diags = need_any_diags .or. need_rqvcuten_500
        need_rqvcuten_600 = MPAS_field_will_be_written('rqvcuten_600hPa')
        need_rqvcuten = need_rqvcuten .or. need_rqvcuten_600
        need_any_diags = need_any_diags .or. need_rqvcuten_600
        need_rqvcuten_700 = MPAS_field_will_be_written('rqvcuten_700hPa')
        need_rqvcuten = need_rqvcuten .or. need_rqvcuten_700
        need_any_diags = need_any_diags .or. need_rqvcuten_700
        need_rqvcuten_825 = MPAS_field_will_be_written('rqvcuten_825hPa')
        need_rqvcuten = need_rqvcuten .or. need_rqvcuten_825
        need_any_diags = need_any_diags .or. need_rqvcuten_825
        need_rqvcuten_850 = MPAS_field_will_be_written('rqvcuten_850hPa')
        need_rqvcuten = need_rqvcuten .or. need_rqvcuten_850
        need_any_diags = need_any_diags .or. need_rqvcuten_850
        need_rqvcuten_875 = MPAS_field_will_be_written('rqvcuten_875hPa')
        need_rqvcuten = need_rqvcuten .or. need_rqvcuten_875
        need_any_diags = need_any_diags .or. need_rqvcuten_875
        need_rqvcuten_900 = MPAS_field_will_be_written('rqvcuten_900hPa')
        need_rqvcuten = need_rqvcuten .or. need_rqvcuten_900
        need_any_diags = need_any_diags .or. need_rqvcuten_900
        need_rqvcuten_925 = MPAS_field_will_be_written('rqvcuten_925hPa')
        need_rqvcuten = need_rqvcuten .or. need_rqvcuten_925
        need_any_diags = need_any_diags .or. need_rqvcuten_925
        need_rqvcuten_950 = MPAS_field_will_be_written('rqvcuten_950hPa')
        need_rqvcuten = need_rqvcuten .or. need_rqvcuten_950
        need_any_diags = need_any_diags .or. need_rqvcuten_950
        need_rqvcuten_975 = MPAS_field_will_be_written('rqvcuten_975hPa')
        need_rqvcuten = need_rqvcuten .or. need_rqvcuten_975
        need_any_diags = need_any_diags .or. need_rqvcuten_975
        need_rqvcuten_1000 = MPAS_field_will_be_written('rqvcuten_1000hPa')
        need_rqvcuten = need_rqvcuten .or. need_rqvcuten_1000
        need_any_diags = need_any_diags .or. need_rqvcuten_1000

        need_rthratenlw_15 = MPAS_field_will_be_written('rthratenlw_15hPa')
        need_rthratenlw = need_rthratenlw .or. need_rthratenlw_15
        need_any_diags = need_any_diags .or. need_rthratenlw_15
        need_rthratenlw_20 = MPAS_field_will_be_written('rthratenlw_20hPa')
        need_rthratenlw = need_rthratenlw .or. need_rthratenlw_20
        need_any_diags = need_any_diags .or. need_rthratenlw_20
        need_rthratenlw_30 = MPAS_field_will_be_written('rthratenlw_30hPa')
        need_rthratenlw = need_rthratenlw .or. need_rthratenlw_30
        need_any_diags = need_any_diags .or. need_rthratenlw_30
        need_rthratenlw_50 = MPAS_field_will_be_written('rthratenlw_50hPa')
        need_rthratenlw = need_rthratenlw .or. need_rthratenlw_50
        need_any_diags = need_any_diags .or. need_rthratenlw_50
        need_rthratenlw_70 = MPAS_field_will_be_written('rthratenlw_70hPa')
        need_rthratenlw = need_rthratenlw .or. need_rthratenlw_70
        need_any_diags = need_any_diags .or. need_rthratenlw_70
        need_rthratenlw_100 = MPAS_field_will_be_written('rthratenlw_100hPa')
        need_rthratenlw = need_rthratenlw .or. need_rthratenlw_100
        need_any_diags = need_any_diags .or. need_rthratenlw_100
        need_rthratenlw_150 = MPAS_field_will_be_written('rthratenlw_150hPa')
        need_rthratenlw = need_rthratenlw .or. need_rthratenlw_150
        need_any_diags = need_any_diags .or. need_rthratenlw_150
        need_rthratenlw_200 = MPAS_field_will_be_written('rthratenlw_200hPa')
        need_rthratenlw = need_rthratenlw .or. need_rthratenlw_200
        need_any_diags = need_any_diags .or. need_rthratenlw_200
        need_rthratenlw_250 = MPAS_field_will_be_written('rthratenlw_250hPa')
        need_rthratenlw = need_rthratenlw .or. need_rthratenlw_250
        need_any_diags = need_any_diags .or. need_rthratenlw_250
        need_rthratenlw_300 = MPAS_field_will_be_written('rthratenlw_300hPa')
        need_rthratenlw = need_rthratenlw .or. need_rthratenlw_300
        need_any_diags = need_any_diags .or. need_rthratenlw_300
        need_rthratenlw_400 = MPAS_field_will_be_written('rthratenlw_400hPa')
        need_rthratenlw = need_rthratenlw .or. need_rthratenlw_400
        need_any_diags = need_any_diags .or. need_rthratenlw_400
        need_rthratenlw_500 = MPAS_field_will_be_written('rthratenlw_500hPa')
        need_rthratenlw = need_rthratenlw .or. need_rthratenlw_500
        need_any_diags = need_any_diags .or. need_rthratenlw_500
        need_rthratenlw_600 = MPAS_field_will_be_written('rthratenlw_600hPa')
        need_rthratenlw = need_rthratenlw .or. need_rthratenlw_600
        need_any_diags = need_any_diags .or. need_rthratenlw_600
        need_rthratenlw_700 = MPAS_field_will_be_written('rthratenlw_700hPa')
        need_rthratenlw = need_rthratenlw .or. need_rthratenlw_700
        need_any_diags = need_any_diags .or. need_rthratenlw_700
        need_rthratenlw_825 = MPAS_field_will_be_written('rthratenlw_825hPa')
        need_rthratenlw = need_rthratenlw .or. need_rthratenlw_825
        need_any_diags = need_any_diags .or. need_rthratenlw_825
        need_rthratenlw_850 = MPAS_field_will_be_written('rthratenlw_850hPa')
        need_rthratenlw = need_rthratenlw .or. need_rthratenlw_850
        need_any_diags = need_any_diags .or. need_rthratenlw_850
        need_rthratenlw_875 = MPAS_field_will_be_written('rthratenlw_875hPa')
        need_rthratenlw = need_rthratenlw .or. need_rthratenlw_875
        need_any_diags = need_any_diags .or. need_rthratenlw_875
        need_rthratenlw_900 = MPAS_field_will_be_written('rthratenlw_900hPa')
        need_rthratenlw = need_rthratenlw .or. need_rthratenlw_900
        need_any_diags = need_any_diags .or. need_rthratenlw_900
        need_rthratenlw_925 = MPAS_field_will_be_written('rthratenlw_925hPa')
        need_rthratenlw = need_rthratenlw .or. need_rthratenlw_925
        need_any_diags = need_any_diags .or. need_rthratenlw_925
        need_rthratenlw_950 = MPAS_field_will_be_written('rthratenlw_950hPa')
        need_rthratenlw = need_rthratenlw .or. need_rthratenlw_950
        need_any_diags = need_any_diags .or. need_rthratenlw_950
        need_rthratenlw_975 = MPAS_field_will_be_written('rthratenlw_975hPa')
        need_rthratenlw = need_rthratenlw .or. need_rthratenlw_975
        need_any_diags = need_any_diags .or. need_rthratenlw_975
        need_rthratenlw_1000 = MPAS_field_will_be_written('rthratenlw_1000hPa')
        need_rthratenlw = need_rthratenlw .or. need_rthratenlw_1000
        need_any_diags = need_any_diags .or. need_rthratenlw_1000

        need_rucuten_15 = MPAS_field_will_be_written('rucuten_15hPa')
        need_rucuten = need_rucuten .or. need_rucuten_15
        need_any_diags = need_any_diags .or. need_rucuten_15
        need_rucuten_20 = MPAS_field_will_be_written('rucuten_20hPa')
        need_rucuten = need_rucuten .or. need_rucuten_20
        need_any_diags = need_any_diags .or. need_rucuten_20
        need_rucuten_30 = MPAS_field_will_be_written('rucuten_30hPa')
        need_rucuten = need_rucuten .or. need_rucuten_30
        need_any_diags = need_any_diags .or. need_rucuten_30
        need_rucuten_50 = MPAS_field_will_be_written('rucuten_50hPa')
        need_rucuten = need_rucuten .or. need_rucuten_50
        need_any_diags = need_any_diags .or. need_rucuten_50
        need_rucuten_70 = MPAS_field_will_be_written('rucuten_70hPa')
        need_rucuten = need_rucuten .or. need_rucuten_70
        need_any_diags = need_any_diags .or. need_rucuten_70
        need_rucuten_100 = MPAS_field_will_be_written('rucuten_100hPa')
        need_rucuten = need_rucuten .or. need_rucuten_100
        need_any_diags = need_any_diags .or. need_rucuten_100
        need_rucuten_150 = MPAS_field_will_be_written('rucuten_150hPa')
        need_rucuten = need_rucuten .or. need_rucuten_150
        need_any_diags = need_any_diags .or. need_rucuten_150
        need_rucuten_200 = MPAS_field_will_be_written('rucuten_200hPa')
        need_rucuten = need_rucuten .or. need_rucuten_200
        need_any_diags = need_any_diags .or. need_rucuten_200
        need_rucuten_250 = MPAS_field_will_be_written('rucuten_250hPa')
        need_rucuten = need_rucuten .or. need_rucuten_250
        need_any_diags = need_any_diags .or. need_rucuten_250
        need_rucuten_300 = MPAS_field_will_be_written('rucuten_300hPa')
        need_rucuten = need_rucuten .or. need_rucuten_300
        need_any_diags = need_any_diags .or. need_rucuten_300
        need_rucuten_400 = MPAS_field_will_be_written('rucuten_400hPa')
        need_rucuten = need_rucuten .or. need_rucuten_400
        need_any_diags = need_any_diags .or. need_rucuten_400
        need_rucuten_500 = MPAS_field_will_be_written('rucuten_500hPa')
        need_rucuten = need_rucuten .or. need_rucuten_500
        need_any_diags = need_any_diags .or. need_rucuten_500
        need_rucuten_600 = MPAS_field_will_be_written('rucuten_600hPa')
        need_rucuten = need_rucuten .or. need_rucuten_600
        need_any_diags = need_any_diags .or. need_rucuten_600
        need_rucuten_700 = MPAS_field_will_be_written('rucuten_700hPa')
        need_rucuten = need_rucuten .or. need_rucuten_700
        need_any_diags = need_any_diags .or. need_rucuten_700
        need_rucuten_825 = MPAS_field_will_be_written('rucuten_825hPa')
        need_rucuten = need_rucuten .or. need_rucuten_825
        need_any_diags = need_any_diags .or. need_rucuten_825
        need_rucuten_850 = MPAS_field_will_be_written('rucuten_850hPa')
        need_rucuten = need_rucuten .or. need_rucuten_850
        need_any_diags = need_any_diags .or. need_rucuten_850
        need_rucuten_875 = MPAS_field_will_be_written('rucuten_875hPa')
        need_rucuten = need_rucuten .or. need_rucuten_875
        need_any_diags = need_any_diags .or. need_rucuten_875
        need_rucuten_900 = MPAS_field_will_be_written('rucuten_900hPa')
        need_rucuten = need_rucuten .or. need_rucuten_900
        need_any_diags = need_any_diags .or. need_rucuten_900
        need_rucuten_925 = MPAS_field_will_be_written('rucuten_925hPa')
        need_rucuten = need_rucuten .or. need_rucuten_925
        need_any_diags = need_any_diags .or. need_rucuten_925
        need_rucuten_950 = MPAS_field_will_be_written('rucuten_950hPa')
        need_rucuten = need_rucuten .or. need_rucuten_950
        need_any_diags = need_any_diags .or. need_rucuten_950
        need_rucuten_975 = MPAS_field_will_be_written('rucuten_975hPa')
        need_rucuten = need_rucuten .or. need_rucuten_975
        need_any_diags = need_any_diags .or. need_rucuten_975
        need_rucuten_1000 = MPAS_field_will_be_written('rucuten_1000hPa')
        need_rucuten = need_rucuten .or. need_rucuten_1000
        need_any_diags = need_any_diags .or. need_rucuten_1000

        need_kzq_15 = MPAS_field_will_be_written('kzq_15hPa')
        need_kzq = need_kzq .or. need_kzq_15
        need_any_diags = need_any_diags .or. need_kzq_15
        need_kzq_20 = MPAS_field_will_be_written('kzq_20hPa')
        need_kzq = need_kzq .or. need_kzq_20
        need_any_diags = need_any_diags .or. need_kzq_20
        need_kzq_30 = MPAS_field_will_be_written('kzq_30hPa')
        need_kzq = need_kzq .or. need_kzq_30
        need_any_diags = need_any_diags .or. need_kzq_30
        need_kzq_50 = MPAS_field_will_be_written('kzq_50hPa')
        need_kzq = need_kzq .or. need_kzq_50
        need_any_diags = need_any_diags .or. need_kzq_50
        need_kzq_70 = MPAS_field_will_be_written('kzq_70hPa')
        need_kzq = need_kzq .or. need_kzq_70
        need_any_diags = need_any_diags .or. need_kzq_70
        need_kzq_100 = MPAS_field_will_be_written('kzq_100hPa')
        need_kzq = need_kzq .or. need_kzq_100
        need_any_diags = need_any_diags .or. need_kzq_100
        need_kzq_150 = MPAS_field_will_be_written('kzq_150hPa')
        need_kzq = need_kzq .or. need_kzq_150
        need_any_diags = need_any_diags .or. need_kzq_150
        need_kzq_200 = MPAS_field_will_be_written('kzq_200hPa')
        need_kzq = need_kzq .or. need_kzq_200
        need_any_diags = need_any_diags .or. need_kzq_200
        need_kzq_250 = MPAS_field_will_be_written('kzq_250hPa')
        need_kzq = need_kzq .or. need_kzq_250
        need_any_diags = need_any_diags .or. need_kzq_250
        need_kzq_300 = MPAS_field_will_be_written('kzq_300hPa')
        need_kzq = need_kzq .or. need_kzq_300
        need_any_diags = need_any_diags .or. need_kzq_300
        need_kzq_400 = MPAS_field_will_be_written('kzq_400hPa')
        need_kzq = need_kzq .or. need_kzq_400
        need_any_diags = need_any_diags .or. need_kzq_400
        need_kzq_500 = MPAS_field_will_be_written('kzq_500hPa')
        need_kzq = need_kzq .or. need_kzq_500
        need_any_diags = need_any_diags .or. need_kzq_500
        need_kzq_600 = MPAS_field_will_be_written('kzq_600hPa')
        need_kzq = need_kzq .or. need_kzq_600
        need_any_diags = need_any_diags .or. need_kzq_600
        need_kzq_700 = MPAS_field_will_be_written('kzq_700hPa')
        need_kzq = need_kzq .or. need_kzq_700
        need_any_diags = need_any_diags .or. need_kzq_700
        need_kzq_825 = MPAS_field_will_be_written('kzq_825hPa')
        need_kzq = need_kzq .or. need_kzq_825
        need_any_diags = need_any_diags .or. need_kzq_825
        need_kzq_850 = MPAS_field_will_be_written('kzq_850hPa')
        need_kzq = need_kzq .or. need_kzq_850
        need_any_diags = need_any_diags .or. need_kzq_850
        need_kzq_875 = MPAS_field_will_be_written('kzq_875hPa')
        need_kzq = need_kzq .or. need_kzq_875
        need_any_diags = need_any_diags .or. need_kzq_875
        need_kzq_900 = MPAS_field_will_be_written('kzq_900hPa')
        need_kzq = need_kzq .or. need_kzq_900
        need_any_diags = need_any_diags .or. need_kzq_900
        need_kzq_925 = MPAS_field_will_be_written('kzq_925hPa')
        need_kzq = need_kzq .or. need_kzq_925
        need_any_diags = need_any_diags .or. need_kzq_925
        need_kzq_950 = MPAS_field_will_be_written('kzq_950hPa')
        need_kzq = need_kzq .or. need_kzq_950
        need_any_diags = need_any_diags .or. need_kzq_950
        need_kzq_975 = MPAS_field_will_be_written('kzq_975hPa')
        need_kzq = need_kzq .or. need_kzq_975
        need_any_diags = need_any_diags .or. need_kzq_975
        need_kzq_1000 = MPAS_field_will_be_written('kzq_1000hPa')
        need_kzq = need_kzq .or. need_kzq_1000
        need_any_diags = need_any_diags .or. need_kzq_1000

        need_kzh_15 = MPAS_field_will_be_written('kzh_15hPa')
        need_kzh = need_kzh .or. need_kzh_15
        need_any_diags = need_any_diags .or. need_kzh_15
        need_kzh_20 = MPAS_field_will_be_written('kzh_20hPa')
        need_kzh = need_kzh .or. need_kzh_20
        need_any_diags = need_any_diags .or. need_kzh_20
        need_kzh_30 = MPAS_field_will_be_written('kzh_30hPa')
        need_kzh = need_kzh .or. need_kzh_30
        need_any_diags = need_any_diags .or. need_kzh_30
        need_kzh_50 = MPAS_field_will_be_written('kzh_50hPa')
        need_kzh = need_kzh .or. need_kzh_50
        need_any_diags = need_any_diags .or. need_kzh_50
        need_kzh_70 = MPAS_field_will_be_written('kzh_70hPa')
        need_kzh = need_kzh .or. need_kzh_70
        need_any_diags = need_any_diags .or. need_kzh_70
        need_kzh_100 = MPAS_field_will_be_written('kzh_100hPa')
        need_kzh = need_kzh .or. need_kzh_100
        need_any_diags = need_any_diags .or. need_kzh_100
        need_kzh_150 = MPAS_field_will_be_written('kzh_150hPa')
        need_kzh = need_kzh .or. need_kzh_150
        need_any_diags = need_any_diags .or. need_kzh_150
        need_kzh_200 = MPAS_field_will_be_written('kzh_200hPa')
        need_kzh = need_kzh .or. need_kzh_200
        need_any_diags = need_any_diags .or. need_kzh_200
        need_kzh_250 = MPAS_field_will_be_written('kzh_250hPa')
        need_kzh = need_kzh .or. need_kzh_250
        need_any_diags = need_any_diags .or. need_kzh_250
        need_kzh_300 = MPAS_field_will_be_written('kzh_300hPa')
        need_kzh = need_kzh .or. need_kzh_300
        need_any_diags = need_any_diags .or. need_kzh_300
        need_kzh_400 = MPAS_field_will_be_written('kzh_400hPa')
        need_kzh = need_kzh .or. need_kzh_400
        need_any_diags = need_any_diags .or. need_kzh_400
        need_kzh_500 = MPAS_field_will_be_written('kzh_500hPa')
        need_kzh = need_kzh .or. need_kzh_500
        need_any_diags = need_any_diags .or. need_kzh_500
        need_kzh_600 = MPAS_field_will_be_written('kzh_600hPa')
        need_kzh = need_kzh .or. need_kzh_600
        need_any_diags = need_any_diags .or. need_kzh_600
        need_kzh_700 = MPAS_field_will_be_written('kzh_700hPa')
        need_kzh = need_kzh .or. need_kzh_700
        need_any_diags = need_any_diags .or. need_kzh_700
        need_kzh_825 = MPAS_field_will_be_written('kzh_825hPa')
        need_kzh = need_kzh .or. need_kzh_825
        need_any_diags = need_any_diags .or. need_kzh_825
        need_kzh_850 = MPAS_field_will_be_written('kzh_850hPa')
        need_kzh = need_kzh .or. need_kzh_850
        need_any_diags = need_any_diags .or. need_kzh_850
        need_kzh_875 = MPAS_field_will_be_written('kzh_875hPa')
        need_kzh = need_kzh .or. need_kzh_875
        need_any_diags = need_any_diags .or. need_kzh_875
        need_kzh_900 = MPAS_field_will_be_written('kzh_900hPa')
        need_kzh = need_kzh .or. need_kzh_900
        need_any_diags = need_any_diags .or. need_kzh_900
        need_kzh_925 = MPAS_field_will_be_written('kzh_925hPa')
        need_kzh = need_kzh .or. need_kzh_925
        need_any_diags = need_any_diags .or. need_kzh_925
        need_kzh_950 = MPAS_field_will_be_written('kzh_950hPa')
        need_kzh = need_kzh .or. need_kzh_950
        need_any_diags = need_any_diags .or. need_kzh_950
        need_kzh_975 = MPAS_field_will_be_written('kzh_975hPa')
        need_kzh = need_kzh .or. need_kzh_975
        need_any_diags = need_any_diags .or. need_kzh_975
        need_kzh_1000 = MPAS_field_will_be_written('kzh_1000hPa')
        need_kzh = need_kzh .or. need_kzh_1000
        need_any_diags = need_any_diags .or. need_kzh_1000

        need_kzm_15 = MPAS_field_will_be_written('kzm_15hPa')
        need_kzm = need_kzm .or. need_kzm_15
        need_any_diags = need_any_diags .or. need_kzm_15
        need_kzm_20 = MPAS_field_will_be_written('kzm_20hPa')
        need_kzm = need_kzm .or. need_kzm_20
        need_any_diags = need_any_diags .or. need_kzm_20
        need_kzm_30 = MPAS_field_will_be_written('kzm_30hPa')
        need_kzm = need_kzm .or. need_kzm_30
        need_any_diags = need_any_diags .or. need_kzm_30
        need_kzm_50 = MPAS_field_will_be_written('kzm_50hPa')
        need_kzm = need_kzm .or. need_kzm_50
        need_any_diags = need_any_diags .or. need_kzm_50
        need_kzm_70 = MPAS_field_will_be_written('kzm_70hPa')
        need_kzm = need_kzm .or. need_kzm_70
        need_any_diags = need_any_diags .or. need_kzm_70
        need_kzm_100 = MPAS_field_will_be_written('kzm_100hPa')
        need_kzm = need_kzm .or. need_kzm_100
        need_any_diags = need_any_diags .or. need_kzm_100
        need_kzm_150 = MPAS_field_will_be_written('kzm_150hPa')
        need_kzm = need_kzm .or. need_kzm_150
        need_any_diags = need_any_diags .or. need_kzm_150
        need_kzm_200 = MPAS_field_will_be_written('kzm_200hPa')
        need_kzm = need_kzm .or. need_kzm_200
        need_any_diags = need_any_diags .or. need_kzm_200
        need_kzm_250 = MPAS_field_will_be_written('kzm_250hPa')
        need_kzm = need_kzm .or. need_kzm_250
        need_any_diags = need_any_diags .or. need_kzm_250
        need_kzm_300 = MPAS_field_will_be_written('kzm_300hPa')
        need_kzm = need_kzm .or. need_kzm_300
        need_any_diags = need_any_diags .or. need_kzm_300
        need_kzm_400 = MPAS_field_will_be_written('kzm_400hPa')
        need_kzm = need_kzm .or. need_kzm_400
        need_any_diags = need_any_diags .or. need_kzm_400
        need_kzm_500 = MPAS_field_will_be_written('kzm_500hPa')
        need_kzm = need_kzm .or. need_kzm_500
        need_any_diags = need_any_diags .or. need_kzm_500
        need_kzm_600 = MPAS_field_will_be_written('kzm_600hPa')
        need_kzm = need_kzm .or. need_kzm_600
        need_any_diags = need_any_diags .or. need_kzm_600
        need_kzm_700 = MPAS_field_will_be_written('kzm_700hPa')
        need_kzm = need_kzm .or. need_kzm_700
        need_any_diags = need_any_diags .or. need_kzm_700
        need_kzm_825 = MPAS_field_will_be_written('kzm_825hPa')
        need_kzm = need_kzm .or. need_kzm_825
        need_any_diags = need_any_diags .or. need_kzm_825
        need_kzm_850 = MPAS_field_will_be_written('kzm_850hPa')
        need_kzm = need_kzm .or. need_kzm_850
        need_any_diags = need_any_diags .or. need_kzm_850
        need_kzm_875 = MPAS_field_will_be_written('kzm_875hPa')
        need_kzm = need_kzm .or. need_kzm_875
        need_any_diags = need_any_diags .or. need_kzm_875
        need_kzm_900 = MPAS_field_will_be_written('kzm_900hPa')
        need_kzm = need_kzm .or. need_kzm_900
        need_any_diags = need_any_diags .or. need_kzm_900
        need_kzm_925 = MPAS_field_will_be_written('kzm_925hPa')
        need_kzm = need_kzm .or. need_kzm_925
        need_any_diags = need_any_diags .or. need_kzm_925
        need_kzm_950 = MPAS_field_will_be_written('kzm_950hPa')
        need_kzm = need_kzm .or. need_kzm_950
        need_any_diags = need_any_diags .or. need_kzm_950
        need_kzm_975 = MPAS_field_will_be_written('kzm_975hPa')
        need_kzm = need_kzm .or. need_kzm_975
        need_any_diags = need_any_diags .or. need_kzm_975
        need_kzm_1000 = MPAS_field_will_be_written('kzm_1000hPa')
        need_kzm = need_kzm .or. need_kzm_1000
        need_any_diags = need_any_diags .or. need_kzm_1000

        need_ni_15 = MPAS_field_will_be_written('ni_15hPa')
        need_ni = need_ni .or. need_ni_15
        need_any_diags = need_any_diags .or. need_ni_15
        need_ni_20 = MPAS_field_will_be_written('ni_20hPa')
        need_ni = need_ni .or. need_ni_20
        need_any_diags = need_any_diags .or. need_ni_20
        need_ni_30 = MPAS_field_will_be_written('ni_30hPa')
        need_ni = need_ni .or. need_ni_30
        need_any_diags = need_any_diags .or. need_ni_30
        need_ni_50 = MPAS_field_will_be_written('ni_50hPa')
        need_ni = need_ni .or. need_ni_50
        need_any_diags = need_any_diags .or. need_ni_50
        need_ni_70 = MPAS_field_will_be_written('ni_70hPa')
        need_ni = need_ni .or. need_ni_70
        need_any_diags = need_any_diags .or. need_ni_70
        need_ni_100 = MPAS_field_will_be_written('ni_100hPa')
        need_ni = need_ni .or. need_ni_100
        need_any_diags = need_any_diags .or. need_ni_100
        need_ni_150 = MPAS_field_will_be_written('ni_150hPa')
        need_ni = need_ni .or. need_ni_150
        need_any_diags = need_any_diags .or. need_ni_150
        need_ni_200 = MPAS_field_will_be_written('ni_200hPa')
        need_ni = need_ni .or. need_ni_200
        need_any_diags = need_any_diags .or. need_ni_200
        need_ni_250 = MPAS_field_will_be_written('ni_250hPa')
        need_ni = need_ni .or. need_ni_250
        need_any_diags = need_any_diags .or. need_ni_250
        need_ni_300 = MPAS_field_will_be_written('ni_300hPa')
        need_ni = need_ni .or. need_ni_300
        need_any_diags = need_any_diags .or. need_ni_300
        need_ni_400 = MPAS_field_will_be_written('ni_400hPa')
        need_ni = need_ni .or. need_ni_400
        need_any_diags = need_any_diags .or. need_ni_400
        need_ni_500 = MPAS_field_will_be_written('ni_500hPa')
        need_ni = need_ni .or. need_ni_500
        need_any_diags = need_any_diags .or. need_ni_500
        need_ni_600 = MPAS_field_will_be_written('ni_600hPa')
        need_ni = need_ni .or. need_ni_600
        need_any_diags = need_any_diags .or. need_ni_600
        need_ni_700 = MPAS_field_will_be_written('ni_700hPa')
        need_ni = need_ni .or. need_ni_700
        need_any_diags = need_any_diags .or. need_ni_700
        need_ni_825 = MPAS_field_will_be_written('ni_825hPa')
        need_ni = need_ni .or. need_ni_825
        need_any_diags = need_any_diags .or. need_ni_825
        need_ni_850 = MPAS_field_will_be_written('ni_850hPa')
        need_ni = need_ni .or. need_ni_850
        need_any_diags = need_any_diags .or. need_ni_850
        need_ni_875 = MPAS_field_will_be_written('ni_875hPa')
        need_ni = need_ni .or. need_ni_875
        need_any_diags = need_any_diags .or. need_ni_875
        need_ni_900 = MPAS_field_will_be_written('ni_900hPa')
        need_ni = need_ni .or. need_ni_900
        need_any_diags = need_any_diags .or. need_ni_900
        need_ni_925 = MPAS_field_will_be_written('ni_925hPa')
        need_ni = need_ni .or. need_ni_925
        need_any_diags = need_any_diags .or. need_ni_925
        need_ni_950 = MPAS_field_will_be_written('ni_950hPa')
        need_ni = need_ni .or. need_ni_950
        need_any_diags = need_any_diags .or. need_ni_950
        need_ni_975 = MPAS_field_will_be_written('ni_975hPa')
        need_ni = need_ni .or. need_ni_975
        need_any_diags = need_any_diags .or. need_ni_975
        need_ni_1000 = MPAS_field_will_be_written('ni_1000hPa')
        need_ni = need_ni .or. need_ni_1000
        need_any_diags = need_any_diags .or. need_ni_1000

        need_nr_15 = MPAS_field_will_be_written('nr_15hPa')
        need_nr = need_nr .or. need_nr_15
        need_any_diags = need_any_diags .or. need_nr_15
        need_nr_20 = MPAS_field_will_be_written('nr_20hPa')
        need_nr = need_nr .or. need_nr_20
        need_any_diags = need_any_diags .or. need_nr_20
        need_nr_30 = MPAS_field_will_be_written('nr_30hPa')
        need_nr = need_nr .or. need_nr_30
        need_any_diags = need_any_diags .or. need_nr_30
        need_nr_50 = MPAS_field_will_be_written('nr_50hPa')
        need_nr = need_nr .or. need_nr_50
        need_any_diags = need_any_diags .or. need_nr_50
        need_nr_70 = MPAS_field_will_be_written('nr_70hPa')
        need_nr = need_nr .or. need_nr_70
        need_any_diags = need_any_diags .or. need_nr_70
        need_nr_100 = MPAS_field_will_be_written('nr_100hPa')
        need_nr = need_nr .or. need_nr_100
        need_any_diags = need_any_diags .or. need_nr_100
        need_nr_150 = MPAS_field_will_be_written('nr_150hPa')
        need_nr = need_nr .or. need_nr_150
        need_any_diags = need_any_diags .or. need_nr_150
        need_nr_200 = MPAS_field_will_be_written('nr_200hPa')
        need_nr = need_nr .or. need_nr_200
        need_any_diags = need_any_diags .or. need_nr_200
        need_nr_250 = MPAS_field_will_be_written('nr_250hPa')
        need_nr = need_nr .or. need_nr_250
        need_any_diags = need_any_diags .or. need_nr_250
        need_nr_300 = MPAS_field_will_be_written('nr_300hPa')
        need_nr = need_nr .or. need_nr_300
        need_any_diags = need_any_diags .or. need_nr_300
        need_nr_400 = MPAS_field_will_be_written('nr_400hPa')
        need_nr = need_nr .or. need_nr_400
        need_any_diags = need_any_diags .or. need_nr_400
        need_nr_500 = MPAS_field_will_be_written('nr_500hPa')
        need_nr = need_nr .or. need_nr_500
        need_any_diags = need_any_diags .or. need_nr_500
        need_nr_600 = MPAS_field_will_be_written('nr_600hPa')
        need_nr = need_nr .or. need_nr_600
        need_any_diags = need_any_diags .or. need_nr_600
        need_nr_700 = MPAS_field_will_be_written('nr_700hPa')
        need_nr = need_nr .or. need_nr_700
        need_any_diags = need_any_diags .or. need_nr_700
        need_nr_825 = MPAS_field_will_be_written('nr_825hPa')
        need_nr = need_nr .or. need_nr_825
        need_any_diags = need_any_diags .or. need_nr_825
        need_nr_850 = MPAS_field_will_be_written('nr_850hPa')
        need_nr = need_nr .or. need_nr_850
        need_any_diags = need_any_diags .or. need_nr_850
        need_nr_875 = MPAS_field_will_be_written('nr_875hPa')
        need_nr = need_nr .or. need_nr_875
        need_any_diags = need_any_diags .or. need_nr_875
        need_nr_900 = MPAS_field_will_be_written('nr_900hPa')
        need_nr = need_nr .or. need_nr_900
        need_any_diags = need_any_diags .or. need_nr_900
        need_nr_925 = MPAS_field_will_be_written('nr_925hPa')
        need_nr = need_nr .or. need_nr_925
        need_any_diags = need_any_diags .or. need_nr_925
        need_nr_950 = MPAS_field_will_be_written('nr_950hPa')
        need_nr = need_nr .or. need_nr_950
        need_any_diags = need_any_diags .or. need_nr_950
        need_nr_975 = MPAS_field_will_be_written('nr_975hPa')
        need_nr = need_nr .or. need_nr_975
        need_any_diags = need_any_diags .or. need_nr_975
        need_nr_1000 = MPAS_field_will_be_written('nr_1000hPa')
        need_nr = need_nr .or. need_nr_1000
        need_any_diags = need_any_diags .or. need_nr_1000

        need_qr_15 = MPAS_field_will_be_written('qr_15hPa')
        need_qr = need_qr .or. need_qr_15
        need_any_diags = need_any_diags .or. need_qr_15
        need_qr_20 = MPAS_field_will_be_written('qr_20hPa')
        need_qr = need_qr .or. need_qr_20
        need_any_diags = need_any_diags .or. need_qr_20
        need_qr_30 = MPAS_field_will_be_written('qr_30hPa')
        need_qr = need_qr .or. need_qr_30
        need_any_diags = need_any_diags .or. need_qr_30
        need_qr_50 = MPAS_field_will_be_written('qr_50hPa')
        need_qr = need_qr .or. need_qr_50
        need_any_diags = need_any_diags .or. need_qr_50
        need_qr_70 = MPAS_field_will_be_written('qr_70hPa')
        need_qr = need_qr .or. need_qr_70
        need_any_diags = need_any_diags .or. need_qr_70
        need_qr_100 = MPAS_field_will_be_written('qr_100hPa')
        need_qr = need_qr .or. need_qr_100
        need_any_diags = need_any_diags .or. need_qr_100
        need_qr_150 = MPAS_field_will_be_written('qr_150hPa')
        need_qr = need_qr .or. need_qr_150
        need_any_diags = need_any_diags .or. need_qr_150
        need_qr_200 = MPAS_field_will_be_written('qr_200hPa')
        need_qr = need_qr .or. need_qr_200
        need_any_diags = need_any_diags .or. need_qr_200
        need_qr_250 = MPAS_field_will_be_written('qr_250hPa')
        need_qr = need_qr .or. need_qr_250
        need_any_diags = need_any_diags .or. need_qr_250
        need_qr_300 = MPAS_field_will_be_written('qr_300hPa')
        need_qr = need_qr .or. need_qr_300
        need_any_diags = need_any_diags .or. need_qr_300
        need_qr_400 = MPAS_field_will_be_written('qr_400hPa')
        need_qr = need_qr .or. need_qr_400
        need_any_diags = need_any_diags .or. need_qr_400
        need_qr_500 = MPAS_field_will_be_written('qr_500hPa')
        need_qr = need_qr .or. need_qr_500
        need_any_diags = need_any_diags .or. need_qr_500
        need_qr_600 = MPAS_field_will_be_written('qr_600hPa')
        need_qr = need_qr .or. need_qr_600
        need_any_diags = need_any_diags .or. need_qr_600
        need_qr_700 = MPAS_field_will_be_written('qr_700hPa')
        need_qr = need_qr .or. need_qr_700
        need_any_diags = need_any_diags .or. need_qr_700
        need_qr_825 = MPAS_field_will_be_written('qr_825hPa')
        need_qr = need_qr .or. need_qr_825
        need_any_diags = need_any_diags .or. need_qr_825
        need_qr_850 = MPAS_field_will_be_written('qr_850hPa')
        need_qr = need_qr .or. need_qr_850
        need_any_diags = need_any_diags .or. need_qr_850
        need_qr_875 = MPAS_field_will_be_written('qr_875hPa')
        need_qr = need_qr .or. need_qr_875
        need_any_diags = need_any_diags .or. need_qr_875
        need_qr_900 = MPAS_field_will_be_written('qr_900hPa')
        need_qr = need_qr .or. need_qr_900
        need_any_diags = need_any_diags .or. need_qr_900
        need_qr_925 = MPAS_field_will_be_written('qr_925hPa')
        need_qr = need_qr .or. need_qr_925
        need_any_diags = need_any_diags .or. need_qr_925
        need_qr_950 = MPAS_field_will_be_written('qr_950hPa')
        need_qr = need_qr .or. need_qr_950
        need_any_diags = need_any_diags .or. need_qr_950
        need_qr_975 = MPAS_field_will_be_written('qr_975hPa')
        need_qr = need_qr .or. need_qr_975
        need_any_diags = need_any_diags .or. need_qr_975
        need_qr_1000 = MPAS_field_will_be_written('qr_1000hPa')
        need_qr = need_qr .or. need_qr_1000
        need_any_diags = need_any_diags .or. need_qr_1000

        need_qs_15 = MPAS_field_will_be_written('qs_15hPa')
        need_qs = need_qs .or. need_qs_15
        need_any_diags = need_any_diags .or. need_qs_15
        need_qs_20 = MPAS_field_will_be_written('qs_20hPa')
        need_qs = need_qs .or. need_qs_20
        need_any_diags = need_any_diags .or. need_qs_20
        need_qs_30 = MPAS_field_will_be_written('qs_30hPa')
        need_qs = need_qs .or. need_qs_30
        need_any_diags = need_any_diags .or. need_qs_30
        need_qs_50 = MPAS_field_will_be_written('qs_50hPa')
        need_qs = need_qs .or. need_qs_50
        need_any_diags = need_any_diags .or. need_qs_50
        need_qs_70 = MPAS_field_will_be_written('qs_70hPa')
        need_qs = need_qs .or. need_qs_70
        need_any_diags = need_any_diags .or. need_qs_70
        need_qs_100 = MPAS_field_will_be_written('qs_100hPa')
        need_qs = need_qs .or. need_qs_100
        need_any_diags = need_any_diags .or. need_qs_100
        need_qs_150 = MPAS_field_will_be_written('qs_150hPa')
        need_qs = need_qs .or. need_qs_150
        need_any_diags = need_any_diags .or. need_qs_150
        need_qs_200 = MPAS_field_will_be_written('qs_200hPa')
        need_qs = need_qs .or. need_qs_200
        need_any_diags = need_any_diags .or. need_qs_200
        need_qs_250 = MPAS_field_will_be_written('qs_250hPa')
        need_qs = need_qs .or. need_qs_250
        need_any_diags = need_any_diags .or. need_qs_250
        need_qs_300 = MPAS_field_will_be_written('qs_300hPa')
        need_qs = need_qs .or. need_qs_300
        need_any_diags = need_any_diags .or. need_qs_300
        need_qs_400 = MPAS_field_will_be_written('qs_400hPa')
        need_qs = need_qs .or. need_qs_400
        need_any_diags = need_any_diags .or. need_qs_400
        need_qs_500 = MPAS_field_will_be_written('qs_500hPa')
        need_qs = need_qs .or. need_qs_500
        need_any_diags = need_any_diags .or. need_qs_500
        need_qs_600 = MPAS_field_will_be_written('qs_600hPa')
        need_qs = need_qs .or. need_qs_600
        need_any_diags = need_any_diags .or. need_qs_600
        need_qs_700 = MPAS_field_will_be_written('qs_700hPa')
        need_qs = need_qs .or. need_qs_700
        need_any_diags = need_any_diags .or. need_qs_700
        need_qs_825 = MPAS_field_will_be_written('qs_825hPa')
        need_qs = need_qs .or. need_qs_825
        need_any_diags = need_any_diags .or. need_qs_825
        need_qs_850 = MPAS_field_will_be_written('qs_850hPa')
        need_qs = need_qs .or. need_qs_850
        need_any_diags = need_any_diags .or. need_qs_850
        need_qs_875 = MPAS_field_will_be_written('qs_875hPa')
        need_qs = need_qs .or. need_qs_875
        need_any_diags = need_any_diags .or. need_qs_875
        need_qs_900 = MPAS_field_will_be_written('qs_900hPa')
        need_qs = need_qs .or. need_qs_900
        need_any_diags = need_any_diags .or. need_qs_900
        need_qs_925 = MPAS_field_will_be_written('qs_925hPa')
        need_qs = need_qs .or. need_qs_925
        need_any_diags = need_any_diags .or. need_qs_925
        need_qs_950 = MPAS_field_will_be_written('qs_950hPa')
        need_qs = need_qs .or. need_qs_950
        need_any_diags = need_any_diags .or. need_qs_950
        need_qs_975 = MPAS_field_will_be_written('qs_975hPa')
        need_qs = need_qs .or. need_qs_975
        need_any_diags = need_any_diags .or. need_qs_975
        need_qs_1000 = MPAS_field_will_be_written('qs_1000hPa')
        need_qs = need_qs .or. need_qs_1000
        need_any_diags = need_any_diags .or. need_qs_1000

        need_qv_15 = MPAS_field_will_be_written('qv_15hPa')
        need_qv = need_qv .or. need_qv_15
        need_any_diags = need_any_diags .or. need_qv_15
        need_qv_20 = MPAS_field_will_be_written('qv_20hPa')
        need_qv = need_qv .or. need_qv_20
        need_any_diags = need_any_diags .or. need_qv_20
        need_qv_30 = MPAS_field_will_be_written('qv_30hPa')
        need_qv = need_qv .or. need_qv_30
        need_any_diags = need_any_diags .or. need_qv_30
        need_qv_50 = MPAS_field_will_be_written('qv_50hPa')
        need_qv = need_qv .or. need_qv_50
        need_any_diags = need_any_diags .or. need_qv_50
        need_qv_70 = MPAS_field_will_be_written('qv_70hPa')
        need_qv = need_qv .or. need_qv_70
        need_any_diags = need_any_diags .or. need_qv_70
        need_qv_100 = MPAS_field_will_be_written('qv_100hPa')
        need_qv = need_qv .or. need_qv_100
        need_any_diags = need_any_diags .or. need_qv_100
        need_qv_150 = MPAS_field_will_be_written('qv_150hPa')
        need_qv = need_qv .or. need_qv_150
        need_any_diags = need_any_diags .or. need_qv_150
        need_qv_200 = MPAS_field_will_be_written('qv_200hPa')
        need_qv = need_qv .or. need_qv_200
        need_any_diags = need_any_diags .or. need_qv_200
        need_qv_250 = MPAS_field_will_be_written('qv_250hPa')
        need_qv = need_qv .or. need_qv_250
        need_any_diags = need_any_diags .or. need_qv_250
        need_qv_300 = MPAS_field_will_be_written('qv_300hPa')
        need_qv = need_qv .or. need_qv_300
        need_any_diags = need_any_diags .or. need_qv_300
        need_qv_400 = MPAS_field_will_be_written('qv_400hPa')
        need_qv = need_qv .or. need_qv_400
        need_any_diags = need_any_diags .or. need_qv_400
        need_qv_500 = MPAS_field_will_be_written('qv_500hPa')
        need_qv = need_qv .or. need_qv_500
        need_any_diags = need_any_diags .or. need_qv_500
        need_qv_600 = MPAS_field_will_be_written('qv_600hPa')
        need_qv = need_qv .or. need_qv_600
        need_any_diags = need_any_diags .or. need_qv_600
        need_qv_700 = MPAS_field_will_be_written('qv_700hPa')
        need_qv = need_qv .or. need_qv_700
        need_any_diags = need_any_diags .or. need_qv_700
        need_qv_825 = MPAS_field_will_be_written('qv_825hPa')
        need_qv = need_qv .or. need_qv_825
        need_any_diags = need_any_diags .or. need_qv_825
        need_qv_850 = MPAS_field_will_be_written('qv_850hPa')
        need_qv = need_qv .or. need_qv_850
        need_any_diags = need_any_diags .or. need_qv_850
        need_qv_875 = MPAS_field_will_be_written('qv_875hPa')
        need_qv = need_qv .or. need_qv_875
        need_any_diags = need_any_diags .or. need_qv_875
        need_qv_900 = MPAS_field_will_be_written('qv_900hPa')
        need_qv = need_qv .or. need_qv_900
        need_any_diags = need_any_diags .or. need_qv_900
        need_qv_925 = MPAS_field_will_be_written('qv_925hPa')
        need_qv = need_qv .or. need_qv_925
        need_any_diags = need_any_diags .or. need_qv_925
        need_qv_950 = MPAS_field_will_be_written('qv_950hPa')
        need_qv = need_qv .or. need_qv_950
        need_any_diags = need_any_diags .or. need_qv_950
        need_qv_975 = MPAS_field_will_be_written('qv_975hPa')
        need_qv = need_qv .or. need_qv_975
        need_any_diags = need_any_diags .or. need_qv_975
        need_qv_1000 = MPAS_field_will_be_written('qv_1000hPa')
        need_qv = need_qv .or. need_qv_1000
        need_any_diags = need_any_diags .or. need_qv_1000

        need_qc_15 = MPAS_field_will_be_written('qc_15hPa')
        need_qc = need_qc .or. need_qc_15
        need_any_diags = need_any_diags .or. need_qc_15
        need_qc_20 = MPAS_field_will_be_written('qc_20hPa')
        need_qc = need_qc .or. need_qc_20
        need_any_diags = need_any_diags .or. need_qc_20
        need_qc_30 = MPAS_field_will_be_written('qc_30hPa')
        need_qc = need_qc .or. need_qc_30
        need_any_diags = need_any_diags .or. need_qc_30
        need_qc_50 = MPAS_field_will_be_written('qc_50hPa')
        need_qc = need_qc .or. need_qc_50
        need_any_diags = need_any_diags .or. need_qc_50
        need_qc_70 = MPAS_field_will_be_written('qc_70hPa')
        need_qc = need_qc .or. need_qc_70
        need_any_diags = need_any_diags .or. need_qc_70
        need_qc_100 = MPAS_field_will_be_written('qc_100hPa')
        need_qc = need_qc .or. need_qc_100
        need_any_diags = need_any_diags .or. need_qc_100
        need_qc_150 = MPAS_field_will_be_written('qc_150hPa')
        need_qc = need_qc .or. need_qc_150
        need_any_diags = need_any_diags .or. need_qc_150
        need_qc_200 = MPAS_field_will_be_written('qc_200hPa')
        need_qc = need_qc .or. need_qc_200
        need_any_diags = need_any_diags .or. need_qc_200
        need_qc_250 = MPAS_field_will_be_written('qc_250hPa')
        need_qc = need_qc .or. need_qc_250
        need_any_diags = need_any_diags .or. need_qc_250
        need_qc_300 = MPAS_field_will_be_written('qc_300hPa')
        need_qc = need_qc .or. need_qc_300
        need_any_diags = need_any_diags .or. need_qc_300
        need_qc_400 = MPAS_field_will_be_written('qc_400hPa')
        need_qc = need_qc .or. need_qc_400
        need_any_diags = need_any_diags .or. need_qc_400
        need_qc_500 = MPAS_field_will_be_written('qc_500hPa')
        need_qc = need_qc .or. need_qc_500
        need_any_diags = need_any_diags .or. need_qc_500
        need_qc_600 = MPAS_field_will_be_written('qc_600hPa')
        need_qc = need_qc .or. need_qc_600
        need_any_diags = need_any_diags .or. need_qc_600
        need_qc_700 = MPAS_field_will_be_written('qc_700hPa')
        need_qc = need_qc .or. need_qc_700
        need_any_diags = need_any_diags .or. need_qc_700
        need_qc_825 = MPAS_field_will_be_written('qc_825hPa')
        need_qc = need_qc .or. need_qc_825
        need_any_diags = need_any_diags .or. need_qc_825
        need_qc_850 = MPAS_field_will_be_written('qc_850hPa')
        need_qc = need_qc .or. need_qc_850
        need_any_diags = need_any_diags .or. need_qc_850
        need_qc_875 = MPAS_field_will_be_written('qc_875hPa')
        need_qc = need_qc .or. need_qc_875
        need_any_diags = need_any_diags .or. need_qc_875
        need_qc_900 = MPAS_field_will_be_written('qc_900hPa')
        need_qc = need_qc .or. need_qc_900
        need_any_diags = need_any_diags .or. need_qc_900
        need_qc_925 = MPAS_field_will_be_written('qc_925hPa')
        need_qc = need_qc .or. need_qc_925
        need_any_diags = need_any_diags .or. need_qc_925
        need_qc_950 = MPAS_field_will_be_written('qc_950hPa')
        need_qc = need_qc .or. need_qc_950
        need_any_diags = need_any_diags .or. need_qc_950
        need_qc_975 = MPAS_field_will_be_written('qc_975hPa')
        need_qc = need_qc .or. need_qc_975
        need_any_diags = need_any_diags .or. need_qc_975
        need_qc_1000 = MPAS_field_will_be_written('qc_1000hPa')
        need_qc = need_qc .or. need_qc_1000
        need_any_diags = need_any_diags .or. need_qc_1000

        need_qg_15 = MPAS_field_will_be_written('qg_15hPa')
        need_qg = need_qg .or. need_qg_15
        need_any_diags = need_any_diags .or. need_qg_15
        need_qg_20 = MPAS_field_will_be_written('qg_20hPa')
        need_qg = need_qg .or. need_qg_20
        need_any_diags = need_any_diags .or. need_qg_20
        need_qg_30 = MPAS_field_will_be_written('qg_30hPa')
        need_qg = need_qg .or. need_qg_30
        need_any_diags = need_any_diags .or. need_qg_30
        need_qg_50 = MPAS_field_will_be_written('qg_50hPa')
        need_qg = need_qg .or. need_qg_50
        need_any_diags = need_any_diags .or. need_qg_50
        need_qg_70 = MPAS_field_will_be_written('qg_70hPa')
        need_qg = need_qg .or. need_qg_70
        need_any_diags = need_any_diags .or. need_qg_70
        need_qg_100 = MPAS_field_will_be_written('qg_100hPa')
        need_qg = need_qg .or. need_qg_100
        need_any_diags = need_any_diags .or. need_qg_100
        need_qg_150 = MPAS_field_will_be_written('qg_150hPa')
        need_qg = need_qg .or. need_qg_150
        need_any_diags = need_any_diags .or. need_qg_150
        need_qg_200 = MPAS_field_will_be_written('qg_200hPa')
        need_qg = need_qg .or. need_qg_200
        need_any_diags = need_any_diags .or. need_qg_200
        need_qg_250 = MPAS_field_will_be_written('qg_250hPa')
        need_qg = need_qg .or. need_qg_250
        need_any_diags = need_any_diags .or. need_qg_250
        need_qg_300 = MPAS_field_will_be_written('qg_300hPa')
        need_qg = need_qg .or. need_qg_300
        need_any_diags = need_any_diags .or. need_qg_300
        need_qg_400 = MPAS_field_will_be_written('qg_400hPa')
        need_qg = need_qg .or. need_qg_400
        need_any_diags = need_any_diags .or. need_qg_400
        need_qg_500 = MPAS_field_will_be_written('qg_500hPa')
        need_qg = need_qg .or. need_qg_500
        need_any_diags = need_any_diags .or. need_qg_500
        need_qg_600 = MPAS_field_will_be_written('qg_600hPa')
        need_qg = need_qg .or. need_qg_600
        need_any_diags = need_any_diags .or. need_qg_600
        need_qg_700 = MPAS_field_will_be_written('qg_700hPa')
        need_qg = need_qg .or. need_qg_700
        need_any_diags = need_any_diags .or. need_qg_700
        need_qg_825 = MPAS_field_will_be_written('qg_825hPa')
        need_qg = need_qg .or. need_qg_825
        need_any_diags = need_any_diags .or. need_qg_825
        need_qg_850 = MPAS_field_will_be_written('qg_850hPa')
        need_qg = need_qg .or. need_qg_850
        need_any_diags = need_any_diags .or. need_qg_850
        need_qg_875 = MPAS_field_will_be_written('qg_875hPa')
        need_qg = need_qg .or. need_qg_875
        need_any_diags = need_any_diags .or. need_qg_875
        need_qg_900 = MPAS_field_will_be_written('qg_900hPa')
        need_qg = need_qg .or. need_qg_900
        need_any_diags = need_any_diags .or. need_qg_900
        need_qg_925 = MPAS_field_will_be_written('qg_925hPa')
        need_qg = need_qg .or. need_qg_925
        need_any_diags = need_any_diags .or. need_qg_925
        need_qg_950 = MPAS_field_will_be_written('qg_950hPa')
        need_qg = need_qg .or. need_qg_950
        need_any_diags = need_any_diags .or. need_qg_950
        need_qg_975 = MPAS_field_will_be_written('qg_975hPa')
        need_qg = need_qg .or. need_qg_975
        need_any_diags = need_any_diags .or. need_qg_975
        need_qg_1000 = MPAS_field_will_be_written('qg_1000hPa')
        need_qg = need_qg .or. need_qg_1000
        need_any_diags = need_any_diags .or. need_qg_1000

        need_qi_15 = MPAS_field_will_be_written('qi_15hPa')
        need_qi = need_qi .or. need_qi_15
        need_any_diags = need_any_diags .or. need_qi_15
        need_qi_20 = MPAS_field_will_be_written('qi_20hPa')
        need_qi = need_qi .or. need_qi_20
        need_any_diags = need_any_diags .or. need_qi_20
        need_qi_30 = MPAS_field_will_be_written('qi_30hPa')
        need_qi = need_qi .or. need_qi_30
        need_any_diags = need_any_diags .or. need_qi_30
        need_qi_50 = MPAS_field_will_be_written('qi_50hPa')
        need_qi = need_qi .or. need_qi_50
        need_any_diags = need_any_diags .or. need_qi_50
        need_qi_70 = MPAS_field_will_be_written('qi_70hPa')
        need_qi = need_qi .or. need_qi_70
        need_any_diags = need_any_diags .or. need_qi_70
        need_qi_100 = MPAS_field_will_be_written('qi_100hPa')
        need_qi = need_qi .or. need_qi_100
        need_any_diags = need_any_diags .or. need_qi_100
        need_qi_150 = MPAS_field_will_be_written('qi_150hPa')
        need_qi = need_qi .or. need_qi_150
        need_any_diags = need_any_diags .or. need_qi_150
        need_qi_200 = MPAS_field_will_be_written('qi_200hPa')
        need_qi = need_qi .or. need_qi_200
        need_any_diags = need_any_diags .or. need_qi_200
        need_qi_250 = MPAS_field_will_be_written('qi_250hPa')
        need_qi = need_qi .or. need_qi_250
        need_any_diags = need_any_diags .or. need_qi_250
        need_qi_300 = MPAS_field_will_be_written('qi_300hPa')
        need_qi = need_qi .or. need_qi_300
        need_any_diags = need_any_diags .or. need_qi_300
        need_qi_400 = MPAS_field_will_be_written('qi_400hPa')
        need_qi = need_qi .or. need_qi_400
        need_any_diags = need_any_diags .or. need_qi_400
        need_qi_500 = MPAS_field_will_be_written('qi_500hPa')
        need_qi = need_qi .or. need_qi_500
        need_any_diags = need_any_diags .or. need_qi_500
        need_qi_600 = MPAS_field_will_be_written('qi_600hPa')
        need_qi = need_qi .or. need_qi_600
        need_any_diags = need_any_diags .or. need_qi_600
        need_qi_700 = MPAS_field_will_be_written('qi_700hPa')
        need_qi = need_qi .or. need_qi_700
        need_any_diags = need_any_diags .or. need_qi_700
        need_qi_825 = MPAS_field_will_be_written('qi_825hPa')
        need_qi = need_qi .or. need_qi_825
        need_any_diags = need_any_diags .or. need_qi_825
        need_qi_850 = MPAS_field_will_be_written('qi_850hPa')
        need_qi = need_qi .or. need_qi_850
        need_any_diags = need_any_diags .or. need_qi_850
        need_qi_875 = MPAS_field_will_be_written('qi_875hPa')
        need_qi = need_qi .or. need_qi_875
        need_any_diags = need_any_diags .or. need_qi_875
        need_qi_900 = MPAS_field_will_be_written('qi_900hPa')
        need_qi = need_qi .or. need_qi_900
        need_any_diags = need_any_diags .or. need_qi_900
        need_qi_925 = MPAS_field_will_be_written('qi_925hPa')
        need_qi = need_qi .or. need_qi_925
        need_any_diags = need_any_diags .or. need_qi_925
        need_qi_950 = MPAS_field_will_be_written('qi_950hPa')
        need_qi = need_qi .or. need_qi_950
        need_any_diags = need_any_diags .or. need_qi_950
        need_qi_975 = MPAS_field_will_be_written('qi_975hPa')
        need_qi = need_qi .or. need_qi_975
        need_any_diags = need_any_diags .or. need_qi_975
        need_qi_1000 = MPAS_field_will_be_written('qi_1000hPa')
        need_qi = need_qi .or. need_qi_1000
        need_any_diags = need_any_diags .or. need_qi_1000

        need_vorticity_50 = MPAS_field_will_be_written('vorticity_50hPa')
        need_vorticity = need_vorticity .or. need_vorticity_50
        need_any_diags = need_any_diags .or. need_vorticity_50
        need_vorticity_100 = MPAS_field_will_be_written('vorticity_100hPa')
        need_vorticity = need_vorticity .or. need_vorticity_100
        need_any_diags = need_any_diags .or. need_vorticity_100
        need_vorticity_200 = MPAS_field_will_be_written('vorticity_200hPa')
        need_vorticity = need_vorticity .or. need_vorticity_200
        need_any_diags = need_any_diags .or. need_vorticity_200
        need_vorticity_250 = MPAS_field_will_be_written('vorticity_250hPa')
        need_vorticity = need_vorticity .or. need_vorticity_250
        need_any_diags = need_any_diags .or. need_vorticity_250
        need_vorticity_500 = MPAS_field_will_be_written('vorticity_500hPa')
        need_vorticity = need_vorticity .or. need_vorticity_500
        need_any_diags = need_any_diags .or. need_vorticity_500
        need_vorticity_700 = MPAS_field_will_be_written('vorticity_700hPa')
        need_vorticity = need_vorticity .or. need_vorticity_700
        need_any_diags = need_any_diags .or. need_vorticity_700
        need_vorticity_850 = MPAS_field_will_be_written('vorticity_850hPa')
        need_vorticity = need_vorticity .or. need_vorticity_850
        need_any_diags = need_any_diags .or. need_vorticity_850
        need_vorticity_925 = MPAS_field_will_be_written('vorticity_925hPa')
        need_vorticity = need_vorticity .or. need_vorticity_925
        need_any_diags = need_any_diags .or. need_vorticity_925
        need_t_isobaric = MPAS_field_will_be_written('t_isobaric')
        need_any_diags = need_any_diags .or. need_t_isobaric
        need_z_isobaric = MPAS_field_will_be_written('z_isobaric')
        need_any_diags = need_any_diags .or. need_z_isobaric
        need_meanT_500_300 = MPAS_field_will_be_written('meanT_500_300')
        need_any_diags = need_any_diags .or. need_meanT_500_300

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
        integer :: iCell,iVert,iVertD,k,kk
        integer, pointer :: nCells, nCellsSolve, nVertLevels, nVertices, vertexDegree, nIsoLevelsT, nIsoLevelsZ
        integer :: nVertLevelsP1
        integer, pointer :: index_qv, index_ni, index_nr, num_scalars, index_qr, index_qs, index_qc, index_qg, index_qi
        integer, dimension(:,:), pointer :: cellsOnVertex

        type (field2DReal), pointer:: pressure_p_field
       
        real (kind=RKIND), dimension(:), pointer :: areaTriangle, latCell
        real (kind=RKIND), dimension(:,:), pointer :: kiteAreasOnVertex
        
        real (kind=RKIND), dimension(:,:), pointer :: exner, height
        real (kind=RKIND), dimension(:,:), pointer :: pressure_b, pressure_p 
        real (kind=RKIND), dimension(:,:), pointer :: relhum, theta_m, vorticity
        real (kind=RKIND), dimension(:,:), pointer :: umeridional, uzonal, vvel
        real (kind=RKIND), dimension(:,:), pointer :: rqvblten, rqvcuten, rthratenlw
        real (kind=RKIND), dimension(:,:), pointer :: rucuten
        real (kind=RKIND), dimension(:,:), pointer :: cldfrac
        real (kind=RKIND), dimension(:,:), pointer :: kzq, kzh, kzm
        real (kind=RKIND), dimension(:,:,:), pointer :: scalars
       
        real (kind=RKIND), dimension(:), pointer :: t_iso_levels
        real (kind=RKIND), dimension(:), pointer :: z_iso_levels
        real (kind=RKIND), dimension(:,:), pointer :: t_isobaric
        real (kind=RKIND), dimension(:,:), pointer :: z_isobaric
        real (kind=RKIND), dimension(:), pointer :: meanT_500_300

        real (kind=RKIND), dimension(:), pointer :: cldfrac_15hPa
        real (kind=RKIND), dimension(:), pointer :: cldfrac_20hPa
        real (kind=RKIND), dimension(:), pointer :: cldfrac_30hPa  
        real (kind=RKIND), dimension(:), pointer :: cldfrac_50hPa
        real (kind=RKIND), dimension(:), pointer :: cldfrac_70hPa
        real (kind=RKIND), dimension(:), pointer :: cldfrac_100hPa
        real (kind=RKIND), dimension(:), pointer :: cldfrac_150hPa
        real (kind=RKIND), dimension(:), pointer :: cldfrac_200hPa
        real (kind=RKIND), dimension(:), pointer :: cldfrac_250hPa
        real (kind=RKIND), dimension(:), pointer :: cldfrac_300hPa
        real (kind=RKIND), dimension(:), pointer :: cldfrac_400hPa
        real (kind=RKIND), dimension(:), pointer :: cldfrac_500hPa
        real (kind=RKIND), dimension(:), pointer :: cldfrac_600hPa
        real (kind=RKIND), dimension(:), pointer :: cldfrac_700hPa
        real (kind=RKIND), dimension(:), pointer :: cldfrac_825hPa
        real (kind=RKIND), dimension(:), pointer :: cldfrac_850hPa
        real (kind=RKIND), dimension(:), pointer :: cldfrac_875hPa
        real (kind=RKIND), dimension(:), pointer :: cldfrac_900hPa
        real (kind=RKIND), dimension(:), pointer :: cldfrac_925hPa
        real (kind=RKIND), dimension(:), pointer :: cldfrac_950hPa
        real (kind=RKIND), dimension(:), pointer :: cldfrac_975hPa
        real (kind=RKIND), dimension(:), pointer :: cldfrac_1000hPa

        real (kind=RKIND), dimension(:), pointer :: temperature_15hPa
        real (kind=RKIND), dimension(:), pointer :: temperature_20hPa
        real (kind=RKIND), dimension(:), pointer :: temperature_30hPa  
        real (kind=RKIND), dimension(:), pointer :: temperature_50hPa
        real (kind=RKIND), dimension(:), pointer :: temperature_70hPa
        real (kind=RKIND), dimension(:), pointer :: temperature_100hPa
        real (kind=RKIND), dimension(:), pointer :: temperature_150hPa
        real (kind=RKIND), dimension(:), pointer :: temperature_200hPa
        real (kind=RKIND), dimension(:), pointer :: temperature_250hPa
        real (kind=RKIND), dimension(:), pointer :: temperature_300hPa
        real (kind=RKIND), dimension(:), pointer :: temperature_400hPa
        real (kind=RKIND), dimension(:), pointer :: temperature_500hPa
        real (kind=RKIND), dimension(:), pointer :: temperature_600hPa
        real (kind=RKIND), dimension(:), pointer :: temperature_700hPa
        real (kind=RKIND), dimension(:), pointer :: temperature_825hPa
        real (kind=RKIND), dimension(:), pointer :: temperature_850hPa
        real (kind=RKIND), dimension(:), pointer :: temperature_875hPa
        real (kind=RKIND), dimension(:), pointer :: temperature_900hPa
        real (kind=RKIND), dimension(:), pointer :: temperature_925hPa
        real (kind=RKIND), dimension(:), pointer :: temperature_950hPa
        real (kind=RKIND), dimension(:), pointer :: temperature_975hPa
        real (kind=RKIND), dimension(:), pointer :: temperature_1000hPa

        real (kind=RKIND), dimension(:), pointer :: relhum_15hPa
        real (kind=RKIND), dimension(:), pointer :: relhum_20hPa
        real (kind=RKIND), dimension(:), pointer :: relhum_30hPa
        real (kind=RKIND), dimension(:), pointer :: relhum_50hPa
        real (kind=RKIND), dimension(:), pointer :: relhum_70hPa
        real (kind=RKIND), dimension(:), pointer :: relhum_100hPa
        real (kind=RKIND), dimension(:), pointer :: relhum_150hPa
        real (kind=RKIND), dimension(:), pointer :: relhum_200hPa
        real (kind=RKIND), dimension(:), pointer :: relhum_250hPa
        real (kind=RKIND), dimension(:), pointer :: relhum_300hPa
        real (kind=RKIND), dimension(:), pointer :: relhum_400hPa
        real (kind=RKIND), dimension(:), pointer :: relhum_500hPa
        real (kind=RKIND), dimension(:), pointer :: relhum_600hPa
        real (kind=RKIND), dimension(:), pointer :: relhum_700hPa
        real (kind=RKIND), dimension(:), pointer :: relhum_825hPa
        real (kind=RKIND), dimension(:), pointer :: relhum_850hPa
        real (kind=RKIND), dimension(:), pointer :: relhum_875hPa
        real (kind=RKIND), dimension(:), pointer :: relhum_900hPa
        real (kind=RKIND), dimension(:), pointer :: relhum_925hPa
        real (kind=RKIND), dimension(:), pointer :: relhum_950hPa
        real (kind=RKIND), dimension(:), pointer :: relhum_975hPa
        real (kind=RKIND), dimension(:), pointer :: relhum_1000hPa        
       
        real (kind=RKIND), dimension(:), pointer :: dewpoint_50hPa
        real (kind=RKIND), dimension(:), pointer :: dewpoint_100hPa
        real (kind=RKIND), dimension(:), pointer :: dewpoint_200hPa
        real (kind=RKIND), dimension(:), pointer :: dewpoint_250hPa
        real (kind=RKIND), dimension(:), pointer :: dewpoint_500hPa
        real (kind=RKIND), dimension(:), pointer :: dewpoint_700hPa
        real (kind=RKIND), dimension(:), pointer :: dewpoint_850hPa
        real (kind=RKIND), dimension(:), pointer :: dewpoint_925hPa
       
        real (kind=RKIND), dimension(:), pointer :: uzonal_15hPa
        real (kind=RKIND), dimension(:), pointer :: uzonal_20hPa
        real (kind=RKIND), dimension(:), pointer :: uzonal_30hPa
        real (kind=RKIND), dimension(:), pointer :: uzonal_50hPa
        real (kind=RKIND), dimension(:), pointer :: uzonal_70hPa
        real (kind=RKIND), dimension(:), pointer :: uzonal_100hPa
        real (kind=RKIND), dimension(:), pointer :: uzonal_150hPa
        real (kind=RKIND), dimension(:), pointer :: uzonal_200hPa
        real (kind=RKIND), dimension(:), pointer :: uzonal_250hPa
        real (kind=RKIND), dimension(:), pointer :: uzonal_300hPa
        real (kind=RKIND), dimension(:), pointer :: uzonal_400hPa
        real (kind=RKIND), dimension(:), pointer :: uzonal_500hPa
        real (kind=RKIND), dimension(:), pointer :: uzonal_600hPa
        real (kind=RKIND), dimension(:), pointer :: uzonal_700hPa
        real (kind=RKIND), dimension(:), pointer :: uzonal_825hPa
        real (kind=RKIND), dimension(:), pointer :: uzonal_850hPa
        real (kind=RKIND), dimension(:), pointer :: uzonal_875hPa
        real (kind=RKIND), dimension(:), pointer :: uzonal_900hPa
        real (kind=RKIND), dimension(:), pointer :: uzonal_925hPa
        real (kind=RKIND), dimension(:), pointer :: uzonal_950hPa
        real (kind=RKIND), dimension(:), pointer :: uzonal_975hPa
        real (kind=RKIND), dimension(:), pointer :: uzonal_1000hPa
       
        real (kind=RKIND), dimension(:), pointer :: umeridional_15hPa
        real (kind=RKIND), dimension(:), pointer :: umeridional_20hPa
        real (kind=RKIND), dimension(:), pointer :: umeridional_30hPa
        real (kind=RKIND), dimension(:), pointer :: umeridional_50hPa
        real (kind=RKIND), dimension(:), pointer :: umeridional_70hPa
        real (kind=RKIND), dimension(:), pointer :: umeridional_100hPa
        real (kind=RKIND), dimension(:), pointer :: umeridional_150hPa
        real (kind=RKIND), dimension(:), pointer :: umeridional_200hPa
        real (kind=RKIND), dimension(:), pointer :: umeridional_250hPa
        real (kind=RKIND), dimension(:), pointer :: umeridional_300hPa
        real (kind=RKIND), dimension(:), pointer :: umeridional_400hPa
        real (kind=RKIND), dimension(:), pointer :: umeridional_500hPa
        real (kind=RKIND), dimension(:), pointer :: umeridional_600hPa
        real (kind=RKIND), dimension(:), pointer :: umeridional_700hPa
        real (kind=RKIND), dimension(:), pointer :: umeridional_825hPa
        real (kind=RKIND), dimension(:), pointer :: umeridional_850hPa
        real (kind=RKIND), dimension(:), pointer :: umeridional_875hPa
        real (kind=RKIND), dimension(:), pointer :: umeridional_900hPa
        real (kind=RKIND), dimension(:), pointer :: umeridional_925hPa
        real (kind=RKIND), dimension(:), pointer :: umeridional_950hPa
        real (kind=RKIND), dimension(:), pointer :: umeridional_975hPa
        real (kind=RKIND), dimension(:), pointer :: umeridional_1000hPa        

        real (kind=RKIND), dimension(:), pointer :: spechum_15hPa
        real (kind=RKIND), dimension(:), pointer :: spechum_20hPa
        real (kind=RKIND), dimension(:), pointer :: spechum_30hPa
        real (kind=RKIND), dimension(:), pointer :: spechum_50hPa
        real (kind=RKIND), dimension(:), pointer :: spechum_70hPa
        real (kind=RKIND), dimension(:), pointer :: spechum_100hPa
        real (kind=RKIND), dimension(:), pointer :: spechum_150hPa
        real (kind=RKIND), dimension(:), pointer :: spechum_200hPa
        real (kind=RKIND), dimension(:), pointer :: spechum_250hPa
        real (kind=RKIND), dimension(:), pointer :: spechum_300hPa
        real (kind=RKIND), dimension(:), pointer :: spechum_400hPa
        real (kind=RKIND), dimension(:), pointer :: spechum_500hPa
        real (kind=RKIND), dimension(:), pointer :: spechum_600hPa
        real (kind=RKIND), dimension(:), pointer :: spechum_700hPa
        real (kind=RKIND), dimension(:), pointer :: spechum_825hPa
        real (kind=RKIND), dimension(:), pointer :: spechum_850hPa
        real (kind=RKIND), dimension(:), pointer :: spechum_875hPa
        real (kind=RKIND), dimension(:), pointer :: spechum_900hPa
        real (kind=RKIND), dimension(:), pointer :: spechum_925hPa
        real (kind=RKIND), dimension(:), pointer :: spechum_950hPa
        real (kind=RKIND), dimension(:), pointer :: spechum_975hPa
        real (kind=RKIND), dimension(:), pointer :: spechum_1000hPa

        real (kind=RKIND), dimension(:), pointer :: height_15hPa
        real (kind=RKIND), dimension(:), pointer :: height_20hPa
        real (kind=RKIND), dimension(:), pointer :: height_30hPa
        real (kind=RKIND), dimension(:), pointer :: height_50hPa
        real (kind=RKIND), dimension(:), pointer :: height_70hPa
        real (kind=RKIND), dimension(:), pointer :: height_100hPa
        real (kind=RKIND), dimension(:), pointer :: height_150hPa
        real (kind=RKIND), dimension(:), pointer :: height_200hPa
        real (kind=RKIND), dimension(:), pointer :: height_250hPa
        real (kind=RKIND), dimension(:), pointer :: height_300hPa
        real (kind=RKIND), dimension(:), pointer :: height_400hPa
        real (kind=RKIND), dimension(:), pointer :: height_500hPa
        real (kind=RKIND), dimension(:), pointer :: height_600hPa
        real (kind=RKIND), dimension(:), pointer :: height_700hPa
        real (kind=RKIND), dimension(:), pointer :: height_825hPa
        real (kind=RKIND), dimension(:), pointer :: height_850hPa
        real (kind=RKIND), dimension(:), pointer :: height_875hPa
        real (kind=RKIND), dimension(:), pointer :: height_900hPa
        real (kind=RKIND), dimension(:), pointer :: height_925hPa
        real (kind=RKIND), dimension(:), pointer :: height_950hPa
        real (kind=RKIND), dimension(:), pointer :: height_975hPa
        real (kind=RKIND), dimension(:), pointer :: height_1000hPa
       
        real (kind=RKIND), dimension(:), pointer :: w_15hPa
        real (kind=RKIND), dimension(:), pointer :: w_20hPa
        real (kind=RKIND), dimension(:), pointer :: w_30hPa        
        real (kind=RKIND), dimension(:), pointer :: w_50hPa
        real (kind=RKIND), dimension(:), pointer :: w_70hPa
        real (kind=RKIND), dimension(:), pointer :: w_100hPa
        real (kind=RKIND), dimension(:), pointer :: w_150hPa
        real (kind=RKIND), dimension(:), pointer :: w_200hPa        
        real (kind=RKIND), dimension(:), pointer :: w_250hPa
        real (kind=RKIND), dimension(:), pointer :: w_300hPa
        real (kind=RKIND), dimension(:), pointer :: w_400hPa
        real (kind=RKIND), dimension(:), pointer :: w_500hPa
        real (kind=RKIND), dimension(:), pointer :: w_600hPa
        real (kind=RKIND), dimension(:), pointer :: w_700hPa
        real (kind=RKIND), dimension(:), pointer :: w_825hPa
        real (kind=RKIND), dimension(:), pointer :: w_850hPa
        real (kind=RKIND), dimension(:), pointer :: w_875hPa
        real (kind=RKIND), dimension(:), pointer :: w_900hPa
        real (kind=RKIND), dimension(:), pointer :: w_925hPa
        real (kind=RKIND), dimension(:), pointer :: w_950hPa
        real (kind=RKIND), dimension(:), pointer :: w_975hPa
        real (kind=RKIND), dimension(:), pointer :: w_1000hPa

        real (kind=RKIND), dimension(:), pointer :: geoph_15hPa
        real (kind=RKIND), dimension(:), pointer :: geoph_20hPa
        real (kind=RKIND), dimension(:), pointer :: geoph_30hPa
        real (kind=RKIND), dimension(:), pointer :: geoph_50hPa
        real (kind=RKIND), dimension(:), pointer :: geoph_70hPa
        real (kind=RKIND), dimension(:), pointer :: geoph_100hPa
        real (kind=RKIND), dimension(:), pointer :: geoph_150hPa
        real (kind=RKIND), dimension(:), pointer :: geoph_200hPa
        real (kind=RKIND), dimension(:), pointer :: geoph_250hPa
        real (kind=RKIND), dimension(:), pointer :: geoph_300hPa
        real (kind=RKIND), dimension(:), pointer :: geoph_400hPa
        real (kind=RKIND), dimension(:), pointer :: geoph_500hPa
        real (kind=RKIND), dimension(:), pointer :: geoph_600hPa
        real (kind=RKIND), dimension(:), pointer :: geoph_700hPa
        real (kind=RKIND), dimension(:), pointer :: geoph_825hPa
        real (kind=RKIND), dimension(:), pointer :: geoph_850hPa
        real (kind=RKIND), dimension(:), pointer :: geoph_875hPa
        real (kind=RKIND), dimension(:), pointer :: geoph_900hPa
        real (kind=RKIND), dimension(:), pointer :: geoph_925hPa
        real (kind=RKIND), dimension(:), pointer :: geoph_950hPa
        real (kind=RKIND), dimension(:), pointer :: geoph_975hPa
        real (kind=RKIND), dimension(:), pointer :: geoph_1000hPa

        real (kind=RKIND), dimension(:), pointer :: rqvblten_15hPa
        real (kind=RKIND), dimension(:), pointer :: rqvblten_20hPa
        real (kind=RKIND), dimension(:), pointer :: rqvblten_30hPa
        real (kind=RKIND), dimension(:), pointer :: rqvblten_50hPa
        real (kind=RKIND), dimension(:), pointer :: rqvblten_70hPa
        real (kind=RKIND), dimension(:), pointer :: rqvblten_100hPa
        real (kind=RKIND), dimension(:), pointer :: rqvblten_150hPa
        real (kind=RKIND), dimension(:), pointer :: rqvblten_200hPa
        real (kind=RKIND), dimension(:), pointer :: rqvblten_250hPa
        real (kind=RKIND), dimension(:), pointer :: rqvblten_300hPa
        real (kind=RKIND), dimension(:), pointer :: rqvblten_400hPa
        real (kind=RKIND), dimension(:), pointer :: rqvblten_500hPa
        real (kind=RKIND), dimension(:), pointer :: rqvblten_600hPa
        real (kind=RKIND), dimension(:), pointer :: rqvblten_700hPa
        real (kind=RKIND), dimension(:), pointer :: rqvblten_825hPa
        real (kind=RKIND), dimension(:), pointer :: rqvblten_850hPa
        real (kind=RKIND), dimension(:), pointer :: rqvblten_875hPa
        real (kind=RKIND), dimension(:), pointer :: rqvblten_900hPa
        real (kind=RKIND), dimension(:), pointer :: rqvblten_925hPa
        real (kind=RKIND), dimension(:), pointer :: rqvblten_950hPa
        real (kind=RKIND), dimension(:), pointer :: rqvblten_975hPa
        real (kind=RKIND), dimension(:), pointer :: rqvblten_1000hPa

        real (kind=RKIND), dimension(:), pointer :: rqvcuten_15hPa
        real (kind=RKIND), dimension(:), pointer :: rqvcuten_20hPa
        real (kind=RKIND), dimension(:), pointer :: rqvcuten_30hPa
        real (kind=RKIND), dimension(:), pointer :: rqvcuten_50hPa
        real (kind=RKIND), dimension(:), pointer :: rqvcuten_70hPa
        real (kind=RKIND), dimension(:), pointer :: rqvcuten_100hPa
        real (kind=RKIND), dimension(:), pointer :: rqvcuten_150hPa
        real (kind=RKIND), dimension(:), pointer :: rqvcuten_200hPa
        real (kind=RKIND), dimension(:), pointer :: rqvcuten_250hPa
        real (kind=RKIND), dimension(:), pointer :: rqvcuten_300hPa
        real (kind=RKIND), dimension(:), pointer :: rqvcuten_400hPa
        real (kind=RKIND), dimension(:), pointer :: rqvcuten_500hPa
        real (kind=RKIND), dimension(:), pointer :: rqvcuten_600hPa
        real (kind=RKIND), dimension(:), pointer :: rqvcuten_700hPa
        real (kind=RKIND), dimension(:), pointer :: rqvcuten_825hPa
        real (kind=RKIND), dimension(:), pointer :: rqvcuten_850hPa
        real (kind=RKIND), dimension(:), pointer :: rqvcuten_875hPa
        real (kind=RKIND), dimension(:), pointer :: rqvcuten_900hPa
        real (kind=RKIND), dimension(:), pointer :: rqvcuten_925hPa
        real (kind=RKIND), dimension(:), pointer :: rqvcuten_950hPa
        real (kind=RKIND), dimension(:), pointer :: rqvcuten_975hPa
        real (kind=RKIND), dimension(:), pointer :: rqvcuten_1000hPa

        real (kind=RKIND), dimension(:), pointer :: rthratenlw_15hPa
        real (kind=RKIND), dimension(:), pointer :: rthratenlw_20hPa
        real (kind=RKIND), dimension(:), pointer :: rthratenlw_30hPa
        real (kind=RKIND), dimension(:), pointer :: rthratenlw_50hPa
        real (kind=RKIND), dimension(:), pointer :: rthratenlw_70hPa
        real (kind=RKIND), dimension(:), pointer :: rthratenlw_100hPa
        real (kind=RKIND), dimension(:), pointer :: rthratenlw_150hPa
        real (kind=RKIND), dimension(:), pointer :: rthratenlw_200hPa
        real (kind=RKIND), dimension(:), pointer :: rthratenlw_250hPa
        real (kind=RKIND), dimension(:), pointer :: rthratenlw_300hPa
        real (kind=RKIND), dimension(:), pointer :: rthratenlw_400hPa
        real (kind=RKIND), dimension(:), pointer :: rthratenlw_500hPa
        real (kind=RKIND), dimension(:), pointer :: rthratenlw_600hPa
        real (kind=RKIND), dimension(:), pointer :: rthratenlw_700hPa
        real (kind=RKIND), dimension(:), pointer :: rthratenlw_825hPa
        real (kind=RKIND), dimension(:), pointer :: rthratenlw_850hPa
        real (kind=RKIND), dimension(:), pointer :: rthratenlw_875hPa
        real (kind=RKIND), dimension(:), pointer :: rthratenlw_900hPa
        real (kind=RKIND), dimension(:), pointer :: rthratenlw_925hPa
        real (kind=RKIND), dimension(:), pointer :: rthratenlw_950hPa
        real (kind=RKIND), dimension(:), pointer :: rthratenlw_975hPa
        real (kind=RKIND), dimension(:), pointer :: rthratenlw_1000hPa

        real (kind=RKIND), dimension(:), pointer :: rucuten_15hPa
        real (kind=RKIND), dimension(:), pointer :: rucuten_20hPa
        real (kind=RKIND), dimension(:), pointer :: rucuten_30hPa
        real (kind=RKIND), dimension(:), pointer :: rucuten_50hPa
        real (kind=RKIND), dimension(:), pointer :: rucuten_70hPa
        real (kind=RKIND), dimension(:), pointer :: rucuten_100hPa
        real (kind=RKIND), dimension(:), pointer :: rucuten_150hPa
        real (kind=RKIND), dimension(:), pointer :: rucuten_200hPa
        real (kind=RKIND), dimension(:), pointer :: rucuten_250hPa
        real (kind=RKIND), dimension(:), pointer :: rucuten_300hPa
        real (kind=RKIND), dimension(:), pointer :: rucuten_400hPa
        real (kind=RKIND), dimension(:), pointer :: rucuten_500hPa
        real (kind=RKIND), dimension(:), pointer :: rucuten_600hPa
        real (kind=RKIND), dimension(:), pointer :: rucuten_700hPa
        real (kind=RKIND), dimension(:), pointer :: rucuten_825hPa
        real (kind=RKIND), dimension(:), pointer :: rucuten_850hPa
        real (kind=RKIND), dimension(:), pointer :: rucuten_875hPa
        real (kind=RKIND), dimension(:), pointer :: rucuten_900hPa
        real (kind=RKIND), dimension(:), pointer :: rucuten_925hPa
        real (kind=RKIND), dimension(:), pointer :: rucuten_950hPa
        real (kind=RKIND), dimension(:), pointer :: rucuten_975hPa
        real (kind=RKIND), dimension(:), pointer :: rucuten_1000hPa        

        real (kind=RKIND), dimension(:), pointer :: kzq_15hPa
        real (kind=RKIND), dimension(:), pointer :: kzq_20hPa
        real (kind=RKIND), dimension(:), pointer :: kzq_30hPa
        real (kind=RKIND), dimension(:), pointer :: kzq_50hPa
        real (kind=RKIND), dimension(:), pointer :: kzq_70hPa
        real (kind=RKIND), dimension(:), pointer :: kzq_100hPa
        real (kind=RKIND), dimension(:), pointer :: kzq_150hPa
        real (kind=RKIND), dimension(:), pointer :: kzq_200hPa
        real (kind=RKIND), dimension(:), pointer :: kzq_250hPa
        real (kind=RKIND), dimension(:), pointer :: kzq_300hPa
        real (kind=RKIND), dimension(:), pointer :: kzq_400hPa
        real (kind=RKIND), dimension(:), pointer :: kzq_500hPa
        real (kind=RKIND), dimension(:), pointer :: kzq_600hPa
        real (kind=RKIND), dimension(:), pointer :: kzq_700hPa
        real (kind=RKIND), dimension(:), pointer :: kzq_825hPa
        real (kind=RKIND), dimension(:), pointer :: kzq_850hPa
        real (kind=RKIND), dimension(:), pointer :: kzq_875hPa
        real (kind=RKIND), dimension(:), pointer :: kzq_900hPa
        real (kind=RKIND), dimension(:), pointer :: kzq_925hPa
        real (kind=RKIND), dimension(:), pointer :: kzq_950hPa
        real (kind=RKIND), dimension(:), pointer :: kzq_975hPa
        real (kind=RKIND), dimension(:), pointer :: kzq_1000hPa

        real (kind=RKIND), dimension(:), pointer :: kzh_15hPa
        real (kind=RKIND), dimension(:), pointer :: kzh_20hPa
        real (kind=RKIND), dimension(:), pointer :: kzh_30hPa
        real (kind=RKIND), dimension(:), pointer :: kzh_50hPa
        real (kind=RKIND), dimension(:), pointer :: kzh_70hPa
        real (kind=RKIND), dimension(:), pointer :: kzh_100hPa
        real (kind=RKIND), dimension(:), pointer :: kzh_150hPa
        real (kind=RKIND), dimension(:), pointer :: kzh_200hPa
        real (kind=RKIND), dimension(:), pointer :: kzh_250hPa
        real (kind=RKIND), dimension(:), pointer :: kzh_300hPa
        real (kind=RKIND), dimension(:), pointer :: kzh_400hPa
        real (kind=RKIND), dimension(:), pointer :: kzh_500hPa
        real (kind=RKIND), dimension(:), pointer :: kzh_600hPa
        real (kind=RKIND), dimension(:), pointer :: kzh_700hPa
        real (kind=RKIND), dimension(:), pointer :: kzh_825hPa
        real (kind=RKIND), dimension(:), pointer :: kzh_850hPa
        real (kind=RKIND), dimension(:), pointer :: kzh_875hPa
        real (kind=RKIND), dimension(:), pointer :: kzh_900hPa
        real (kind=RKIND), dimension(:), pointer :: kzh_925hPa
        real (kind=RKIND), dimension(:), pointer :: kzh_950hPa
        real (kind=RKIND), dimension(:), pointer :: kzh_975hPa
        real (kind=RKIND), dimension(:), pointer :: kzh_1000hPa

        real (kind=RKIND), dimension(:), pointer :: kzm_15hPa
        real (kind=RKIND), dimension(:), pointer :: kzm_20hPa
        real (kind=RKIND), dimension(:), pointer :: kzm_30hPa
        real (kind=RKIND), dimension(:), pointer :: kzm_50hPa
        real (kind=RKIND), dimension(:), pointer :: kzm_70hPa
        real (kind=RKIND), dimension(:), pointer :: kzm_100hPa
        real (kind=RKIND), dimension(:), pointer :: kzm_150hPa
        real (kind=RKIND), dimension(:), pointer :: kzm_200hPa
        real (kind=RKIND), dimension(:), pointer :: kzm_250hPa
        real (kind=RKIND), dimension(:), pointer :: kzm_300hPa
        real (kind=RKIND), dimension(:), pointer :: kzm_400hPa
        real (kind=RKIND), dimension(:), pointer :: kzm_500hPa
        real (kind=RKIND), dimension(:), pointer :: kzm_600hPa
        real (kind=RKIND), dimension(:), pointer :: kzm_700hPa
        real (kind=RKIND), dimension(:), pointer :: kzm_825hPa
        real (kind=RKIND), dimension(:), pointer :: kzm_850hPa
        real (kind=RKIND), dimension(:), pointer :: kzm_875hPa
        real (kind=RKIND), dimension(:), pointer :: kzm_900hPa
        real (kind=RKIND), dimension(:), pointer :: kzm_925hPa
        real (kind=RKIND), dimension(:), pointer :: kzm_950hPa
        real (kind=RKIND), dimension(:), pointer :: kzm_975hPa
        real (kind=RKIND), dimension(:), pointer :: kzm_1000hPa

        real (kind=RKIND), dimension(:), pointer :: ni_15hPa
        real (kind=RKIND), dimension(:), pointer :: ni_20hPa
        real (kind=RKIND), dimension(:), pointer :: ni_30hPa
        real (kind=RKIND), dimension(:), pointer :: ni_50hPa
        real (kind=RKIND), dimension(:), pointer :: ni_70hPa
        real (kind=RKIND), dimension(:), pointer :: ni_100hPa
        real (kind=RKIND), dimension(:), pointer :: ni_150hPa
        real (kind=RKIND), dimension(:), pointer :: ni_200hPa
        real (kind=RKIND), dimension(:), pointer :: ni_250hPa
        real (kind=RKIND), dimension(:), pointer :: ni_300hPa
        real (kind=RKIND), dimension(:), pointer :: ni_400hPa
        real (kind=RKIND), dimension(:), pointer :: ni_500hPa
        real (kind=RKIND), dimension(:), pointer :: ni_600hPa
        real (kind=RKIND), dimension(:), pointer :: ni_700hPa
        real (kind=RKIND), dimension(:), pointer :: ni_825hPa
        real (kind=RKIND), dimension(:), pointer :: ni_850hPa
        real (kind=RKIND), dimension(:), pointer :: ni_875hPa
        real (kind=RKIND), dimension(:), pointer :: ni_900hPa
        real (kind=RKIND), dimension(:), pointer :: ni_925hPa
        real (kind=RKIND), dimension(:), pointer :: ni_950hPa
        real (kind=RKIND), dimension(:), pointer :: ni_975hPa
        real (kind=RKIND), dimension(:), pointer :: ni_1000hPa

        real (kind=RKIND), dimension(:), pointer :: nr_15hPa
        real (kind=RKIND), dimension(:), pointer :: nr_20hPa
        real (kind=RKIND), dimension(:), pointer :: nr_30hPa
        real (kind=RKIND), dimension(:), pointer :: nr_50hPa
        real (kind=RKIND), dimension(:), pointer :: nr_70hPa
        real (kind=RKIND), dimension(:), pointer :: nr_100hPa
        real (kind=RKIND), dimension(:), pointer :: nr_150hPa
        real (kind=RKIND), dimension(:), pointer :: nr_200hPa
        real (kind=RKIND), dimension(:), pointer :: nr_250hPa
        real (kind=RKIND), dimension(:), pointer :: nr_300hPa
        real (kind=RKIND), dimension(:), pointer :: nr_400hPa
        real (kind=RKIND), dimension(:), pointer :: nr_500hPa
        real (kind=RKIND), dimension(:), pointer :: nr_600hPa
        real (kind=RKIND), dimension(:), pointer :: nr_700hPa
        real (kind=RKIND), dimension(:), pointer :: nr_825hPa
        real (kind=RKIND), dimension(:), pointer :: nr_850hPa
        real (kind=RKIND), dimension(:), pointer :: nr_875hPa
        real (kind=RKIND), dimension(:), pointer :: nr_900hPa
        real (kind=RKIND), dimension(:), pointer :: nr_925hPa
        real (kind=RKIND), dimension(:), pointer :: nr_950hPa
        real (kind=RKIND), dimension(:), pointer :: nr_975hPa
        real (kind=RKIND), dimension(:), pointer :: nr_1000hPa

        real (kind=RKIND), dimension(:), pointer :: qr_15hPa
        real (kind=RKIND), dimension(:), pointer :: qr_20hPa
        real (kind=RKIND), dimension(:), pointer :: qr_30hPa
        real (kind=RKIND), dimension(:), pointer :: qr_50hPa
        real (kind=RKIND), dimension(:), pointer :: qr_70hPa
        real (kind=RKIND), dimension(:), pointer :: qr_100hPa
        real (kind=RKIND), dimension(:), pointer :: qr_150hPa
        real (kind=RKIND), dimension(:), pointer :: qr_200hPa
        real (kind=RKIND), dimension(:), pointer :: qr_250hPa
        real (kind=RKIND), dimension(:), pointer :: qr_300hPa
        real (kind=RKIND), dimension(:), pointer :: qr_400hPa
        real (kind=RKIND), dimension(:), pointer :: qr_500hPa
        real (kind=RKIND), dimension(:), pointer :: qr_600hPa
        real (kind=RKIND), dimension(:), pointer :: qr_700hPa
        real (kind=RKIND), dimension(:), pointer :: qr_825hPa
        real (kind=RKIND), dimension(:), pointer :: qr_850hPa
        real (kind=RKIND), dimension(:), pointer :: qr_875hPa
        real (kind=RKIND), dimension(:), pointer :: qr_900hPa
        real (kind=RKIND), dimension(:), pointer :: qr_925hPa
        real (kind=RKIND), dimension(:), pointer :: qr_950hPa
        real (kind=RKIND), dimension(:), pointer :: qr_975hPa
        real (kind=RKIND), dimension(:), pointer :: qr_1000hPa

        real (kind=RKIND), dimension(:), pointer :: qs_15hPa
        real (kind=RKIND), dimension(:), pointer :: qs_20hPa
        real (kind=RKIND), dimension(:), pointer :: qs_30hPa
        real (kind=RKIND), dimension(:), pointer :: qs_50hPa
        real (kind=RKIND), dimension(:), pointer :: qs_70hPa
        real (kind=RKIND), dimension(:), pointer :: qs_100hPa
        real (kind=RKIND), dimension(:), pointer :: qs_150hPa
        real (kind=RKIND), dimension(:), pointer :: qs_200hPa
        real (kind=RKIND), dimension(:), pointer :: qs_250hPa
        real (kind=RKIND), dimension(:), pointer :: qs_300hPa
        real (kind=RKIND), dimension(:), pointer :: qs_400hPa
        real (kind=RKIND), dimension(:), pointer :: qs_500hPa
        real (kind=RKIND), dimension(:), pointer :: qs_600hPa
        real (kind=RKIND), dimension(:), pointer :: qs_700hPa
        real (kind=RKIND), dimension(:), pointer :: qs_825hPa
        real (kind=RKIND), dimension(:), pointer :: qs_850hPa
        real (kind=RKIND), dimension(:), pointer :: qs_875hPa
        real (kind=RKIND), dimension(:), pointer :: qs_900hPa
        real (kind=RKIND), dimension(:), pointer :: qs_925hPa
        real (kind=RKIND), dimension(:), pointer :: qs_950hPa
        real (kind=RKIND), dimension(:), pointer :: qs_975hPa
        real (kind=RKIND), dimension(:), pointer :: qs_1000hPa

        real (kind=RKIND), dimension(:), pointer :: qv_15hPa
        real (kind=RKIND), dimension(:), pointer :: qv_20hPa
        real (kind=RKIND), dimension(:), pointer :: qv_30hPa
        real (kind=RKIND), dimension(:), pointer :: qv_50hPa
        real (kind=RKIND), dimension(:), pointer :: qv_70hPa
        real (kind=RKIND), dimension(:), pointer :: qv_100hPa
        real (kind=RKIND), dimension(:), pointer :: qv_150hPa
        real (kind=RKIND), dimension(:), pointer :: qv_200hPa
        real (kind=RKIND), dimension(:), pointer :: qv_250hPa
        real (kind=RKIND), dimension(:), pointer :: qv_300hPa
        real (kind=RKIND), dimension(:), pointer :: qv_400hPa
        real (kind=RKIND), dimension(:), pointer :: qv_500hPa
        real (kind=RKIND), dimension(:), pointer :: qv_600hPa
        real (kind=RKIND), dimension(:), pointer :: qv_700hPa
        real (kind=RKIND), dimension(:), pointer :: qv_825hPa
        real (kind=RKIND), dimension(:), pointer :: qv_850hPa
        real (kind=RKIND), dimension(:), pointer :: qv_875hPa
        real (kind=RKIND), dimension(:), pointer :: qv_900hPa
        real (kind=RKIND), dimension(:), pointer :: qv_925hPa
        real (kind=RKIND), dimension(:), pointer :: qv_950hPa
        real (kind=RKIND), dimension(:), pointer :: qv_975hPa
        real (kind=RKIND), dimension(:), pointer :: qv_1000hPa

        real (kind=RKIND), dimension(:), pointer :: qc_15hPa
        real (kind=RKIND), dimension(:), pointer :: qc_20hPa
        real (kind=RKIND), dimension(:), pointer :: qc_30hPa
        real (kind=RKIND), dimension(:), pointer :: qc_50hPa
        real (kind=RKIND), dimension(:), pointer :: qc_70hPa
        real (kind=RKIND), dimension(:), pointer :: qc_100hPa
        real (kind=RKIND), dimension(:), pointer :: qc_150hPa
        real (kind=RKIND), dimension(:), pointer :: qc_200hPa
        real (kind=RKIND), dimension(:), pointer :: qc_250hPa
        real (kind=RKIND), dimension(:), pointer :: qc_300hPa
        real (kind=RKIND), dimension(:), pointer :: qc_400hPa
        real (kind=RKIND), dimension(:), pointer :: qc_500hPa
        real (kind=RKIND), dimension(:), pointer :: qc_600hPa
        real (kind=RKIND), dimension(:), pointer :: qc_700hPa
        real (kind=RKIND), dimension(:), pointer :: qc_825hPa
        real (kind=RKIND), dimension(:), pointer :: qc_850hPa
        real (kind=RKIND), dimension(:), pointer :: qc_875hPa
        real (kind=RKIND), dimension(:), pointer :: qc_900hPa
        real (kind=RKIND), dimension(:), pointer :: qc_925hPa
        real (kind=RKIND), dimension(:), pointer :: qc_950hPa
        real (kind=RKIND), dimension(:), pointer :: qc_975hPa
        real (kind=RKIND), dimension(:), pointer :: qc_1000hPa

        real (kind=RKIND), dimension(:), pointer :: qg_15hPa
        real (kind=RKIND), dimension(:), pointer :: qg_20hPa
        real (kind=RKIND), dimension(:), pointer :: qg_30hPa
        real (kind=RKIND), dimension(:), pointer :: qg_50hPa
        real (kind=RKIND), dimension(:), pointer :: qg_70hPa
        real (kind=RKIND), dimension(:), pointer :: qg_100hPa
        real (kind=RKIND), dimension(:), pointer :: qg_150hPa
        real (kind=RKIND), dimension(:), pointer :: qg_200hPa
        real (kind=RKIND), dimension(:), pointer :: qg_250hPa
        real (kind=RKIND), dimension(:), pointer :: qg_300hPa
        real (kind=RKIND), dimension(:), pointer :: qg_400hPa
        real (kind=RKIND), dimension(:), pointer :: qg_500hPa
        real (kind=RKIND), dimension(:), pointer :: qg_600hPa
        real (kind=RKIND), dimension(:), pointer :: qg_700hPa
        real (kind=RKIND), dimension(:), pointer :: qg_825hPa
        real (kind=RKIND), dimension(:), pointer :: qg_850hPa
        real (kind=RKIND), dimension(:), pointer :: qg_875hPa
        real (kind=RKIND), dimension(:), pointer :: qg_900hPa
        real (kind=RKIND), dimension(:), pointer :: qg_925hPa
        real (kind=RKIND), dimension(:), pointer :: qg_950hPa
        real (kind=RKIND), dimension(:), pointer :: qg_975hPa
        real (kind=RKIND), dimension(:), pointer :: qg_1000hPa

        real (kind=RKIND), dimension(:), pointer :: qi_15hPa
        real (kind=RKIND), dimension(:), pointer :: qi_20hPa
        real (kind=RKIND), dimension(:), pointer :: qi_30hPa
        real (kind=RKIND), dimension(:), pointer :: qi_50hPa
        real (kind=RKIND), dimension(:), pointer :: qi_70hPa
        real (kind=RKIND), dimension(:), pointer :: qi_100hPa
        real (kind=RKIND), dimension(:), pointer :: qi_150hPa
        real (kind=RKIND), dimension(:), pointer :: qi_200hPa
        real (kind=RKIND), dimension(:), pointer :: qi_250hPa
        real (kind=RKIND), dimension(:), pointer :: qi_300hPa
        real (kind=RKIND), dimension(:), pointer :: qi_400hPa
        real (kind=RKIND), dimension(:), pointer :: qi_500hPa
        real (kind=RKIND), dimension(:), pointer :: qi_600hPa
        real (kind=RKIND), dimension(:), pointer :: qi_700hPa
        real (kind=RKIND), dimension(:), pointer :: qi_825hPa
        real (kind=RKIND), dimension(:), pointer :: qi_850hPa
        real (kind=RKIND), dimension(:), pointer :: qi_875hPa
        real (kind=RKIND), dimension(:), pointer :: qi_900hPa
        real (kind=RKIND), dimension(:), pointer :: qi_925hPa
        real (kind=RKIND), dimension(:), pointer :: qi_950hPa
        real (kind=RKIND), dimension(:), pointer :: qi_975hPa
        real (kind=RKIND), dimension(:), pointer :: qi_1000hPa

        real (kind=RKIND), dimension(:), pointer :: vorticity_50hPa
        real (kind=RKIND), dimension(:), pointer :: vorticity_100hPa
        real (kind=RKIND), dimension(:), pointer :: vorticity_200hPa
        real (kind=RKIND), dimension(:), pointer :: vorticity_250hPa
        real (kind=RKIND), dimension(:), pointer :: vorticity_500hPa
        real (kind=RKIND), dimension(:), pointer :: vorticity_700hPa
        real (kind=RKIND), dimension(:), pointer :: vorticity_850hPa
        real (kind=RKIND), dimension(:), pointer :: vorticity_925hPa
       
        real (kind=RKIND) :: evp
       
       !--------------------
       
        real (kind=RKIND), dimension(:), pointer :: mslp
       
        real (kind=RKIND), dimension(:,:), allocatable :: pressure, pressureCp1, pressure2, pressure_v, temperature
        real (kind=RKIND), dimension(:,:), allocatable :: dewpoint, spechum, geoph
       
       !local interpolated fields:
        integer :: nIntP
        real (kind=RKIND) :: w1,w2,z0,z1,z2
        real (kind=RKIND), dimension(:,:), allocatable :: field_in,press_in
        real (kind=RKIND), dimension(:,:), allocatable :: field_interp,press_interp
        
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
        call mpas_pool_get_array(tend_physics, 'rqvblten', rqvblten)
        call mpas_pool_get_array(tend_physics, 'rthratenlw', rthratenlw)
        call mpas_pool_get_array(tend_physics, 'rucuten', rucuten)

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
       
        call mpas_pool_get_array(diag, 't_iso_levels', t_iso_levels)
        call mpas_pool_get_array(diag, 'z_iso_levels', z_iso_levels)
        call mpas_pool_get_array(diag, 't_isobaric', t_isobaric)
        call mpas_pool_get_array(diag, 'z_isobaric', z_isobaric)
        call mpas_pool_get_array(diag, 'meanT_500_300', meanT_500_300)
      
        call mpas_pool_get_array(diag, 'cldfrac_15hPa'  , cldfrac_15hPa)
        call mpas_pool_get_array(diag, 'cldfrac_20hPa'  , cldfrac_20hPa)
        call mpas_pool_get_array(diag, 'cldfrac_30hPa'  , cldfrac_30hPa)
        call mpas_pool_get_array(diag, 'cldfrac_50hPa'  , cldfrac_50hPa)
        call mpas_pool_get_array(diag, 'cldfrac_70hPa'  , cldfrac_70hPa)
        call mpas_pool_get_array(diag, 'cldfrac_100hPa' , cldfrac_100hPa)
        call mpas_pool_get_array(diag, 'cldfrac_150hPa' , cldfrac_150hPa)
        call mpas_pool_get_array(diag, 'cldfrac_200hPa' , cldfrac_200hPa)
        call mpas_pool_get_array(diag, 'cldfrac_250hPa' , cldfrac_250hPa)
        call mpas_pool_get_array(diag, 'cldfrac_300hPa' , cldfrac_300hPa)
        call mpas_pool_get_array(diag, 'cldfrac_400hPa' , cldfrac_400hPa)
        call mpas_pool_get_array(diag, 'cldfrac_500hPa' , cldfrac_500hPa)
        call mpas_pool_get_array(diag, 'cldfrac_600hPa' , cldfrac_600hPa)
        call mpas_pool_get_array(diag, 'cldfrac_700hPa' , cldfrac_700hPa)
        call mpas_pool_get_array(diag, 'cldfrac_825hPa' , cldfrac_825hPa)
        call mpas_pool_get_array(diag, 'cldfrac_850hPa' , cldfrac_850hPa)
        call mpas_pool_get_array(diag, 'cldfrac_875hPa' , cldfrac_875hPa)
        call mpas_pool_get_array(diag, 'cldfrac_900hPa' , cldfrac_900hPa)
        call mpas_pool_get_array(diag, 'cldfrac_925hPa' , cldfrac_925hPa)
        call mpas_pool_get_array(diag, 'cldfrac_950hPa' , cldfrac_950hPa)
        call mpas_pool_get_array(diag, 'cldfrac_975hPa' , cldfrac_975hPa)
        call mpas_pool_get_array(diag, 'cldfrac_1000hPa', cldfrac_1000hPa)


        call mpas_pool_get_array(diag, 'temperature_15hPa'  , temperature_15hPa)
        call mpas_pool_get_array(diag, 'temperature_20hPa'  , temperature_20hPa)
        call mpas_pool_get_array(diag, 'temperature_30hPa'  , temperature_30hPa)
        call mpas_pool_get_array(diag, 'temperature_50hPa'  , temperature_50hPa)
        call mpas_pool_get_array(diag, 'temperature_70hPa'  , temperature_70hPa)
        call mpas_pool_get_array(diag, 'temperature_100hPa' , temperature_100hPa)
        call mpas_pool_get_array(diag, 'temperature_150hPa' , temperature_150hPa)
        call mpas_pool_get_array(diag, 'temperature_200hPa' , temperature_200hPa)
        call mpas_pool_get_array(diag, 'temperature_250hPa' , temperature_250hPa)
        call mpas_pool_get_array(diag, 'temperature_300hPa' , temperature_300hPa)
        call mpas_pool_get_array(diag, 'temperature_400hPa' , temperature_400hPa)
        call mpas_pool_get_array(diag, 'temperature_500hPa' , temperature_500hPa)
        call mpas_pool_get_array(diag, 'temperature_600hPa' , temperature_600hPa)
        call mpas_pool_get_array(diag, 'temperature_700hPa' , temperature_700hPa)
        call mpas_pool_get_array(diag, 'temperature_825hPa' , temperature_825hPa)
        call mpas_pool_get_array(diag, 'temperature_850hPa' , temperature_850hPa)
        call mpas_pool_get_array(diag, 'temperature_875hPa' , temperature_875hPa)
        call mpas_pool_get_array(diag, 'temperature_900hPa' , temperature_900hPa)
        call mpas_pool_get_array(diag, 'temperature_925hPa' , temperature_925hPa)
        call mpas_pool_get_array(diag, 'temperature_950hPa' , temperature_950hPa)
        call mpas_pool_get_array(diag, 'temperature_975hPa' , temperature_975hPa)
        call mpas_pool_get_array(diag, 'temperature_1000hPa', temperature_1000hPa)        

        call mpas_pool_get_array(diag, 'relhum_15hPa', relhum_15hPa)
        call mpas_pool_get_array(diag, 'relhum_20hPa', relhum_20hPa)
        call mpas_pool_get_array(diag, 'relhum_30hPa', relhum_30hPa)
        call mpas_pool_get_array(diag, 'relhum_50hPa', relhum_50hPa)
        call mpas_pool_get_array(diag, 'relhum_70hPa', relhum_70hPa)
        call mpas_pool_get_array(diag, 'relhum_100hPa', relhum_100hPa)
        call mpas_pool_get_array(diag, 'relhum_150hPa', relhum_150hPa)
        call mpas_pool_get_array(diag, 'relhum_200hPa', relhum_200hPa)
        call mpas_pool_get_array(diag, 'relhum_250hPa', relhum_250hPa)
        call mpas_pool_get_array(diag, 'relhum_300hPa', relhum_300hPa)
        call mpas_pool_get_array(diag, 'relhum_400hPa', relhum_400hPa)
        call mpas_pool_get_array(diag, 'relhum_500hPa', relhum_500hPa)
        call mpas_pool_get_array(diag, 'relhum_600hPa', relhum_600hPa)
        call mpas_pool_get_array(diag, 'relhum_700hPa', relhum_700hPa)
        call mpas_pool_get_array(diag, 'relhum_825hPa', relhum_825hPa)
        call mpas_pool_get_array(diag, 'relhum_850hPa', relhum_850hPa)
        call mpas_pool_get_array(diag, 'relhum_875hPa', relhum_875hPa)
        call mpas_pool_get_array(diag, 'relhum_900hPa', relhum_900hPa)
        call mpas_pool_get_array(diag, 'relhum_925hPa', relhum_925hPa)
        call mpas_pool_get_array(diag, 'relhum_950hPa', relhum_950hPa)
        call mpas_pool_get_array(diag, 'relhum_975hPa', relhum_975hPa)
        call mpas_pool_get_array(diag, 'relhum_1000hPa', relhum_1000hPa)        
       
        call mpas_pool_get_array(diag, 'dewpoint_50hPa', dewpoint_50hPa)
        call mpas_pool_get_array(diag, 'dewpoint_100hPa', dewpoint_100hPa)
        call mpas_pool_get_array(diag, 'dewpoint_200hPa', dewpoint_200hPa)
        call mpas_pool_get_array(diag, 'dewpoint_250hPa', dewpoint_250hPa)
        call mpas_pool_get_array(diag, 'dewpoint_500hPa', dewpoint_500hPa)
        call mpas_pool_get_array(diag, 'dewpoint_700hPa', dewpoint_700hPa)
        call mpas_pool_get_array(diag, 'dewpoint_850hPa', dewpoint_850hPa)
        call mpas_pool_get_array(diag, 'dewpoint_925hPa', dewpoint_925hPa)
       
        call mpas_pool_get_array(diag, 'uzonal_15hPa'   , uzonal_15hPa)
        call mpas_pool_get_array(diag, 'uzonal_20hPa'   , uzonal_20hPa)
        call mpas_pool_get_array(diag, 'uzonal_30hPa'   , uzonal_30hPa)
        call mpas_pool_get_array(diag, 'uzonal_50hPa'   , uzonal_50hPa)
        call mpas_pool_get_array(diag, 'uzonal_70hPa'   , uzonal_70hPa)
        call mpas_pool_get_array(diag, 'uzonal_100hPa'  , uzonal_100hPa)
        call mpas_pool_get_array(diag, 'uzonal_150hPa'  , uzonal_150hPa)
        call mpas_pool_get_array(diag, 'uzonal_200hPa'  , uzonal_200hPa)
        call mpas_pool_get_array(diag, 'uzonal_250hPa'  , uzonal_250hPa)
        call mpas_pool_get_array(diag, 'uzonal_300hPa'  , uzonal_300hPa)
        call mpas_pool_get_array(diag, 'uzonal_400hPa'  , uzonal_400hPa)
        call mpas_pool_get_array(diag, 'uzonal_500hPa'  , uzonal_500hPa)
        call mpas_pool_get_array(diag, 'uzonal_600hPa'  , uzonal_600hPa)
        call mpas_pool_get_array(diag, 'uzonal_700hPa'  , uzonal_700hPa)
        call mpas_pool_get_array(diag, 'uzonal_825hPa'  , uzonal_825hPa)
        call mpas_pool_get_array(diag, 'uzonal_850hPa'  , uzonal_850hPa)
        call mpas_pool_get_array(diag, 'uzonal_875hPa'  , uzonal_875hPa)
        call mpas_pool_get_array(diag, 'uzonal_900hPa'  , uzonal_900hPa)
        call mpas_pool_get_array(diag, 'uzonal_925hPa'  , uzonal_925hPa)
        call mpas_pool_get_array(diag, 'uzonal_950hPa'  , uzonal_950hPa)
        call mpas_pool_get_array(diag, 'uzonal_975hPa'  , uzonal_975hPa)
        call mpas_pool_get_array(diag, 'uzonal_1000hPa' , uzonal_1000hPa)
       
        call mpas_pool_get_array(diag, 'umeridional_15hPa', umeridional_15hPa)
        call mpas_pool_get_array(diag, 'umeridional_20hPa', umeridional_20hPa)
        call mpas_pool_get_array(diag, 'umeridional_30hPa', umeridional_30hPa)
        call mpas_pool_get_array(diag, 'umeridional_50hPa', umeridional_50hPa)
        call mpas_pool_get_array(diag, 'umeridional_70hPa', umeridional_70hPa)
        call mpas_pool_get_array(diag, 'umeridional_100hPa', umeridional_100hPa)
        call mpas_pool_get_array(diag, 'umeridional_150hPa', umeridional_150hPa)
        call mpas_pool_get_array(diag, 'umeridional_200hPa', umeridional_200hPa)
        call mpas_pool_get_array(diag, 'umeridional_250hPa', umeridional_250hPa)
        call mpas_pool_get_array(diag, 'umeridional_300hPa', umeridional_300hPa)
        call mpas_pool_get_array(diag, 'umeridional_400hPa', umeridional_400hPa)
        call mpas_pool_get_array(diag, 'umeridional_500hPa', umeridional_500hPa)
        call mpas_pool_get_array(diag, 'umeridional_600hPa', umeridional_600hPa)
        call mpas_pool_get_array(diag, 'umeridional_700hPa', umeridional_700hPa)
        call mpas_pool_get_array(diag, 'umeridional_825hPa', umeridional_825hPa)
        call mpas_pool_get_array(diag, 'umeridional_850hPa', umeridional_850hPa)
        call mpas_pool_get_array(diag, 'umeridional_875hPa', umeridional_875hPa)
        call mpas_pool_get_array(diag, 'umeridional_900hPa', umeridional_900hPa)
        call mpas_pool_get_array(diag, 'umeridional_925hPa', umeridional_925hPa)
        call mpas_pool_get_array(diag, 'umeridional_950hPa', umeridional_950hPa)
        call mpas_pool_get_array(diag, 'umeridional_975hPa', umeridional_975hPa)
        call mpas_pool_get_array(diag, 'umeridional_1000hPa', umeridional_1000hPa)
       
        call mpas_pool_get_array(diag, 'spechum_15hPa', spechum_15hPa)
        call mpas_pool_get_array(diag, 'spechum_20hPa', spechum_20hPa)
        call mpas_pool_get_array(diag, 'spechum_30hPa', spechum_30hPa)
        call mpas_pool_get_array(diag, 'spechum_50hPa', spechum_50hPa)
        call mpas_pool_get_array(diag, 'spechum_70hPa', spechum_70hPa)
        call mpas_pool_get_array(diag, 'spechum_100hPa', spechum_100hPa)
        call mpas_pool_get_array(diag, 'spechum_150hPa', spechum_150hPa)
        call mpas_pool_get_array(diag, 'spechum_200hPa', spechum_200hPa)
        call mpas_pool_get_array(diag, 'spechum_250hPa', spechum_250hPa)
        call mpas_pool_get_array(diag, 'spechum_300hPa', spechum_300hPa)
        call mpas_pool_get_array(diag, 'spechum_400hPa', spechum_400hPa)
        call mpas_pool_get_array(diag, 'spechum_500hPa', spechum_500hPa)
        call mpas_pool_get_array(diag, 'spechum_600hPa', spechum_600hPa)
        call mpas_pool_get_array(diag, 'spechum_700hPa', spechum_700hPa)
        call mpas_pool_get_array(diag, 'spechum_825hPa', spechum_825hPa)
        call mpas_pool_get_array(diag, 'spechum_850hPa', spechum_850hPa)
        call mpas_pool_get_array(diag, 'spechum_875hPa', spechum_875hPa)
        call mpas_pool_get_array(diag, 'spechum_900hPa', spechum_900hPa)
        call mpas_pool_get_array(diag, 'spechum_925hPa', spechum_925hPa)
        call mpas_pool_get_array(diag, 'spechum_950hPa', spechum_950hPa)
        call mpas_pool_get_array(diag, 'spechum_975hPa', spechum_975hPa)
        call mpas_pool_get_array(diag, 'spechum_1000hPa', spechum_1000hPa)

        call mpas_pool_get_array(diag, 'height_15hPa', height_15hPa)
        call mpas_pool_get_array(diag, 'height_20hPa', height_20hPa)
        call mpas_pool_get_array(diag, 'height_30hPa', height_30hPa)
        call mpas_pool_get_array(diag, 'height_50hPa', height_50hPa)
        call mpas_pool_get_array(diag, 'height_70hPa', height_70hPa)
        call mpas_pool_get_array(diag, 'height_100hPa', height_100hPa)
        call mpas_pool_get_array(diag, 'height_150hPa', height_150hPa)
        call mpas_pool_get_array(diag, 'height_200hPa', height_200hPa)
        call mpas_pool_get_array(diag, 'height_250hPa', height_250hPa)
        call mpas_pool_get_array(diag, 'height_300hPa', height_300hPa)
        call mpas_pool_get_array(diag, 'height_400hPa', height_400hPa)
        call mpas_pool_get_array(diag, 'height_500hPa', height_500hPa)
        call mpas_pool_get_array(diag, 'height_600hPa', height_600hPa)
        call mpas_pool_get_array(diag, 'height_700hPa', height_700hPa)
        call mpas_pool_get_array(diag, 'height_825hPa', height_825hPa)
        call mpas_pool_get_array(diag, 'height_850hPa', height_850hPa)
        call mpas_pool_get_array(diag, 'height_875hPa', height_875hPa)
        call mpas_pool_get_array(diag, 'height_900hPa', height_900hPa)
        call mpas_pool_get_array(diag, 'height_925hPa', height_925hPa)
        call mpas_pool_get_array(diag, 'height_950hPa', height_950hPa)
        call mpas_pool_get_array(diag, 'height_975hPa', height_975hPa)
        call mpas_pool_get_array(diag, 'height_1000hPa', height_1000hPa)
       
        call mpas_pool_get_array(diag, 'w_15hPa', w_15hPa)
        call mpas_pool_get_array(diag, 'w_20hPa', w_20hPa)
        call mpas_pool_get_array(diag, 'w_30hPa', w_30hPa)                
        call mpas_pool_get_array(diag, 'w_50hPa', w_50hPa)
        call mpas_pool_get_array(diag, 'w_70hPa', w_70hPa)
        call mpas_pool_get_array(diag, 'w_100hPa', w_100hPa)
        call mpas_pool_get_array(diag, 'w_150hPa', w_150hPa)
        call mpas_pool_get_array(diag, 'w_200hPa', w_200hPa)        
        call mpas_pool_get_array(diag, 'w_250hPa', w_250hPa)
        call mpas_pool_get_array(diag, 'w_300hPa', w_300hPa)
        call mpas_pool_get_array(diag, 'w_400hPa', w_400hPa)
        call mpas_pool_get_array(diag, 'w_500hPa', w_500hPa)
        call mpas_pool_get_array(diag, 'w_600hPa', w_600hPa)
        call mpas_pool_get_array(diag, 'w_700hPa', w_700hPa)
        call mpas_pool_get_array(diag, 'w_825hPa', w_825hPa)
        call mpas_pool_get_array(diag, 'w_850hPa', w_850hPa)
        call mpas_pool_get_array(diag, 'w_875hPa', w_875hPa)
        call mpas_pool_get_array(diag, 'w_900hPa', w_900hPa)
        call mpas_pool_get_array(diag, 'w_925hPa', w_925hPa)
        call mpas_pool_get_array(diag, 'w_950hPa', w_950hPa)
        call mpas_pool_get_array(diag, 'w_975hPa', w_975hPa)
        call mpas_pool_get_array(diag, 'w_1000hPa', w_1000hPa)
       
        call mpas_pool_get_array(diag, 'geoph_15hPa', geoph_15hPa)
        call mpas_pool_get_array(diag, 'geoph_20hPa', geoph_20hPa)
        call mpas_pool_get_array(diag, 'geoph_30hPa', geoph_30hPa)
        call mpas_pool_get_array(diag, 'geoph_50hPa', geoph_50hPa)
        call mpas_pool_get_array(diag, 'geoph_70hPa', geoph_70hPa)
        call mpas_pool_get_array(diag, 'geoph_100hPa', geoph_100hPa)
        call mpas_pool_get_array(diag, 'geoph_150hPa', geoph_150hPa)
        call mpas_pool_get_array(diag, 'geoph_200hPa', geoph_200hPa)
        call mpas_pool_get_array(diag, 'geoph_250hPa', geoph_250hPa)
        call mpas_pool_get_array(diag, 'geoph_300hPa', geoph_300hPa)
        call mpas_pool_get_array(diag, 'geoph_400hPa', geoph_400hPa)
        call mpas_pool_get_array(diag, 'geoph_500hPa', geoph_500hPa)
        call mpas_pool_get_array(diag, 'geoph_600hPa', geoph_600hPa)
        call mpas_pool_get_array(diag, 'geoph_700hPa', geoph_700hPa)
        call mpas_pool_get_array(diag, 'geoph_825hPa', geoph_825hPa)
        call mpas_pool_get_array(diag, 'geoph_850hPa', geoph_850hPa)
        call mpas_pool_get_array(diag, 'geoph_875hPa', geoph_875hPa)
        call mpas_pool_get_array(diag, 'geoph_900hPa', geoph_900hPa)
        call mpas_pool_get_array(diag, 'geoph_925hPa', geoph_925hPa)
        call mpas_pool_get_array(diag, 'geoph_950hPa', geoph_950hPa)
        call mpas_pool_get_array(diag, 'geoph_975hPa', geoph_975hPa)
        call mpas_pool_get_array(diag, 'geoph_1000hPa', geoph_1000hPa)


        call mpas_pool_get_array(diag, 'rqvblten_15hPa', rqvblten_15hPa)
        call mpas_pool_get_array(diag, 'rqvblten_20hPa', rqvblten_20hPa)
        call mpas_pool_get_array(diag, 'rqvblten_30hPa', rqvblten_30hPa)
        call mpas_pool_get_array(diag, 'rqvblten_50hPa', rqvblten_50hPa)
        call mpas_pool_get_array(diag, 'rqvblten_70hPa', rqvblten_70hPa)
        call mpas_pool_get_array(diag, 'rqvblten_100hPa', rqvblten_100hPa)
        call mpas_pool_get_array(diag, 'rqvblten_150hPa', rqvblten_150hPa)
        call mpas_pool_get_array(diag, 'rqvblten_200hPa', rqvblten_200hPa)
        call mpas_pool_get_array(diag, 'rqvblten_250hPa', rqvblten_250hPa)
        call mpas_pool_get_array(diag, 'rqvblten_300hPa', rqvblten_300hPa)
        call mpas_pool_get_array(diag, 'rqvblten_400hPa', rqvblten_400hPa)
        call mpas_pool_get_array(diag, 'rqvblten_500hPa', rqvblten_500hPa)
        call mpas_pool_get_array(diag, 'rqvblten_600hPa', rqvblten_600hPa)
        call mpas_pool_get_array(diag, 'rqvblten_700hPa', rqvblten_700hPa)
        call mpas_pool_get_array(diag, 'rqvblten_825hPa', rqvblten_825hPa)
        call mpas_pool_get_array(diag, 'rqvblten_850hPa', rqvblten_850hPa)
        call mpas_pool_get_array(diag, 'rqvblten_875hPa', rqvblten_875hPa)
        call mpas_pool_get_array(diag, 'rqvblten_900hPa', rqvblten_900hPa)
        call mpas_pool_get_array(diag, 'rqvblten_925hPa', rqvblten_925hPa)
        call mpas_pool_get_array(diag, 'rqvblten_950hPa', rqvblten_950hPa)
        call mpas_pool_get_array(diag, 'rqvblten_975hPa', rqvblten_975hPa)
        call mpas_pool_get_array(diag, 'rqvblten_1000hPa', rqvblten_1000hPa)

        call mpas_pool_get_array(diag, 'rqvcuten_15hPa', rqvcuten_15hPa)
        call mpas_pool_get_array(diag, 'rqvcuten_20hPa', rqvcuten_20hPa)
        call mpas_pool_get_array(diag, 'rqvcuten_30hPa', rqvcuten_30hPa)
        call mpas_pool_get_array(diag, 'rqvcuten_50hPa', rqvcuten_50hPa)
        call mpas_pool_get_array(diag, 'rqvcuten_70hPa', rqvcuten_70hPa)
        call mpas_pool_get_array(diag, 'rqvcuten_100hPa', rqvcuten_100hPa)
        call mpas_pool_get_array(diag, 'rqvcuten_150hPa', rqvcuten_150hPa)
        call mpas_pool_get_array(diag, 'rqvcuten_200hPa', rqvcuten_200hPa)
        call mpas_pool_get_array(diag, 'rqvcuten_250hPa', rqvcuten_250hPa)
        call mpas_pool_get_array(diag, 'rqvcuten_300hPa', rqvcuten_300hPa)
        call mpas_pool_get_array(diag, 'rqvcuten_400hPa', rqvcuten_400hPa)
        call mpas_pool_get_array(diag, 'rqvcuten_500hPa', rqvcuten_500hPa)
        call mpas_pool_get_array(diag, 'rqvcuten_600hPa', rqvcuten_600hPa)
        call mpas_pool_get_array(diag, 'rqvcuten_700hPa', rqvcuten_700hPa)
        call mpas_pool_get_array(diag, 'rqvcuten_825hPa', rqvcuten_825hPa)
        call mpas_pool_get_array(diag, 'rqvcuten_850hPa', rqvcuten_850hPa)
        call mpas_pool_get_array(diag, 'rqvcuten_875hPa', rqvcuten_875hPa)
        call mpas_pool_get_array(diag, 'rqvcuten_900hPa', rqvcuten_900hPa)
        call mpas_pool_get_array(diag, 'rqvcuten_925hPa', rqvcuten_925hPa)
        call mpas_pool_get_array(diag, 'rqvcuten_950hPa', rqvcuten_950hPa)
        call mpas_pool_get_array(diag, 'rqvcuten_975hPa', rqvcuten_975hPa)
        call mpas_pool_get_array(diag, 'rqvcuten_1000hPa', rqvcuten_1000hPa)

        call mpas_pool_get_array(diag, 'rthratenlw_15hPa', rthratenlw_15hPa)
        call mpas_pool_get_array(diag, 'rthratenlw_20hPa', rthratenlw_20hPa)
        call mpas_pool_get_array(diag, 'rthratenlw_30hPa', rthratenlw_30hPa)
        call mpas_pool_get_array(diag, 'rthratenlw_50hPa', rthratenlw_50hPa)
        call mpas_pool_get_array(diag, 'rthratenlw_70hPa', rthratenlw_70hPa)
        call mpas_pool_get_array(diag, 'rthratenlw_100hPa', rthratenlw_100hPa)
        call mpas_pool_get_array(diag, 'rthratenlw_150hPa', rthratenlw_150hPa)
        call mpas_pool_get_array(diag, 'rthratenlw_200hPa', rthratenlw_200hPa)
        call mpas_pool_get_array(diag, 'rthratenlw_250hPa', rthratenlw_250hPa)
        call mpas_pool_get_array(diag, 'rthratenlw_300hPa', rthratenlw_300hPa)
        call mpas_pool_get_array(diag, 'rthratenlw_400hPa', rthratenlw_400hPa)
        call mpas_pool_get_array(diag, 'rthratenlw_500hPa', rthratenlw_500hPa)
        call mpas_pool_get_array(diag, 'rthratenlw_600hPa', rthratenlw_600hPa)
        call mpas_pool_get_array(diag, 'rthratenlw_700hPa', rthratenlw_700hPa)
        call mpas_pool_get_array(diag, 'rthratenlw_825hPa', rthratenlw_825hPa)
        call mpas_pool_get_array(diag, 'rthratenlw_850hPa', rthratenlw_850hPa)
        call mpas_pool_get_array(diag, 'rthratenlw_875hPa', rthratenlw_875hPa)
        call mpas_pool_get_array(diag, 'rthratenlw_900hPa', rthratenlw_900hPa)
        call mpas_pool_get_array(diag, 'rthratenlw_925hPa', rthratenlw_925hPa)
        call mpas_pool_get_array(diag, 'rthratenlw_950hPa', rthratenlw_950hPa)
        call mpas_pool_get_array(diag, 'rthratenlw_975hPa', rthratenlw_975hPa)
        call mpas_pool_get_array(diag, 'rthratenlw_1000hPa', rthratenlw_1000hPa)

        call mpas_pool_get_array(diag, 'rucuten_15hPa', rucuten_15hPa)
        call mpas_pool_get_array(diag, 'rucuten_20hPa', rucuten_20hPa)
        call mpas_pool_get_array(diag, 'rucuten_30hPa', rucuten_30hPa)
        call mpas_pool_get_array(diag, 'rucuten_50hPa', rucuten_50hPa)
        call mpas_pool_get_array(diag, 'rucuten_70hPa', rucuten_70hPa)
        call mpas_pool_get_array(diag, 'rucuten_100hPa', rucuten_100hPa)
        call mpas_pool_get_array(diag, 'rucuten_150hPa', rucuten_150hPa)
        call mpas_pool_get_array(diag, 'rucuten_200hPa', rucuten_200hPa)
        call mpas_pool_get_array(diag, 'rucuten_250hPa', rucuten_250hPa)
        call mpas_pool_get_array(diag, 'rucuten_300hPa', rucuten_300hPa)
        call mpas_pool_get_array(diag, 'rucuten_400hPa', rucuten_400hPa)
        call mpas_pool_get_array(diag, 'rucuten_500hPa', rucuten_500hPa)
        call mpas_pool_get_array(diag, 'rucuten_600hPa', rucuten_600hPa)
        call mpas_pool_get_array(diag, 'rucuten_700hPa', rucuten_700hPa)
        call mpas_pool_get_array(diag, 'rucuten_825hPa', rucuten_825hPa)
        call mpas_pool_get_array(diag, 'rucuten_850hPa', rucuten_850hPa)
        call mpas_pool_get_array(diag, 'rucuten_875hPa', rucuten_875hPa)
        call mpas_pool_get_array(diag, 'rucuten_900hPa', rucuten_900hPa)
        call mpas_pool_get_array(diag, 'rucuten_925hPa', rucuten_925hPa)
        call mpas_pool_get_array(diag, 'rucuten_950hPa', rucuten_950hPa)
        call mpas_pool_get_array(diag, 'rucuten_975hPa', rucuten_975hPa)
        call mpas_pool_get_array(diag, 'rucuten_1000hPa', rucuten_1000hPa)        

        call mpas_pool_get_array(diag, 'kzq_15hPa', kzq_15hPa)
        call mpas_pool_get_array(diag, 'kzq_20hPa', kzq_20hPa)
        call mpas_pool_get_array(diag, 'kzq_30hPa', kzq_30hPa)
        call mpas_pool_get_array(diag, 'kzq_50hPa', kzq_50hPa)
        call mpas_pool_get_array(diag, 'kzq_70hPa', kzq_70hPa)
        call mpas_pool_get_array(diag, 'kzq_100hPa', kzq_100hPa)
        call mpas_pool_get_array(diag, 'kzq_150hPa', kzq_150hPa)
        call mpas_pool_get_array(diag, 'kzq_200hPa', kzq_200hPa)
        call mpas_pool_get_array(diag, 'kzq_250hPa', kzq_250hPa)
        call mpas_pool_get_array(diag, 'kzq_300hPa', kzq_300hPa)
        call mpas_pool_get_array(diag, 'kzq_400hPa', kzq_400hPa)
        call mpas_pool_get_array(diag, 'kzq_500hPa', kzq_500hPa)
        call mpas_pool_get_array(diag, 'kzq_600hPa', kzq_600hPa)
        call mpas_pool_get_array(diag, 'kzq_700hPa', kzq_700hPa)
        call mpas_pool_get_array(diag, 'kzq_825hPa', kzq_825hPa)
        call mpas_pool_get_array(diag, 'kzq_850hPa', kzq_850hPa)
        call mpas_pool_get_array(diag, 'kzq_875hPa', kzq_875hPa)
        call mpas_pool_get_array(diag, 'kzq_900hPa', kzq_900hPa)
        call mpas_pool_get_array(diag, 'kzq_925hPa', kzq_925hPa)
        call mpas_pool_get_array(diag, 'kzq_950hPa', kzq_950hPa)
        call mpas_pool_get_array(diag, 'kzq_975hPa', kzq_975hPa)
        call mpas_pool_get_array(diag, 'kzq_1000hPa', kzq_1000hPa)

        call mpas_pool_get_array(diag, 'kzh_15hPa', kzh_15hPa)
        call mpas_pool_get_array(diag, 'kzh_20hPa', kzh_20hPa)
        call mpas_pool_get_array(diag, 'kzh_30hPa', kzh_30hPa)
        call mpas_pool_get_array(diag, 'kzh_50hPa', kzh_50hPa)
        call mpas_pool_get_array(diag, 'kzh_70hPa', kzh_70hPa)
        call mpas_pool_get_array(diag, 'kzh_100hPa', kzh_100hPa)
        call mpas_pool_get_array(diag, 'kzh_150hPa', kzh_150hPa)
        call mpas_pool_get_array(diag, 'kzh_200hPa', kzh_200hPa)
        call mpas_pool_get_array(diag, 'kzh_250hPa', kzh_250hPa)
        call mpas_pool_get_array(diag, 'kzh_300hPa', kzh_300hPa)
        call mpas_pool_get_array(diag, 'kzh_400hPa', kzh_400hPa)
        call mpas_pool_get_array(diag, 'kzh_500hPa', kzh_500hPa)
        call mpas_pool_get_array(diag, 'kzh_600hPa', kzh_600hPa)
        call mpas_pool_get_array(diag, 'kzh_700hPa', kzh_700hPa)
        call mpas_pool_get_array(diag, 'kzh_825hPa', kzh_825hPa)
        call mpas_pool_get_array(diag, 'kzh_850hPa', kzh_850hPa)
        call mpas_pool_get_array(diag, 'kzh_875hPa', kzh_875hPa)
        call mpas_pool_get_array(diag, 'kzh_900hPa', kzh_900hPa)
        call mpas_pool_get_array(diag, 'kzh_925hPa', kzh_925hPa)
        call mpas_pool_get_array(diag, 'kzh_950hPa', kzh_950hPa)
        call mpas_pool_get_array(diag, 'kzh_975hPa', kzh_975hPa)
        call mpas_pool_get_array(diag, 'kzh_1000hPa', kzh_1000hPa)

        call mpas_pool_get_array(diag, 'kzm_15hPa', kzm_15hPa)
        call mpas_pool_get_array(diag, 'kzm_20hPa', kzm_20hPa)
        call mpas_pool_get_array(diag, 'kzm_30hPa', kzm_30hPa)
        call mpas_pool_get_array(diag, 'kzm_50hPa', kzm_50hPa)
        call mpas_pool_get_array(diag, 'kzm_70hPa', kzm_70hPa)
        call mpas_pool_get_array(diag, 'kzm_100hPa', kzm_100hPa)
        call mpas_pool_get_array(diag, 'kzm_150hPa', kzm_150hPa)
        call mpas_pool_get_array(diag, 'kzm_200hPa', kzm_200hPa)
        call mpas_pool_get_array(diag, 'kzm_250hPa', kzm_250hPa)
        call mpas_pool_get_array(diag, 'kzm_300hPa', kzm_300hPa)
        call mpas_pool_get_array(diag, 'kzm_400hPa', kzm_400hPa)
        call mpas_pool_get_array(diag, 'kzm_500hPa', kzm_500hPa)
        call mpas_pool_get_array(diag, 'kzm_600hPa', kzm_600hPa)
        call mpas_pool_get_array(diag, 'kzm_700hPa', kzm_700hPa)
        call mpas_pool_get_array(diag, 'kzm_825hPa', kzm_825hPa)
        call mpas_pool_get_array(diag, 'kzm_850hPa', kzm_850hPa)
        call mpas_pool_get_array(diag, 'kzm_875hPa', kzm_875hPa)
        call mpas_pool_get_array(diag, 'kzm_900hPa', kzm_900hPa)
        call mpas_pool_get_array(diag, 'kzm_925hPa', kzm_925hPa)
        call mpas_pool_get_array(diag, 'kzm_950hPa', kzm_950hPa)
        call mpas_pool_get_array(diag, 'kzm_975hPa', kzm_975hPa)
        call mpas_pool_get_array(diag, 'kzm_1000hPa', kzm_1000hPa)

        call mpas_pool_get_array(diag, 'ni_15hPa', ni_15hPa)
        call mpas_pool_get_array(diag, 'ni_20hPa', ni_20hPa)
        call mpas_pool_get_array(diag, 'ni_30hPa', ni_30hPa)
        call mpas_pool_get_array(diag, 'ni_50hPa', ni_50hPa)
        call mpas_pool_get_array(diag, 'ni_70hPa', ni_70hPa)
        call mpas_pool_get_array(diag, 'ni_100hPa', ni_100hPa)
        call mpas_pool_get_array(diag, 'ni_150hPa', ni_150hPa)
        call mpas_pool_get_array(diag, 'ni_200hPa', ni_200hPa)
        call mpas_pool_get_array(diag, 'ni_250hPa', ni_250hPa)
        call mpas_pool_get_array(diag, 'ni_300hPa', ni_300hPa)
        call mpas_pool_get_array(diag, 'ni_400hPa', ni_400hPa)
        call mpas_pool_get_array(diag, 'ni_500hPa', ni_500hPa)
        call mpas_pool_get_array(diag, 'ni_600hPa', ni_600hPa)
        call mpas_pool_get_array(diag, 'ni_700hPa', ni_700hPa)
        call mpas_pool_get_array(diag, 'ni_825hPa', ni_825hPa)
        call mpas_pool_get_array(diag, 'ni_850hPa', ni_850hPa)
        call mpas_pool_get_array(diag, 'ni_875hPa', ni_875hPa)
        call mpas_pool_get_array(diag, 'ni_900hPa', ni_900hPa)
        call mpas_pool_get_array(diag, 'ni_925hPa', ni_925hPa)
        call mpas_pool_get_array(diag, 'ni_950hPa', ni_950hPa)
        call mpas_pool_get_array(diag, 'ni_975hPa', ni_975hPa)
        call mpas_pool_get_array(diag, 'ni_1000hPa', ni_1000hPa)

        call mpas_pool_get_array(diag, 'nr_15hPa', nr_15hPa)
        call mpas_pool_get_array(diag, 'nr_20hPa', nr_20hPa)
        call mpas_pool_get_array(diag, 'nr_30hPa', nr_30hPa)
        call mpas_pool_get_array(diag, 'nr_50hPa', nr_50hPa)
        call mpas_pool_get_array(diag, 'nr_70hPa', nr_70hPa)
        call mpas_pool_get_array(diag, 'nr_100hPa', nr_100hPa)
        call mpas_pool_get_array(diag, 'nr_150hPa', nr_150hPa)
        call mpas_pool_get_array(diag, 'nr_200hPa', nr_200hPa)
        call mpas_pool_get_array(diag, 'nr_250hPa', nr_250hPa)
        call mpas_pool_get_array(diag, 'nr_300hPa', nr_300hPa)
        call mpas_pool_get_array(diag, 'nr_400hPa', nr_400hPa)
        call mpas_pool_get_array(diag, 'nr_500hPa', nr_500hPa)
        call mpas_pool_get_array(diag, 'nr_600hPa', nr_600hPa)
        call mpas_pool_get_array(diag, 'nr_700hPa', nr_700hPa)
        call mpas_pool_get_array(diag, 'nr_825hPa', nr_825hPa)
        call mpas_pool_get_array(diag, 'nr_850hPa', nr_850hPa)
        call mpas_pool_get_array(diag, 'nr_875hPa', nr_875hPa)
        call mpas_pool_get_array(diag, 'nr_900hPa', nr_900hPa)
        call mpas_pool_get_array(diag, 'nr_925hPa', nr_925hPa)
        call mpas_pool_get_array(diag, 'nr_950hPa', nr_950hPa)
        call mpas_pool_get_array(diag, 'nr_975hPa', nr_975hPa)
        call mpas_pool_get_array(diag, 'nr_1000hPa', nr_1000hPa)

        call mpas_pool_get_array(diag, 'qr_15hPa', qr_15hPa)
        call mpas_pool_get_array(diag, 'qr_20hPa', qr_20hPa)
        call mpas_pool_get_array(diag, 'qr_30hPa', qr_30hPa)
        call mpas_pool_get_array(diag, 'qr_50hPa', qr_50hPa)
        call mpas_pool_get_array(diag, 'qr_70hPa', qr_70hPa)
        call mpas_pool_get_array(diag, 'qr_100hPa', qr_100hPa)
        call mpas_pool_get_array(diag, 'qr_150hPa', qr_150hPa)
        call mpas_pool_get_array(diag, 'qr_200hPa', qr_200hPa)
        call mpas_pool_get_array(diag, 'qr_250hPa', qr_250hPa)
        call mpas_pool_get_array(diag, 'qr_300hPa', qr_300hPa)
        call mpas_pool_get_array(diag, 'qr_400hPa', qr_400hPa)
        call mpas_pool_get_array(diag, 'qr_500hPa', qr_500hPa)
        call mpas_pool_get_array(diag, 'qr_600hPa', qr_600hPa)
        call mpas_pool_get_array(diag, 'qr_700hPa', qr_700hPa)
        call mpas_pool_get_array(diag, 'qr_825hPa', qr_825hPa)
        call mpas_pool_get_array(diag, 'qr_850hPa', qr_850hPa)
        call mpas_pool_get_array(diag, 'qr_875hPa', qr_875hPa)
        call mpas_pool_get_array(diag, 'qr_900hPa', qr_900hPa)
        call mpas_pool_get_array(diag, 'qr_925hPa', qr_925hPa)
        call mpas_pool_get_array(diag, 'qr_950hPa', qr_950hPa)
        call mpas_pool_get_array(diag, 'qr_975hPa', qr_975hPa)
        call mpas_pool_get_array(diag, 'qr_1000hPa', qr_1000hPa)

        call mpas_pool_get_array(diag, 'qs_15hPa', qs_15hPa)
        call mpas_pool_get_array(diag, 'qs_20hPa', qs_20hPa)
        call mpas_pool_get_array(diag, 'qs_30hPa', qs_30hPa)
        call mpas_pool_get_array(diag, 'qs_50hPa', qs_50hPa)
        call mpas_pool_get_array(diag, 'qs_70hPa', qs_70hPa)
        call mpas_pool_get_array(diag, 'qs_100hPa', qs_100hPa)
        call mpas_pool_get_array(diag, 'qs_150hPa', qs_150hPa)
        call mpas_pool_get_array(diag, 'qs_200hPa', qs_200hPa)
        call mpas_pool_get_array(diag, 'qs_250hPa', qs_250hPa)
        call mpas_pool_get_array(diag, 'qs_300hPa', qs_300hPa)
        call mpas_pool_get_array(diag, 'qs_400hPa', qs_400hPa)
        call mpas_pool_get_array(diag, 'qs_500hPa', qs_500hPa)
        call mpas_pool_get_array(diag, 'qs_600hPa', qs_600hPa)
        call mpas_pool_get_array(diag, 'qs_700hPa', qs_700hPa)
        call mpas_pool_get_array(diag, 'qs_825hPa', qs_825hPa)
        call mpas_pool_get_array(diag, 'qs_850hPa', qs_850hPa)
        call mpas_pool_get_array(diag, 'qs_875hPa', qs_875hPa)
        call mpas_pool_get_array(diag, 'qs_900hPa', qs_900hPa)
        call mpas_pool_get_array(diag, 'qs_925hPa', qs_925hPa)
        call mpas_pool_get_array(diag, 'qs_950hPa', qs_950hPa)
        call mpas_pool_get_array(diag, 'qs_975hPa', qs_975hPa)
        call mpas_pool_get_array(diag, 'qs_1000hPa', qs_1000hPa)

        call mpas_pool_get_array(diag, 'qv_15hPa', qv_15hPa)
        call mpas_pool_get_array(diag, 'qv_20hPa', qv_20hPa)
        call mpas_pool_get_array(diag, 'qv_30hPa', qv_30hPa)
        call mpas_pool_get_array(diag, 'qv_50hPa', qv_50hPa)
        call mpas_pool_get_array(diag, 'qv_70hPa', qv_70hPa)
        call mpas_pool_get_array(diag, 'qv_100hPa', qv_100hPa)
        call mpas_pool_get_array(diag, 'qv_150hPa', qv_150hPa)
        call mpas_pool_get_array(diag, 'qv_200hPa', qv_200hPa)
        call mpas_pool_get_array(diag, 'qv_250hPa', qv_250hPa)
        call mpas_pool_get_array(diag, 'qv_300hPa', qv_300hPa)
        call mpas_pool_get_array(diag, 'qv_400hPa', qv_400hPa)
        call mpas_pool_get_array(diag, 'qv_500hPa', qv_500hPa)
        call mpas_pool_get_array(diag, 'qv_600hPa', qv_600hPa)
        call mpas_pool_get_array(diag, 'qv_700hPa', qv_700hPa)
        call mpas_pool_get_array(diag, 'qv_825hPa', qv_825hPa)
        call mpas_pool_get_array(diag, 'qv_850hPa', qv_850hPa)
        call mpas_pool_get_array(diag, 'qv_875hPa', qv_875hPa)
        call mpas_pool_get_array(diag, 'qv_900hPa', qv_900hPa)
        call mpas_pool_get_array(diag, 'qv_925hPa', qv_925hPa)
        call mpas_pool_get_array(diag, 'qv_950hPa', qv_950hPa)
        call mpas_pool_get_array(diag, 'qv_975hPa', qv_975hPa)
        call mpas_pool_get_array(diag, 'qv_1000hPa', qv_1000hPa)

        call mpas_pool_get_array(diag, 'qc_15hPa', qc_15hPa)
        call mpas_pool_get_array(diag, 'qc_20hPa', qc_20hPa)
        call mpas_pool_get_array(diag, 'qc_30hPa', qc_30hPa)
        call mpas_pool_get_array(diag, 'qc_50hPa', qc_50hPa)
        call mpas_pool_get_array(diag, 'qc_70hPa', qc_70hPa)
        call mpas_pool_get_array(diag, 'qc_100hPa', qc_100hPa)
        call mpas_pool_get_array(diag, 'qc_150hPa', qc_150hPa)
        call mpas_pool_get_array(diag, 'qc_200hPa', qc_200hPa)
        call mpas_pool_get_array(diag, 'qc_250hPa', qc_250hPa)
        call mpas_pool_get_array(diag, 'qc_300hPa', qc_300hPa)
        call mpas_pool_get_array(diag, 'qc_400hPa', qc_400hPa)
        call mpas_pool_get_array(diag, 'qc_500hPa', qc_500hPa)
        call mpas_pool_get_array(diag, 'qc_600hPa', qc_600hPa)
        call mpas_pool_get_array(diag, 'qc_700hPa', qc_700hPa)
        call mpas_pool_get_array(diag, 'qc_825hPa', qc_825hPa)
        call mpas_pool_get_array(diag, 'qc_850hPa', qc_850hPa)
        call mpas_pool_get_array(diag, 'qc_875hPa', qc_875hPa)
        call mpas_pool_get_array(diag, 'qc_900hPa', qc_900hPa)
        call mpas_pool_get_array(diag, 'qc_925hPa', qc_925hPa)
        call mpas_pool_get_array(diag, 'qc_950hPa', qc_950hPa)
        call mpas_pool_get_array(diag, 'qc_975hPa', qc_975hPa)
        call mpas_pool_get_array(diag, 'qc_1000hPa', qc_1000hPa)

        call mpas_pool_get_array(diag, 'qg_15hPa', qg_15hPa)
        call mpas_pool_get_array(diag, 'qg_20hPa', qg_20hPa)
        call mpas_pool_get_array(diag, 'qg_30hPa', qg_30hPa)
        call mpas_pool_get_array(diag, 'qg_50hPa', qg_50hPa)
        call mpas_pool_get_array(diag, 'qg_70hPa', qg_70hPa)
        call mpas_pool_get_array(diag, 'qg_100hPa', qg_100hPa)
        call mpas_pool_get_array(diag, 'qg_150hPa', qg_150hPa)
        call mpas_pool_get_array(diag, 'qg_200hPa', qg_200hPa)
        call mpas_pool_get_array(diag, 'qg_250hPa', qg_250hPa)
        call mpas_pool_get_array(diag, 'qg_300hPa', qg_300hPa)
        call mpas_pool_get_array(diag, 'qg_400hPa', qg_400hPa)
        call mpas_pool_get_array(diag, 'qg_500hPa', qg_500hPa)
        call mpas_pool_get_array(diag, 'qg_600hPa', qg_600hPa)
        call mpas_pool_get_array(diag, 'qg_700hPa', qg_700hPa)
        call mpas_pool_get_array(diag, 'qg_825hPa', qg_825hPa)
        call mpas_pool_get_array(diag, 'qg_850hPa', qg_850hPa)
        call mpas_pool_get_array(diag, 'qg_875hPa', qg_875hPa)
        call mpas_pool_get_array(diag, 'qg_900hPa', qg_900hPa)
        call mpas_pool_get_array(diag, 'qg_925hPa', qg_925hPa)
        call mpas_pool_get_array(diag, 'qg_950hPa', qg_950hPa)
        call mpas_pool_get_array(diag, 'qg_975hPa', qg_975hPa)
        call mpas_pool_get_array(diag, 'qg_1000hPa', qg_1000hPa)

        call mpas_pool_get_array(diag, 'qi_15hPa', qi_15hPa)
        call mpas_pool_get_array(diag, 'qi_20hPa', qi_20hPa)
        call mpas_pool_get_array(diag, 'qi_30hPa', qi_30hPa)
        call mpas_pool_get_array(diag, 'qi_50hPa', qi_50hPa)
        call mpas_pool_get_array(diag, 'qi_70hPa', qi_70hPa)
        call mpas_pool_get_array(diag, 'qi_100hPa', qi_100hPa)
        call mpas_pool_get_array(diag, 'qi_150hPa', qi_150hPa)
        call mpas_pool_get_array(diag, 'qi_200hPa', qi_200hPa)
        call mpas_pool_get_array(diag, 'qi_250hPa', qi_250hPa)
        call mpas_pool_get_array(diag, 'qi_300hPa', qi_300hPa)
        call mpas_pool_get_array(diag, 'qi_400hPa', qi_400hPa)
        call mpas_pool_get_array(diag, 'qi_500hPa', qi_500hPa)
        call mpas_pool_get_array(diag, 'qi_600hPa', qi_600hPa)
        call mpas_pool_get_array(diag, 'qi_700hPa', qi_700hPa)
        call mpas_pool_get_array(diag, 'qi_825hPa', qi_825hPa)
        call mpas_pool_get_array(diag, 'qi_850hPa', qi_850hPa)
        call mpas_pool_get_array(diag, 'qi_875hPa', qi_875hPa)
        call mpas_pool_get_array(diag, 'qi_900hPa', qi_900hPa)
        call mpas_pool_get_array(diag, 'qi_925hPa', qi_925hPa)
        call mpas_pool_get_array(diag, 'qi_950hPa', qi_950hPa)
        call mpas_pool_get_array(diag, 'qi_975hPa', qi_975hPa)
        call mpas_pool_get_array(diag, 'qi_1000hPa', qi_1000hPa)

        call mpas_pool_get_array(diag, 'vorticity_50hPa', vorticity_50hPa)
        call mpas_pool_get_array(diag, 'vorticity_100hPa', vorticity_100hPa)
        call mpas_pool_get_array(diag, 'vorticity_200hPa', vorticity_200hPa)
        call mpas_pool_get_array(diag, 'vorticity_250hPa', vorticity_250hPa)
        call mpas_pool_get_array(diag, 'vorticity_500hPa', vorticity_500hPa)
        call mpas_pool_get_array(diag, 'vorticity_700hPa', vorticity_700hPa)
        call mpas_pool_get_array(diag, 'vorticity_850hPa', vorticity_850hPa)
        call mpas_pool_get_array(diag, 'vorticity_925hPa', vorticity_925hPa)
       
        call mpas_pool_get_array(diag, 'mslp', mslp)
       
        if(.not.allocated(pressure)    ) allocate(pressure(nVertLevels,nCells)      )
        if(.not.allocated(pressureCp1) ) allocate(pressureCp1(nVertLevels,nCells+1) )
        if(.not.allocated(pressure2)   ) allocate(pressure2(nVertLevelsP1,nCells)   )
        if(.not.allocated(pressure_v)  ) allocate(pressure_v(nVertLevels,nVertices) )
        if(.not.allocated(temperature) ) allocate(temperature(nVertLevels,nCells)   )
        if(.not.allocated(spechum)  ) allocate(spechum(nVertLevels,nCells)    ) !LFR
        if(.not.allocated(geoph)    ) allocate(geoph(nVertLevelsP1,nCells)    ) !EGK
        if(.not.allocated(dewpoint) ) allocate(dewpoint(nVertLevels,nCells)   )
       
        if (need_t_isobaric) then
            t_iso_levels(1) = 30000.0
            t_iso_levels(2) = 35000.0
            t_iso_levels(3) = 40000.0
            t_iso_levels(4) = 45000.0
            t_iso_levels(5) = 50000.0
        end if
       
        if (need_z_isobaric) then
            z_iso_levels(1)  = 30000.0
            z_iso_levels(2)  = 35000.0
            z_iso_levels(3)  = 40000.0
            z_iso_levels(4)  = 45000.0
            z_iso_levels(5)  = 50000.0
            z_iso_levels(6)  = 55000.0
            z_iso_levels(7)  = 60000.0
            z_iso_levels(8)  = 65000.0
            z_iso_levels(9)  = 70000.0
            z_iso_levels(10) = 75000.0
            z_iso_levels(11) = 80000.0
            z_iso_levels(12) = 85000.0
            z_iso_levels(13) = 90000.0
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
       
        if (NEED_TEMP .or. NEED_RELHUM .or. NEED_SPECHUM .or. NEED_DEWPOINT .or. need_mslp) then
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
        nIntP = 22
        if(.not.allocated(field_interp)) allocate(field_interp(nCells,nIntP) )
        if(.not.allocated(press_interp)) allocate(press_interp(nCells,nIntP) )
        do iCell = 1, nCells
           press_interp(iCell,1)  = 15.0_RKIND
           press_interp(iCell,2)  = 20.0_RKIND
           press_interp(iCell,3)  = 30.0_RKIND
           press_interp(iCell,4)  = 50.0_RKIND
           press_interp(iCell,5)  = 70.0_RKIND
           press_interp(iCell,6)  = 100.0_RKIND
           press_interp(iCell,7)  = 150.0_RKIND
           press_interp(iCell,8)  = 200.0_RKIND
           press_interp(iCell,9)  = 250.0_RKIND
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
       
        if(.not.allocated(press_in)) allocate(press_in(nCells,nVertLevels))
        do iCell = 1, nCells
        do k = 1, nVertLevels
           kk = nVertLevels+1-k
           press_in(iCell,kk) = pressure(k,iCell)
        enddo
        enddo
       
        if(.not.allocated(field_in)) allocate(field_in(nCells,nVertLevels))

        if (NEED_CLDFRAC) then
           !... temperature:
            do iCell = 1, nCells
               do k = 1, nVertLevels
                  kk = nVertLevels+1-k
                  field_in(iCell,kk) = cldfrac(k,iCell)
                  !write(88,*) "LFR-DBG: k,iCell,cldFrac ",k,iCell,cldfrac(k,iCell)
               enddo
            enddo
            call interp_tofixed_pressure(nCells,nVertLevels,nIntP,press_in,field_in,press_interp,field_interp)

            cldfrac_15hPa(1:nCells)   = field_interp(1:nCells,1)
            cldfrac_20hPa(1:nCells)   = field_interp(1:nCells,2)
            cldfrac_30hPa(1:nCells)   = field_interp(1:nCells,3)
            cldfrac_50hPa(1:nCells)   = field_interp(1:nCells,4)
            cldfrac_70hPa(1:nCells)   = field_interp(1:nCells,5)
            cldfrac_100hPa(1:nCells)  = field_interp(1:nCells,6)
            cldfrac_150hPa(1:nCells)  = field_interp(1:nCells,7)
            cldfrac_200hPa(1:nCells)  = field_interp(1:nCells,8)
            cldfrac_250hPa(1:nCells)  = field_interp(1:nCells,9)
            cldfrac_300hPa(1:nCells)  = field_interp(1:nCells,10)
            cldfrac_400hPa(1:nCells)  = field_interp(1:nCells,11)
            cldfrac_500hPa(1:nCells)  = field_interp(1:nCells,12)
            cldfrac_600hPa(1:nCells)  = field_interp(1:nCells,13)
            cldfrac_700hPa(1:nCells)  = field_interp(1:nCells,14)
            cldfrac_825hPa(1:nCells)  = field_interp(1:nCells,15)
            cldfrac_850hPa(1:nCells)  = field_interp(1:nCells,16)
            cldfrac_875hPa(1:nCells)  = field_interp(1:nCells,17)
            cldfrac_900hPa(1:nCells)  = field_interp(1:nCells,18)
            cldfrac_925hPa(1:nCells)  = field_interp(1:nCells,19)
            cldfrac_950hPa(1:nCells)  = field_interp(1:nCells,20)
            cldfrac_975hPa(1:nCells)  = field_interp(1:nCells,21)
            cldfrac_1000hPa(1:nCells) = field_interp(1:nCells,22)
            ! call mpas_log_write('--- end interpolate temperature:')
        end if


        if (NEED_TEMP) then
           !... temperature:
            do iCell = 1, nCells
            do k = 1, nVertLevels
               kk = nVertLevels+1-k
               field_in(iCell,kk) = temperature(k,iCell)
            enddo
            enddo
            call interp_tofixed_pressure(nCells,nVertLevels,nIntP,press_in,field_in,press_interp,field_interp)

            temperature_15hPa(1:nCells)   = field_interp(1:nCells,1)
            temperature_20hPa(1:nCells)   = field_interp(1:nCells,2)
            temperature_30hPa(1:nCells)   = field_interp(1:nCells,3)
            temperature_50hPa(1:nCells)   = field_interp(1:nCells,4)
            temperature_70hPa(1:nCells)   = field_interp(1:nCells,5)
            temperature_100hPa(1:nCells)  = field_interp(1:nCells,6)
            temperature_150hPa(1:nCells)  = field_interp(1:nCells,7)
            temperature_200hPa(1:nCells)  = field_interp(1:nCells,8)
            temperature_250hPa(1:nCells)  = field_interp(1:nCells,9)
            temperature_300hPa(1:nCells)  = field_interp(1:nCells,10)
            temperature_400hPa(1:nCells)  = field_interp(1:nCells,11)
            temperature_500hPa(1:nCells)  = field_interp(1:nCells,12)
            temperature_600hPa(1:nCells)  = field_interp(1:nCells,13)
            temperature_700hPa(1:nCells)  = field_interp(1:nCells,14)
            temperature_825hPa(1:nCells)  = field_interp(1:nCells,15)
            temperature_850hPa(1:nCells)  = field_interp(1:nCells,16)
            temperature_875hPa(1:nCells)  = field_interp(1:nCells,17)
            temperature_900hPa(1:nCells)  = field_interp(1:nCells,18)
            temperature_925hPa(1:nCells)  = field_interp(1:nCells,19)
            temperature_950hPa(1:nCells)  = field_interp(1:nCells,20)
            temperature_975hPa(1:nCells)  = field_interp(1:nCells,21)
            temperature_1000hPa(1:nCells) = field_interp(1:nCells,22)
            ! call mpas_log_write('--- end interpolate temperature:')
        end if

!        print *, "mp_thompson3=", need_mp_thompson
!        print *, "cu_grell_freitas3=", need_cu_grell_freitas
!        print *, "bl_mynn3=", need_bl_mynn

        if (need_mp_thompson .eqv. .true. .or. need_bl_mynn .eqv. .true.) then
            if (NEED_NI) then
               !... ni:
                do iCell = 1, nCells
                do k = 1, nVertLevels
                   kk = nVertLevels+1-k
                   field_in(iCell,kk) = scalars(index_ni,k,iCell)
                enddo
                enddo
                call interp_tofixed_pressure(nCells,nVertLevels,nIntP,press_in,field_in,press_interp,field_interp)
    
                ni_15hPa(1:nCells) = field_interp(1:nCells,1)
                ni_20hPa(1:nCells) = field_interp(1:nCells,2)
                ni_30hPa(1:nCells) = field_interp(1:nCells,3)
                ni_50hPa(1:nCells) = field_interp(1:nCells,4)
                ni_70hPa(1:nCells) = field_interp(1:nCells,5)
                ni_100hPa(1:nCells) = field_interp(1:nCells,6)
                ni_150hPa(1:nCells) = field_interp(1:nCells,7)
                ni_200hPa(1:nCells) = field_interp(1:nCells,8)
                ni_250hPa(1:nCells) = field_interp(1:nCells,9)
                ni_300hPa(1:nCells) = field_interp(1:nCells,10)
                ni_400hPa(1:nCells) = field_interp(1:nCells,11)
                ni_500hPa(1:nCells) = field_interp(1:nCells,12)
                ni_600hPa(1:nCells) = field_interp(1:nCells,13)
                ni_700hPa(1:nCells) = field_interp(1:nCells,14)
                ni_825hPa(1:nCells) = field_interp(1:nCells,15)
                ni_850hPa(1:nCells) = field_interp(1:nCells,16)
                ni_875hPa(1:nCells) = field_interp(1:nCells,17)
                ni_900hPa(1:nCells) = field_interp(1:nCells,18)
                ni_925hPa(1:nCells) = field_interp(1:nCells,19)
                ni_950hPa(1:nCells) = field_interp(1:nCells,20)
                ni_975hPa(1:nCells) = field_interp(1:nCells,21)
                ni_1000hPa(1:nCells) = field_interp(1:nCells,22)
                ! call mpas_log_write('--- end interpolate ni:')
            end if
        end if

        if (need_mp_thompson .eqv. .true.) then
            if (NEED_NR) then
               !... nr:
                do iCell = 1, nCells
                do k = 1, nVertLevels
                   kk = nVertLevels+1-k
                   field_in(iCell,kk) = scalars(index_nr,k,iCell)
                enddo
                enddo
                call interp_tofixed_pressure(nCells,nVertLevels,nIntP,press_in,field_in,press_interp,field_interp)
    
                nr_15hPa(1:nCells) = field_interp(1:nCells,1)
                nr_20hPa(1:nCells) = field_interp(1:nCells,2)
                nr_30hPa(1:nCells) = field_interp(1:nCells,3)
                nr_50hPa(1:nCells) = field_interp(1:nCells,4)
                nr_70hPa(1:nCells) = field_interp(1:nCells,5)
                nr_100hPa(1:nCells) = field_interp(1:nCells,6)
                nr_150hPa(1:nCells) = field_interp(1:nCells,7)
                nr_200hPa(1:nCells) = field_interp(1:nCells,8)
                nr_250hPa(1:nCells) = field_interp(1:nCells,9)
                nr_300hPa(1:nCells) = field_interp(1:nCells,10)
                nr_400hPa(1:nCells) = field_interp(1:nCells,11)
                nr_500hPa(1:nCells) = field_interp(1:nCells,12)
                nr_600hPa(1:nCells) = field_interp(1:nCells,13)
                nr_700hPa(1:nCells) = field_interp(1:nCells,14)
                nr_825hPa(1:nCells) = field_interp(1:nCells,15)
                nr_850hPa(1:nCells) = field_interp(1:nCells,16)
                nr_875hPa(1:nCells) = field_interp(1:nCells,17)
                nr_900hPa(1:nCells) = field_interp(1:nCells,18)
                nr_925hPa(1:nCells) = field_interp(1:nCells,19)
                nr_950hPa(1:nCells) = field_interp(1:nCells,20)
                nr_975hPa(1:nCells) = field_interp(1:nCells,21)
                nr_1000hPa(1:nCells) = field_interp(1:nCells,22)
                ! call mpas_log_write('--- end interpolate nr:')
            end if
        end if

        if (NEED_QR) then
           !... qr:
            do iCell = 1, nCells
            do k = 1, nVertLevels
               kk = nVertLevels+1-k
               field_in(iCell,kk) = scalars(index_qr,k,iCell)
            enddo
            enddo
            call interp_tofixed_pressure(nCells,nVertLevels,nIntP,press_in,field_in,press_interp,field_interp)

            qr_15hPa(1:nCells) = field_interp(1:nCells,1)
            qr_20hPa(1:nCells) = field_interp(1:nCells,2)
            qr_30hPa(1:nCells) = field_interp(1:nCells,3)
            qr_50hPa(1:nCells) = field_interp(1:nCells,4)
            qr_70hPa(1:nCells) = field_interp(1:nCells,5)
            qr_100hPa(1:nCells) = field_interp(1:nCells,6)
            qr_150hPa(1:nCells) = field_interp(1:nCells,7)
            qr_200hPa(1:nCells) = field_interp(1:nCells,8)
            qr_250hPa(1:nCells) = field_interp(1:nCells,9)
            qr_300hPa(1:nCells) = field_interp(1:nCells,10)
            qr_400hPa(1:nCells) = field_interp(1:nCells,11)
            qr_500hPa(1:nCells) = field_interp(1:nCells,12)
            qr_600hPa(1:nCells) = field_interp(1:nCells,13)
            qr_700hPa(1:nCells) = field_interp(1:nCells,14)
            qr_825hPa(1:nCells) = field_interp(1:nCells,15)
            qr_850hPa(1:nCells) = field_interp(1:nCells,16)
            qr_875hPa(1:nCells) = field_interp(1:nCells,17)
            qr_900hPa(1:nCells) = field_interp(1:nCells,18)
            qr_925hPa(1:nCells) = field_interp(1:nCells,19)
            qr_950hPa(1:nCells) = field_interp(1:nCells,20)
            qr_975hPa(1:nCells) = field_interp(1:nCells,21)
            qr_1000hPa(1:nCells) = field_interp(1:nCells,22)
            ! call mpas_log_write('--- end interpolate qr:')
        end if

        if (NEED_QS) then
           !... qs:
            do iCell = 1, nCells
            do k = 1, nVertLevels
               kk = nVertLevels+1-k
               field_in(iCell,kk) = scalars(index_qs,k,iCell)
            enddo
            enddo
            call interp_tofixed_pressure(nCells,nVertLevels,nIntP,press_in,field_in,press_interp,field_interp)

            qs_15hPa(1:nCells) = field_interp(1:nCells,1)
            qs_20hPa(1:nCells) = field_interp(1:nCells,2)
            qs_30hPa(1:nCells) = field_interp(1:nCells,3)
            qs_50hPa(1:nCells) = field_interp(1:nCells,4)
            qs_70hPa(1:nCells) = field_interp(1:nCells,5)
            qs_100hPa(1:nCells) = field_interp(1:nCells,6)
            qs_150hPa(1:nCells) = field_interp(1:nCells,7)
            qs_200hPa(1:nCells) = field_interp(1:nCells,8)
            qs_250hPa(1:nCells) = field_interp(1:nCells,9)
            qs_300hPa(1:nCells) = field_interp(1:nCells,10)
            qs_400hPa(1:nCells) = field_interp(1:nCells,11)
            qs_500hPa(1:nCells) = field_interp(1:nCells,12)
            qs_600hPa(1:nCells) = field_interp(1:nCells,13)
            qs_700hPa(1:nCells) = field_interp(1:nCells,14)
            qs_825hPa(1:nCells) = field_interp(1:nCells,15)
            qs_850hPa(1:nCells) = field_interp(1:nCells,16)
            qs_875hPa(1:nCells) = field_interp(1:nCells,17)
            qs_900hPa(1:nCells) = field_interp(1:nCells,18)
            qs_925hPa(1:nCells) = field_interp(1:nCells,19)
            qs_950hPa(1:nCells) = field_interp(1:nCells,20)
            qs_975hPa(1:nCells) = field_interp(1:nCells,21)
            qs_1000hPa(1:nCells) = field_interp(1:nCells,22)
            ! call mpas_log_write('--- end interpolate qs:')
        end if

        if (NEED_QV) then
           !... qv:
            do iCell = 1, nCells
            do k = 1, nVertLevels
               kk = nVertLevels+1-k
               field_in(iCell,kk) = scalars(index_qv,k,iCell)
            enddo
            enddo
            call interp_tofixed_pressure(nCells,nVertLevels,nIntP,press_in,field_in,press_interp,field_interp)

            qv_15hPa(1:nCells) = field_interp(1:nCells,1)
            qv_20hPa(1:nCells) = field_interp(1:nCells,2)
            qv_30hPa(1:nCells) = field_interp(1:nCells,3)
            qv_50hPa(1:nCells) = field_interp(1:nCells,4)
            qv_70hPa(1:nCells) = field_interp(1:nCells,5)
            qv_100hPa(1:nCells) = field_interp(1:nCells,6)
            qv_150hPa(1:nCells) = field_interp(1:nCells,7)
            qv_200hPa(1:nCells) = field_interp(1:nCells,8)
            qv_250hPa(1:nCells) = field_interp(1:nCells,9)
            qv_300hPa(1:nCells) = field_interp(1:nCells,10)
            qv_400hPa(1:nCells) = field_interp(1:nCells,11)
            qv_500hPa(1:nCells) = field_interp(1:nCells,12)
            qv_600hPa(1:nCells) = field_interp(1:nCells,13)
            qv_700hPa(1:nCells) = field_interp(1:nCells,14)
            qv_825hPa(1:nCells) = field_interp(1:nCells,15)
            qv_850hPa(1:nCells) = field_interp(1:nCells,16)
            qv_875hPa(1:nCells) = field_interp(1:nCells,17)
            qv_900hPa(1:nCells) = field_interp(1:nCells,18)
            qv_925hPa(1:nCells) = field_interp(1:nCells,19)
            qv_950hPa(1:nCells) = field_interp(1:nCells,20)
            qv_975hPa(1:nCells) = field_interp(1:nCells,21)
            qv_1000hPa(1:nCells) = field_interp(1:nCells,22)
            ! call mpas_log_write('--- end interpolate qv:')
        end if

        if (NEED_QC) then
           !... qc:
            do iCell = 1, nCells
            do k = 1, nVertLevels
               kk = nVertLevels+1-k
               field_in(iCell,kk) = scalars(index_qc,k,iCell)
!               field_in(iCell,kk) = 1.
            enddo
            enddo
            call interp_tofixed_pressure(nCells,nVertLevels,nIntP,press_in,field_in,press_interp,field_interp)

            qc_15hPa(1:nCells) = field_interp(1:nCells,1)
            qc_20hPa(1:nCells) = field_interp(1:nCells,2)
            qc_30hPa(1:nCells) = field_interp(1:nCells,3)
            qc_50hPa(1:nCells) = field_interp(1:nCells,4)
            qc_70hPa(1:nCells) = field_interp(1:nCells,5)
            qc_100hPa(1:nCells) = field_interp(1:nCells,6)
            qc_150hPa(1:nCells) = field_interp(1:nCells,7)
            qc_200hPa(1:nCells) = field_interp(1:nCells,8)
            qc_250hPa(1:nCells) = field_interp(1:nCells,9)
            qc_300hPa(1:nCells) = field_interp(1:nCells,10)
            qc_400hPa(1:nCells) = field_interp(1:nCells,11)
            qc_500hPa(1:nCells) = field_interp(1:nCells,12)
            qc_600hPa(1:nCells) = field_interp(1:nCells,13)
            qc_700hPa(1:nCells) = field_interp(1:nCells,14)
            qc_825hPa(1:nCells) = field_interp(1:nCells,15)
            qc_850hPa(1:nCells) = field_interp(1:nCells,16)
            qc_875hPa(1:nCells) = field_interp(1:nCells,17)
            qc_900hPa(1:nCells) = field_interp(1:nCells,18)
            qc_925hPa(1:nCells) = field_interp(1:nCells,19)
            qc_950hPa(1:nCells) = field_interp(1:nCells,20)
            qc_975hPa(1:nCells) = field_interp(1:nCells,21)
            qc_1000hPa(1:nCells) = field_interp(1:nCells,22)
            ! call mpas_log_write('--- end interpolate qc:')
        end if

        if (NEED_QG) then
           !... qg:
            do iCell = 1, nCells
            do k = 1, nVertLevels
               kk = nVertLevels+1-k
               field_in(iCell,kk) = scalars(index_qg,k,iCell)
            enddo
            enddo
            call interp_tofixed_pressure(nCells,nVertLevels,nIntP,press_in,field_in,press_interp,field_interp)

            qg_15hPa(1:nCells) = field_interp(1:nCells,1)
            qg_20hPa(1:nCells) = field_interp(1:nCells,2)
            qg_30hPa(1:nCells) = field_interp(1:nCells,3)
            qg_50hPa(1:nCells) = field_interp(1:nCells,4)
            qg_70hPa(1:nCells) = field_interp(1:nCells,5)
            qg_100hPa(1:nCells) = field_interp(1:nCells,6)
            qg_150hPa(1:nCells) = field_interp(1:nCells,7)
            qg_200hPa(1:nCells) = field_interp(1:nCells,8)
            qg_250hPa(1:nCells) = field_interp(1:nCells,9)
            qg_300hPa(1:nCells) = field_interp(1:nCells,10)
            qg_400hPa(1:nCells) = field_interp(1:nCells,11)
            qg_500hPa(1:nCells) = field_interp(1:nCells,12)
            qg_600hPa(1:nCells) = field_interp(1:nCells,13)
            qg_700hPa(1:nCells) = field_interp(1:nCells,14)
            qg_825hPa(1:nCells) = field_interp(1:nCells,15)
            qg_850hPa(1:nCells) = field_interp(1:nCells,16)
            qg_875hPa(1:nCells) = field_interp(1:nCells,17)
            qg_900hPa(1:nCells) = field_interp(1:nCells,18)
            qg_925hPa(1:nCells) = field_interp(1:nCells,19)
            qg_950hPa(1:nCells) = field_interp(1:nCells,20)
            qg_975hPa(1:nCells) = field_interp(1:nCells,21)
            qg_1000hPa(1:nCells) = field_interp(1:nCells,22)
            ! call mpas_log_write('--- end interpolate qg:')
        end if

        if (NEED_QI) then
           !... qi:
            do iCell = 1, nCells
            do k = 1, nVertLevels
               kk = nVertLevels+1-k
               field_in(iCell,kk) = scalars(index_qi,k,iCell)
            enddo
            enddo
            call interp_tofixed_pressure(nCells,nVertLevels,nIntP,press_in,field_in,press_interp,field_interp)

            qi_15hPa(1:nCells) = field_interp(1:nCells,1)
            qi_20hPa(1:nCells) = field_interp(1:nCells,2)
            qi_30hPa(1:nCells) = field_interp(1:nCells,3)
            qi_50hPa(1:nCells) = field_interp(1:nCells,4)
            qi_70hPa(1:nCells) = field_interp(1:nCells,5)
            qi_100hPa(1:nCells) = field_interp(1:nCells,6)
            qi_150hPa(1:nCells) = field_interp(1:nCells,7)
            qi_200hPa(1:nCells) = field_interp(1:nCells,8)
            qi_250hPa(1:nCells) = field_interp(1:nCells,9)
            qi_300hPa(1:nCells) = field_interp(1:nCells,10)
            qi_400hPa(1:nCells) = field_interp(1:nCells,11)
            qi_500hPa(1:nCells) = field_interp(1:nCells,12)
            qi_600hPa(1:nCells) = field_interp(1:nCells,13)
            qi_700hPa(1:nCells) = field_interp(1:nCells,14)
            qi_825hPa(1:nCells) = field_interp(1:nCells,15)
            qi_850hPa(1:nCells) = field_interp(1:nCells,16)
            qi_875hPa(1:nCells) = field_interp(1:nCells,17)
            qi_900hPa(1:nCells) = field_interp(1:nCells,18)
            qi_925hPa(1:nCells) = field_interp(1:nCells,19)
            qi_950hPa(1:nCells) = field_interp(1:nCells,20)
            qi_975hPa(1:nCells) = field_interp(1:nCells,21)
            qi_1000hPa(1:nCells) = field_interp(1:nCells,22)
            ! call mpas_log_write('--- end interpolate qi:')
        end if

        if (NEED_RQVBLTEN) then
           !... rqvblten:
            do iCell = 1, nCells
            do k = 1, nVertLevels
               kk = nVertLevels+1-k
               field_in(iCell,kk) = rqvblten(k,iCell)
            enddo
            enddo
            call interp_tofixed_pressure(nCells,nVertLevels,nIntP,press_in,field_in,press_interp,field_interp)

            rqvblten_15hPa(1:nCells) = field_interp(1:nCells,1)
            rqvblten_20hPa(1:nCells) = field_interp(1:nCells,2)
            rqvblten_30hPa(1:nCells) = field_interp(1:nCells,3)
            rqvblten_50hPa(1:nCells) = field_interp(1:nCells,4)
            rqvblten_70hPa(1:nCells) = field_interp(1:nCells,5)
            rqvblten_100hPa(1:nCells) = field_interp(1:nCells,6)
            rqvblten_150hPa(1:nCells) = field_interp(1:nCells,7)
            rqvblten_200hPa(1:nCells) = field_interp(1:nCells,8)
            rqvblten_250hPa(1:nCells) = field_interp(1:nCells,9)
            rqvblten_300hPa(1:nCells) = field_interp(1:nCells,10)
            rqvblten_400hPa(1:nCells) = field_interp(1:nCells,11)
            rqvblten_500hPa(1:nCells) = field_interp(1:nCells,12)
            rqvblten_600hPa(1:nCells) = field_interp(1:nCells,13)
            rqvblten_700hPa(1:nCells) = field_interp(1:nCells,14)
            rqvblten_825hPa(1:nCells) = field_interp(1:nCells,15)
            rqvblten_850hPa(1:nCells) = field_interp(1:nCells,16)
            rqvblten_875hPa(1:nCells) = field_interp(1:nCells,17)
            rqvblten_900hPa(1:nCells) = field_interp(1:nCells,18)
            rqvblten_925hPa(1:nCells) = field_interp(1:nCells,19)
            rqvblten_950hPa(1:nCells) = field_interp(1:nCells,20)
            rqvblten_975hPa(1:nCells) = field_interp(1:nCells,21)
            rqvblten_1000hPa(1:nCells) = field_interp(1:nCells,22)
            ! call mpas_log_write('--- end interpolate rqvblten:')
        end if

        if (NEED_RQVCUTEN) then
           !... rqvcuten:
            do iCell = 1, nCells
            do k = 1, nVertLevels
               kk = nVertLevels+1-k
               field_in(iCell,kk) = rqvcuten(k,iCell)
            enddo
            enddo
            call interp_tofixed_pressure(nCells,nVertLevels,nIntP,press_in,field_in,press_interp,field_interp)

            rqvcuten_15hPa(1:nCells) = field_interp(1:nCells,1)
            rqvcuten_20hPa(1:nCells) = field_interp(1:nCells,2)
            rqvcuten_30hPa(1:nCells) = field_interp(1:nCells,3)
            rqvcuten_50hPa(1:nCells) = field_interp(1:nCells,4)
            rqvcuten_70hPa(1:nCells) = field_interp(1:nCells,5)
            rqvcuten_100hPa(1:nCells) = field_interp(1:nCells,6)
            rqvcuten_150hPa(1:nCells) = field_interp(1:nCells,7)
            rqvcuten_200hPa(1:nCells) = field_interp(1:nCells,8)
            rqvcuten_250hPa(1:nCells) = field_interp(1:nCells,9)
            rqvcuten_300hPa(1:nCells) = field_interp(1:nCells,10)
            rqvcuten_400hPa(1:nCells) = field_interp(1:nCells,11)
            rqvcuten_500hPa(1:nCells) = field_interp(1:nCells,12)
            rqvcuten_600hPa(1:nCells) = field_interp(1:nCells,13)
            rqvcuten_700hPa(1:nCells) = field_interp(1:nCells,14)
            rqvcuten_825hPa(1:nCells) = field_interp(1:nCells,15)
            rqvcuten_850hPa(1:nCells) = field_interp(1:nCells,16)
            rqvcuten_875hPa(1:nCells) = field_interp(1:nCells,17)
            rqvcuten_900hPa(1:nCells) = field_interp(1:nCells,18)
            rqvcuten_925hPa(1:nCells) = field_interp(1:nCells,19)
            rqvcuten_950hPa(1:nCells) = field_interp(1:nCells,20)
            rqvcuten_975hPa(1:nCells) = field_interp(1:nCells,21)
            rqvcuten_1000hPa(1:nCells) = field_interp(1:nCells,22)
            ! call mpas_log_write('--- end interpolate rqvcuten:')
        end if

        if (NEED_RTHRATENLW) then
           !... rthratenlw:
            do iCell = 1, nCells
            do k = 1, nVertLevels
               kk = nVertLevels+1-k
               field_in(iCell,kk) = rthratenlw(k,iCell)
            enddo
            enddo
            call interp_tofixed_pressure(nCells,nVertLevels,nIntP,press_in,field_in,press_interp,field_interp)

            rthratenlw_15hPa(1:nCells) = field_interp(1:nCells,1)
            rthratenlw_20hPa(1:nCells) = field_interp(1:nCells,2)
            rthratenlw_30hPa(1:nCells) = field_interp(1:nCells,3)
            rthratenlw_50hPa(1:nCells) = field_interp(1:nCells,4)
            rthratenlw_70hPa(1:nCells) = field_interp(1:nCells,5)
            rthratenlw_100hPa(1:nCells) = field_interp(1:nCells,6)
            rthratenlw_150hPa(1:nCells) = field_interp(1:nCells,7)
            rthratenlw_200hPa(1:nCells) = field_interp(1:nCells,8)
            rthratenlw_250hPa(1:nCells) = field_interp(1:nCells,9)
            rthratenlw_300hPa(1:nCells) = field_interp(1:nCells,10)
            rthratenlw_400hPa(1:nCells) = field_interp(1:nCells,11)
            rthratenlw_500hPa(1:nCells) = field_interp(1:nCells,12)
            rthratenlw_600hPa(1:nCells) = field_interp(1:nCells,13)
            rthratenlw_700hPa(1:nCells) = field_interp(1:nCells,14)
            rthratenlw_825hPa(1:nCells) = field_interp(1:nCells,15)
            rthratenlw_850hPa(1:nCells) = field_interp(1:nCells,16)
            rthratenlw_875hPa(1:nCells) = field_interp(1:nCells,17)
            rthratenlw_900hPa(1:nCells) = field_interp(1:nCells,18)
            rthratenlw_925hPa(1:nCells) = field_interp(1:nCells,19)
            rthratenlw_950hPa(1:nCells) = field_interp(1:nCells,20)
            rthratenlw_975hPa(1:nCells) = field_interp(1:nCells,21)
            rthratenlw_1000hPa(1:nCells) = field_interp(1:nCells,22)
            ! call mpas_log_write('--- end interpolate rthratenlw:')
        end if

        if (NEED_RUCUTEN) then
         !... rucuten:
          do iCell = 1, nCells
          do k = 1, nVertLevels
             kk = nVertLevels+1-k
             field_in(iCell,kk) = rucuten(k,iCell)
          enddo
          enddo
          call interp_tofixed_pressure(nCells,nVertLevels,nIntP,press_in,field_in,press_interp,field_interp)

          rucuten_15hPa(1:nCells) = field_interp(1:nCells,1)
          rucuten_20hPa(1:nCells) = field_interp(1:nCells,2)
          rucuten_30hPa(1:nCells) = field_interp(1:nCells,3)
          rucuten_50hPa(1:nCells) = field_interp(1:nCells,4)
          rucuten_70hPa(1:nCells) = field_interp(1:nCells,5)
          rucuten_100hPa(1:nCells) = field_interp(1:nCells,6)
          rucuten_150hPa(1:nCells) = field_interp(1:nCells,7)
          rucuten_200hPa(1:nCells) = field_interp(1:nCells,8)
          rucuten_250hPa(1:nCells) = field_interp(1:nCells,9)
          rucuten_300hPa(1:nCells) = field_interp(1:nCells,10)
          rucuten_400hPa(1:nCells) = field_interp(1:nCells,11)
          rucuten_500hPa(1:nCells) = field_interp(1:nCells,12)
          rucuten_600hPa(1:nCells) = field_interp(1:nCells,13)
          rucuten_700hPa(1:nCells) = field_interp(1:nCells,14)
          rucuten_825hPa(1:nCells) = field_interp(1:nCells,15)
          rucuten_850hPa(1:nCells) = field_interp(1:nCells,16)
          rucuten_875hPa(1:nCells) = field_interp(1:nCells,17)
          rucuten_900hPa(1:nCells) = field_interp(1:nCells,18)
          rucuten_925hPa(1:nCells) = field_interp(1:nCells,19)
          rucuten_950hPa(1:nCells) = field_interp(1:nCells,20)
          rucuten_975hPa(1:nCells) = field_interp(1:nCells,21)
          rucuten_1000hPa(1:nCells) = field_interp(1:nCells,22)
          ! call mpas_log_write('--- end interpolate rucuten:')
        end if

        if (NEED_RELHUM) then
           !... relative humidity:
            do iCell = 1, nCells
            do k = 1, nVertLevels
               kk = nVertLevels+1-k
               field_in(iCell,kk) = relhum(k,iCell)
            enddo
            enddo
            call interp_tofixed_pressure(nCells,nVertLevels,nIntP,press_in,field_in,press_interp,field_interp)
            relhum_15hPa(1:nCells) = field_interp(1:nCells,1)
            relhum_20hPa(1:nCells) = field_interp(1:nCells,1)
            relhum_30hPa(1:nCells) = field_interp(1:nCells,1)
            relhum_50hPa(1:nCells) = field_interp(1:nCells,1)
            relhum_70hPa(1:nCells) = field_interp(1:nCells,1)
            relhum_100hPa(1:nCells) = field_interp(1:nCells,2)
            relhum_150hPa(1:nCells) = field_interp(1:nCells,3)
            relhum_200hPa(1:nCells) = field_interp(1:nCells,1)
            relhum_250hPa(1:nCells) = field_interp(1:nCells,4)
            relhum_300hPa(1:nCells) = field_interp(1:nCells,5)
            relhum_400hPa(1:nCells) = field_interp(1:nCells,6)
            relhum_500hPa(1:nCells) = field_interp(1:nCells,7)
            relhum_600hPa(1:nCells) = field_interp(1:nCells,8)
            relhum_700hPa(1:nCells) = field_interp(1:nCells,9)
            relhum_825hPa(1:nCells) = field_interp(1:nCells,10)
            relhum_850hPa(1:nCells) = field_interp(1:nCells,11)
            relhum_875hPa(1:nCells) = field_interp(1:nCells,12)
            relhum_900hPa(1:nCells) = field_interp(1:nCells,13)
            relhum_925hPa(1:nCells) = field_interp(1:nCells,14)
            relhum_950hPa(1:nCells) = field_interp(1:nCells,15)
            relhum_975hPa(1:nCells) = field_interp(1:nCells,16)
            relhum_1000hPa(1:nCells) = field_interp(1:nCells,17)            
           ! call mpas_log_write('--- end interpolate relative humidity:')
        end if
       
        if (NEED_KZQ) then
           !... kzq:
            do iCell = 1, nCells
            do k = 1, nVertLevels
               kk = nVertLevels+1-k
               field_in(iCell,kk) = kzq(k,iCell)
            enddo
            enddo
            call interp_tofixed_pressure(nCells,nVertLevels,nIntP,press_in,field_in,press_interp,field_interp)

            kzq_15hPa(1:nCells) = field_interp(1:nCells,1)
            kzq_20hPa(1:nCells) = field_interp(1:nCells,2)
            kzq_30hPa(1:nCells) = field_interp(1:nCells,3)
            kzq_50hPa(1:nCells) = field_interp(1:nCells,4)
            kzq_70hPa(1:nCells) = field_interp(1:nCells,5)
            kzq_100hPa(1:nCells) = field_interp(1:nCells,6)
            kzq_150hPa(1:nCells) = field_interp(1:nCells,7)
            kzq_200hPa(1:nCells) = field_interp(1:nCells,8)
            kzq_250hPa(1:nCells) = field_interp(1:nCells,9)
            kzq_300hPa(1:nCells) = field_interp(1:nCells,10)
            kzq_400hPa(1:nCells) = field_interp(1:nCells,11)
            kzq_500hPa(1:nCells) = field_interp(1:nCells,12)
            kzq_600hPa(1:nCells) = field_interp(1:nCells,13)
            kzq_700hPa(1:nCells) = field_interp(1:nCells,14)
            kzq_825hPa(1:nCells) = field_interp(1:nCells,15)
            kzq_850hPa(1:nCells) = field_interp(1:nCells,16)
            kzq_875hPa(1:nCells) = field_interp(1:nCells,17)
            kzq_900hPa(1:nCells) = field_interp(1:nCells,18)
            kzq_925hPa(1:nCells) = field_interp(1:nCells,19)
            kzq_950hPa(1:nCells) = field_interp(1:nCells,20)
            kzq_975hPa(1:nCells) = field_interp(1:nCells,21)
            kzq_1000hPa(1:nCells) = field_interp(1:nCells,22)
            ! call mpas_log_write('--- end interpolate kzq:')
        end if

        if (NEED_KZH) then
           !... kzh:
            do iCell = 1, nCells
            do k = 1, nVertLevels
               kk = nVertLevels+1-k
               field_in(iCell,kk) = kzh(k,iCell)
            enddo
            enddo
            call interp_tofixed_pressure(nCells,nVertLevels,nIntP,press_in,field_in,press_interp,field_interp)

            kzh_15hPa(1:nCells) = field_interp(1:nCells,1)
            kzh_20hPa(1:nCells) = field_interp(1:nCells,2)
            kzh_30hPa(1:nCells) = field_interp(1:nCells,3)
            kzh_50hPa(1:nCells) = field_interp(1:nCells,4)
            kzh_70hPa(1:nCells) = field_interp(1:nCells,5)
            kzh_100hPa(1:nCells) = field_interp(1:nCells,6)
            kzh_150hPa(1:nCells) = field_interp(1:nCells,7)
            kzh_200hPa(1:nCells) = field_interp(1:nCells,8)
            kzh_250hPa(1:nCells) = field_interp(1:nCells,9)
            kzh_300hPa(1:nCells) = field_interp(1:nCells,10)
            kzh_400hPa(1:nCells) = field_interp(1:nCells,11)
            kzh_500hPa(1:nCells) = field_interp(1:nCells,12)
            kzh_600hPa(1:nCells) = field_interp(1:nCells,13)
            kzh_700hPa(1:nCells) = field_interp(1:nCells,14)
            kzh_825hPa(1:nCells) = field_interp(1:nCells,15)
            kzh_850hPa(1:nCells) = field_interp(1:nCells,16)
            kzh_875hPa(1:nCells) = field_interp(1:nCells,17)
            kzh_900hPa(1:nCells) = field_interp(1:nCells,18)
            kzh_925hPa(1:nCells) = field_interp(1:nCells,19)
            kzh_950hPa(1:nCells) = field_interp(1:nCells,20)
            kzh_975hPa(1:nCells) = field_interp(1:nCells,21)
            kzh_1000hPa(1:nCells) = field_interp(1:nCells,22)
            ! call mpas_log_write('--- end interpolate kzh:')
        end if

        if (NEED_KZM) then
           !... kzm:
            do iCell = 1, nCells
            do k = 1, nVertLevels
               kk = nVertLevels+1-k
               field_in(iCell,kk) = kzm(k,iCell)
            enddo
            enddo
            call interp_tofixed_pressure(nCells,nVertLevels,nIntP,press_in,field_in,press_interp,field_interp)

            kzm_15hPa(1:nCells) = field_interp(1:nCells,1)
            kzm_20hPa(1:nCells) = field_interp(1:nCells,2)
            kzm_30hPa(1:nCells) = field_interp(1:nCells,3)
            kzm_50hPa(1:nCells) = field_interp(1:nCells,4)
            kzm_70hPa(1:nCells) = field_interp(1:nCells,5)
            kzm_100hPa(1:nCells) = field_interp(1:nCells,6)
            kzm_150hPa(1:nCells) = field_interp(1:nCells,7)
            kzm_200hPa(1:nCells) = field_interp(1:nCells,8)
            kzm_250hPa(1:nCells) = field_interp(1:nCells,9)
            kzm_300hPa(1:nCells) = field_interp(1:nCells,10)
            kzm_400hPa(1:nCells) = field_interp(1:nCells,11)
            kzm_500hPa(1:nCells) = field_interp(1:nCells,12)
            kzm_600hPa(1:nCells) = field_interp(1:nCells,13)
            kzm_700hPa(1:nCells) = field_interp(1:nCells,14)
            kzm_825hPa(1:nCells) = field_interp(1:nCells,15)
            kzm_850hPa(1:nCells) = field_interp(1:nCells,16)
            kzm_875hPa(1:nCells) = field_interp(1:nCells,17)
            kzm_900hPa(1:nCells) = field_interp(1:nCells,18)
            kzm_925hPa(1:nCells) = field_interp(1:nCells,19)
            kzm_950hPa(1:nCells) = field_interp(1:nCells,20)
            kzm_975hPa(1:nCells) = field_interp(1:nCells,21)
            kzm_1000hPa(1:nCells) = field_interp(1:nCells,22)
            ! call mpas_log_write('--- end interpolate kzm:')
        end if


        if (NEED_DEWPOINT) then
           !... dewpoint
            do iCell = 1, nCells
            do k = 1, nVertLevels
               kk = nVertLevels+1-k
               field_in(iCell,kk) = dewpoint(k,iCell)
            enddo
            enddo
            call interp_tofixed_pressure(nCells,nVertLevels,nIntP,press_in,field_in,press_interp,field_interp)
            dewpoint_50hPa(1:nCells) = field_interp(1:nCells,1)
            dewpoint_100hPa(1:nCells) = field_interp(1:nCells,2)
            dewpoint_200hPa(1:nCells) = field_interp(1:nCells,3)
            dewpoint_250hPa(1:nCells) = field_interp(1:nCells,4)
            dewpoint_500hPa(1:nCells) = field_interp(1:nCells,5)
            dewpoint_700hPa(1:nCells) = field_interp(1:nCells,6)
            dewpoint_850hPa(1:nCells) = field_interp(1:nCells,7)
            dewpoint_925hPa(1:nCells) = field_interp(1:nCells,8)
           ! call mpas_log_write('--- end interpolate relative humidity:')
        end if
       
        if (NEED_UZONAL) then
           !... u zonal wind:
            do iCell = 1, nCells
            do k = 1, nVertLevels
               kk = nVertLevels+1-k
               field_in(iCell,kk) = uzonal(k,iCell)
            enddo
            enddo
            call interp_tofixed_pressure(nCells,nVertLevels,nIntP,press_in,field_in,press_interp,field_interp)
            uzonal_15hPa(1:nCells)    = field_interp(1:nCells,1)
            uzonal_20hPa(1:nCells)    = field_interp(1:nCells,2)
            uzonal_30hPa(1:nCells)    = field_interp(1:nCells,3)
            uzonal_50hPa(1:nCells)    = field_interp(1:nCells,4)
            uzonal_70hPa(1:nCells)    = field_interp(1:nCells,5)
            uzonal_100hPa(1:nCells)   = field_interp(1:nCells,6)
            uzonal_150hPa(1:nCells)   = field_interp(1:nCells,7)
            uzonal_200hPa(1:nCells)   = field_interp(1:nCells,8)
            uzonal_250hPa(1:nCells)   = field_interp(1:nCells,9)
            uzonal_300hPa(1:nCells)   = field_interp(1:nCells,10)
            uzonal_400hPa(1:nCells)   = field_interp(1:nCells,11)
            uzonal_500hPa(1:nCells)   = field_interp(1:nCells,12)
            uzonal_600hPa(1:nCells)   = field_interp(1:nCells,13)
            uzonal_700hPa(1:nCells)   = field_interp(1:nCells,14)
            uzonal_825hPa(1:nCells)   = field_interp(1:nCells,15)
            uzonal_850hPa(1:nCells)   = field_interp(1:nCells,16)
            uzonal_875hPa(1:nCells)   = field_interp(1:nCells,17)
            uzonal_900hPa(1:nCells)   = field_interp(1:nCells,18)
            uzonal_925hPa(1:nCells)   = field_interp(1:nCells,19)
            uzonal_950hPa(1:nCells)   = field_interp(1:nCells,20)
            uzonal_975hPa(1:nCells)   = field_interp(1:nCells,21)
            uzonal_1000hPa(1:nCells)  = field_interp(1:nCells,22)

                       ! call mpas_log_write('--- end interpolate zonal wind:')
        end if
       
        if (NEED_UMERIDIONAL) then
           !... u meridional wind:
            do iCell = 1, nCells
            do k = 1, nVertLevels
               kk = nVertLevels+1-k
               field_in(iCell,kk) = umeridional(k,iCell)
            enddo
            enddo
            call interp_tofixed_pressure(nCells,nVertLevels,nIntP,press_in,field_in,press_interp,field_interp)
            umeridional_15hPa(1:nCells) = field_interp(1:nCells,1)
            umeridional_20hPa(1:nCells) = field_interp(1:nCells,2)
            umeridional_30hPa(1:nCells) = field_interp(1:nCells,3)
            umeridional_50hPa(1:nCells) = field_interp(1:nCells,4)
            umeridional_70hPa(1:nCells) = field_interp(1:nCells,5)
            umeridional_100hPa(1:nCells) = field_interp(1:nCells,6)
            umeridional_150hPa(1:nCells) = field_interp(1:nCells,7)
            umeridional_200hPa(1:nCells) = field_interp(1:nCells,8)
            umeridional_250hPa(1:nCells) = field_interp(1:nCells,9)
            umeridional_300hPa(1:nCells) = field_interp(1:nCells,10)
            umeridional_400hPa(1:nCells) = field_interp(1:nCells,11)
            umeridional_500hPa(1:nCells) = field_interp(1:nCells,12)
            umeridional_600hPa(1:nCells) = field_interp(1:nCells,13)
            umeridional_700hPa(1:nCells) = field_interp(1:nCells,14)
            umeridional_825hPa(1:nCells) = field_interp(1:nCells,15)
            umeridional_850hPa(1:nCells) = field_interp(1:nCells,16)
            umeridional_875hPa(1:nCells) = field_interp(1:nCells,17)
            umeridional_900hPa(1:nCells) = field_interp(1:nCells,18)
            umeridional_925hPa(1:nCells) = field_interp(1:nCells,19)
            umeridional_950hPa(1:nCells) = field_interp(1:nCells,20)
            umeridional_975hPa(1:nCells) = field_interp(1:nCells,21)
            umeridional_1000hPa(1:nCells) = field_interp(1:nCells,22)
           ! call mpas_log_write('--- end interpolate meridional wind:')
        end if
       
        if (NEED_SPECHUM) then
           !... spechum:
            do iCell = 1, nCells
            do k = 1, nVertLevels
               kk = nVertLevels+1-k
               field_in(iCell,kk) = spechum(k,iCell)
            enddo
            enddo
            call interp_tofixed_pressure(nCells,nVertLevels,nIntP,press_in,field_in,press_interp,field_interp)
            spechum_15hPa(1:nCells) = field_interp(1:nCells,1)
            spechum_20hPa(1:nCells) = field_interp(1:nCells,2)
            spechum_30hPa(1:nCells) = field_interp(1:nCells,3)
            spechum_50hPa(1:nCells) = field_interp(1:nCells,4)
            spechum_70hPa(1:nCells) = field_interp(1:nCells,5)
            spechum_100hPa(1:nCells) = field_interp(1:nCells,6)
            spechum_150hPa(1:nCells) = field_interp(1:nCells,7)
            spechum_200hPa(1:nCells) = field_interp(1:nCells,8)
            spechum_250hPa(1:nCells) = field_interp(1:nCells,9)
            spechum_300hPa(1:nCells) = field_interp(1:nCells,10)
            spechum_400hPa(1:nCells) = field_interp(1:nCells,11)
            spechum_500hPa(1:nCells) = field_interp(1:nCells,12)
            spechum_600hPa(1:nCells) = field_interp(1:nCells,13)
            spechum_700hPa(1:nCells) = field_interp(1:nCells,14)
            spechum_825hPa(1:nCells) = field_interp(1:nCells,15)
            spechum_850hPa(1:nCells) = field_interp(1:nCells,16)
            spechum_875hPa(1:nCells) = field_interp(1:nCells,17)
            spechum_900hPa(1:nCells) = field_interp(1:nCells,18)
            spechum_925hPa(1:nCells) = field_interp(1:nCells,19)
            spechum_950hPa(1:nCells) = field_interp(1:nCells,20)
            spechum_975hPa(1:nCells) = field_interp(1:nCells,21)
            spechum_1000hPa(1:nCells) = field_interp(1:nCells,22)
           ! call mpas_log_write('--- end interpolate spechum:')
        end if



        if(allocated(field_in)) deallocate(field_in)
        if(allocated(press_in)) deallocate(press_in)
       
        if (NEED_W .or. NEED_HEIGHT .or. NEED_GEOPH) then
           !interpolation to fixed pressure levels for fields located at cells centers and at vertical
           !velocity points:
            if(.not.allocated(press_in)) allocate(press_in(nCells,nVertLevelsP1))
            do iCell = 1, nCells
            do k = 1, nVertLevelsP1
               kk = nVertLevelsP1+1-k
               press_in(iCell,kk) = pressure2(k,iCell)
            enddo
            enddo
       
            if(.not.allocated(field_in)) allocate(field_in(nCells,nVertLevelsP1))
            !... height:
            do iCell = 1, nCells
            do k = 1, nVertLevelsP1
               kk = nVertLevelsP1+1-k
               field_in(iCell,kk) = height(k,iCell)
            enddo
            enddo
            call interp_tofixed_pressure(nCells,nVertLevelsP1,nIntP,press_in,field_in,press_interp,field_interp)
            height_15hPa(1:nCells) = field_interp(1:nCells,1)
            height_20hPa(1:nCells) = field_interp(1:nCells,2)
            height_30hPa(1:nCells) = field_interp(1:nCells,3)
            height_50hPa(1:nCells) = field_interp(1:nCells,4)
            height_70hPa(1:nCells) = field_interp(1:nCells,5)
            height_100hPa(1:nCells) = field_interp(1:nCells,6)
            height_150hPa(1:nCells) = field_interp(1:nCells,7)
            height_200hPa(1:nCells) = field_interp(1:nCells,8)
            height_250hPa(1:nCells) = field_interp(1:nCells,9)
            height_300hPa(1:nCells) = field_interp(1:nCells,10)
            height_400hPa(1:nCells) = field_interp(1:nCells,11)
            height_500hPa(1:nCells) = field_interp(1:nCells,12)
            height_600hPa(1:nCells) = field_interp(1:nCells,13)
            height_700hPa(1:nCells) = field_interp(1:nCells,14)
            height_825hPa(1:nCells) = field_interp(1:nCells,15)
            height_850hPa(1:nCells) = field_interp(1:nCells,16)
            height_875hPa(1:nCells) = field_interp(1:nCells,17)
            height_900hPa(1:nCells) = field_interp(1:nCells,18)
            height_925hPa(1:nCells) = field_interp(1:nCells,19)
            height_950hPa(1:nCells) = field_interp(1:nCells,20)
            height_975hPa(1:nCells) = field_interp(1:nCells,21)
            height_1000hPa(1:nCells) = field_interp(1:nCells,22)
           ! call mpas_log_write('--- end interpolate height:')
        
           !... vertical velocity
            do iCell = 1, nCells
            do k = 1, nVertLevelsP1
               kk = nVertLevelsP1+1-k
               field_in(iCell,kk) = vvel(k,iCell)
            enddo
            enddo
            call interp_tofixed_pressure(nCells,nVertLevelsP1,nIntP,press_in,field_in,press_interp,field_interp)
            w_15hPa(1:nCells) = field_interp(1:nCells,1)
            w_20hPa(1:nCells) = field_interp(1:nCells,2)
            w_30hPa(1:nCells) = field_interp(1:nCells,3)
            w_50hPa(1:nCells) = field_interp(1:nCells,4)
            w_70hPa(1:nCells) = field_interp(1:nCells,5)
            w_100hPa(1:nCells) = field_interp(1:nCells,6)
            w_150hPa(1:nCells) = field_interp(1:nCells,7)
            w_200hPa(1:nCells) = field_interp(1:nCells,8)
            w_250hPa(1:nCells) = field_interp(1:nCells,9)
            w_300hPa(1:nCells) = field_interp(1:nCells,10)
            w_400hPa(1:nCells) = field_interp(1:nCells,11)
            w_500hPa(1:nCells) = field_interp(1:nCells,12)
            w_600hPa(1:nCells) = field_interp(1:nCells,13)
            w_700hPa(1:nCells) = field_interp(1:nCells,14)
            w_825hPa(1:nCells) = field_interp(1:nCells,15)
            w_850hPa(1:nCells) = field_interp(1:nCells,16)
            w_875hPa(1:nCells) = field_interp(1:nCells,17)
            w_900hPa(1:nCells) = field_interp(1:nCells,18)
            w_925hPa(1:nCells) = field_interp(1:nCells,19)
            w_950hPa(1:nCells) = field_interp(1:nCells,20)
            w_975hPa(1:nCells) = field_interp(1:nCells,21)
            w_1000hPa(1:nCells) = field_interp(1:nCells,22)
           ! call mpas_log_write('--- end interpolate vertical velocity:')

           !... geoph:

            if(allocated(field_in)) deallocate(field_in)
            if(.not.allocated(field_in)) allocate(field_in(nCells,nVertLevelsP1))

            call compute_geopotential_height(nCells,nVertLevelsP1,latCell,height,geoph) !EGK

            do iCell = 1, nCells
            do k = 1, nVertLevelsP1
               kk = nVertLevelsP1+1-k
               field_in(iCell,kk) = geoph(k,iCell)
            enddo
            enddo
            call interp_tofixed_pressure(nCells,nVertLevelsP1,nIntP,press_in,field_in,press_interp,field_interp)
            geoph_15hPa(1:nCells) = field_interp(1:nCells,1)
            geoph_20hPa(1:nCells) = field_interp(1:nCells,2)
            geoph_30hPa(1:nCells) = field_interp(1:nCells,3)
            geoph_50hPa(1:nCells) = field_interp(1:nCells,4)
            geoph_70hPa(1:nCells) = field_interp(1:nCells,5)
            geoph_100hPa(1:nCells) = field_interp(1:nCells,6)
            geoph_150hPa(1:nCells) = field_interp(1:nCells,7)
            geoph_200hPa(1:nCells) = field_interp(1:nCells,8)
            geoph_250hPa(1:nCells) = field_interp(1:nCells,9)
            geoph_300hPa(1:nCells) = field_interp(1:nCells,10)
            geoph_400hPa(1:nCells) = field_interp(1:nCells,11)
            geoph_500hPa(1:nCells) = field_interp(1:nCells,12)
            geoph_600hPa(1:nCells) = field_interp(1:nCells,13)
            geoph_700hPa(1:nCells) = field_interp(1:nCells,14)
            geoph_825hPa(1:nCells) = field_interp(1:nCells,15)
            geoph_850hPa(1:nCells) = field_interp(1:nCells,16)
            geoph_875hPa(1:nCells) = field_interp(1:nCells,17)
            geoph_900hPa(1:nCells) = field_interp(1:nCells,18)
            geoph_925hPa(1:nCells) = field_interp(1:nCells,19)
            geoph_950hPa(1:nCells) = field_interp(1:nCells,20)
            geoph_975hPa(1:nCells) = field_interp(1:nCells,21)
            geoph_1000hPa(1:nCells) = field_interp(1:nCells,22)

           ! call mpas_log_write('--- end interpolate geoph:')

            if(allocated(field_in)) deallocate(field_in)
            if(allocated(press_in)) deallocate(press_in)
        end if

        if(allocated(field_interp)) deallocate(field_interp)
        if(allocated(press_interp)) deallocate(press_interp)
       
        if (NEED_VORTICITY) then
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
            vorticity_50hPa(1:nVertices) = field_interp(1:nVertices,1)
            vorticity_100hPa(1:nVertices) = field_interp(1:nVertices,2)
            vorticity_200hPa(1:nVertices) = field_interp(1:nVertices,3)
            vorticity_250hPa(1:nVertices) = field_interp(1:nVertices,4)
            vorticity_500hPa(1:nVertices) = field_interp(1:nVertices,5)
            vorticity_700hPa(1:nVertices) = field_interp(1:nVertices,6)
            vorticity_850hPa(1:nVertices) = field_interp(1:nVertices,7)
            vorticity_925hPa(1:nVertices) = field_interp(1:nVertices,8)
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
            allocate(field_interp(nCells, nIsoLevelsT))
            allocate(press_interp(nCells, nIsoLevelsT))
     
            do k=1,nIsoLevelsT
               press_interp(:,k) = t_iso_levels(k)
            end do
     
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
                call interp_tofixed_pressure(nCells, nVertLevels, nIsoLevelsT, press_in, field_in, press_interp, field_interp)
         
                do k=1,nIsoLevelsT
                   t_isobaric(k,1:nCells) = field_interp(1:nCells,k)
                end do
            end if
     
     
            !!!!!!!!!!! Calculate mean temperature in 500 hPa - 300 hPa layer !!!!!!!!!!!
     
            if (need_meanT_500_300) then
                call compute_layer_mean(meanT_500_300, 50000.0_RKIND, 30000.0_RKIND, field_in, press_in)
            end if
     
     
            deallocate(field_in)
            deallocate(field_interp)
            deallocate(press_in)
            deallocate(press_interp)
        end if
     
     
        !!!!!!!!!!! Additional height levels for vortex tracking !!!!!!!!!!!
        if (need_z_isobaric) then
            allocate(field_in(nCells, nVertLevelsP1))
            allocate(press_in(nCells, nVertLevelsP1))
            allocate(field_interp(nCells, nIsoLevelsZ))
            allocate(press_interp(nCells, nIsoLevelsZ))
     
            do k=1,nIsoLevelsZ
               press_interp(:,k) = z_iso_levels(k)
            end do
     
            do iCell=1,nCells
            do k=1,nVertLevelsP1
               kk = nVertLevelsP1+1-k
               field_in(iCell,kk) = height(k,iCell)
            end do
            end do
     
            do iCell=1,nCells
            do k=1,nVertLevelsP1
               kk = nVertLevelsP1+1-k
               press_in(iCell,kk) = pressure2(k,iCell) * 100.0
            end do
            end do
     
            call interp_tofixed_pressure(nCells, nVertLevelsP1, nIsoLevelsZ, press_in, field_in, press_interp, field_interp)
     
            do k=1,nIsoLevelsZ
               z_isobaric(k,1:nCells) = field_interp(1:nCells,k)
            end do
     
            deallocate(field_in)
            deallocate(field_interp)
            deallocate(press_in)
            deallocate(press_interp)
        end if
    
        if(allocated(temperature) ) deallocate(temperature )
        if(allocated(pressure2)   ) deallocate(pressure2   )
        if(allocated(pressure)    ) deallocate(pressure    )
        if(allocated(dewpoint)    ) deallocate(dewpoint )
        if(allocated(spechum)     ) deallocate(spechum  )
        if(allocated(geoph)       ) deallocate(geoph    )

    end subroutine interp_diagnostics


!***********************************************************************
   !
   !  routine compute_geopotential_height
   !
   !> \brief   Convert geometric height to geopotential height
   !>          Adapted from convert_gpsro_bufr.f90 in DART/observations/gps/. 
   !> \author  Soyoung Ha
   !> \date    17 Feb 2017
   !
   !> \details
   !>  Given latitude (in degree), convert geometric height (in meter)
   !>  into geopotential height (in meter).
   !>  
   !>  Input:   
   !>  ncol -- nCells
   !>  nlev -- nIsobaricLevels
   !>  lat  -- latitude [radian]   
   !>  H    -- geometric height [m]
   !>  
   !>  Output:   
   !>  GPH  -- geopotential height [m]
   !>  
   !----------------------------------------------------------------------- 
    subroutine compute_geopotential_height(ncol, nlev, lat, H, GPH)

       implicit none

    !  Input and Output arguments:
       integer,          intent(in):: ncol, nlev
       real(kind=RKIND), intent(in),  dimension(ncol)       :: lat
       real(kind=RKIND), intent(in),  dimension(nlev,ncol)  :: H
       real(kind=RKIND), intent(in out), dimension(nlev,ncol)  :: GPH

    !  Local variables
       integer :: k, iCell
       real(kind=RKIND), dimension(ncol) :: sin2, termr, termg
       real(kind=RKIND) :: semi_major_axis, semi_minor_axis, grav_polar, grav_equator
       real(kind=RKIND) :: earth_omega, grav_constant, flattening, somigliana
       real(kind=RKIND) :: grav_ratio, grav, eccentricity

    !  Parameters below from WGS-84 model software inside GPS receivers.
       parameter(semi_major_axis = 6378.1370e3_RKIND)    ! (m)
       parameter(semi_minor_axis = 6356.7523142e3_RKIND) ! (m)
       parameter(grav_polar = 9.8321849378_RKIND)        ! (m/s2)
       parameter(grav_equator = 9.7803253359_RKIND)      ! (m/s2)

       parameter(earth_omega = 7.292115e-5_RKIND)        ! (rad/s)  
       parameter(grav = 9.80665_RKIND)                 ! (m/s2) WMO std g at 45 deg lat

       parameter(grav_constant = 3.986004418e14_RKIND)   ! (m3/s2)
       parameter(eccentricity = 0.081819_RKIND)        ! unitless

    !  Derived geophysical constants
       parameter(flattening = (semi_major_axis-semi_minor_axis) / semi_major_axis)

       parameter(somigliana = (semi_minor_axis/semi_major_axis)*(grav_polar/grav_equator)-1.0_RKIND)

       parameter(grav_ratio = (earth_omega*earth_omega * &
                               semi_major_axis*semi_major_axis * semi_minor_axis)/grav_constant)


       sin2(:)  = sin(lat(:))**2
       termg(:) = grav_equator * ( (1.0_RKIND+somigliana*sin2(:)) / &
                      sqrt(1.0_RKIND - eccentricity**2 * sin2(:)) )
       termr(:) = semi_major_axis / (1.0_RKIND + flattening + grav_ratio - 2.0_RKIND*flattening*sin2(:))

       do iCell = 1, ncol
          do k = 1, nlev
             GPH(k,iCell) = (termg(iCell)/grav)*((termr(iCell)*H(k,iCell))/(termr(iCell)+H(k,iCell)))
          end do
       end do

    end subroutine compute_geopotential_height


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
