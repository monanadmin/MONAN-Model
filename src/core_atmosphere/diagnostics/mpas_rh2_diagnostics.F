! 2-m relative humidity diagnostic (rh2)
!
! rh2 [%] is computed from t2m, q2 (mixing ratio, as produced by sf_mynn)
! and surface pressure, using Bolton (1980) saturation vapor pressure over water.
!
module mpas_rh2_diagnostics

    use mpas_derived_types, only : MPAS_pool_type, MPAS_clock_type
    use mpas_kind_types, only : RKIND

    type (MPAS_pool_type), pointer :: mesh
    type (MPAS_pool_type), pointer :: diag
#ifdef DO_PHYSICS
    type (MPAS_pool_type), pointer :: diag_physics
#endif

    type (MPAS_clock_type), pointer :: clock

    public :: rh2_diagnostics_setup, &
              rh2_diagnostics_compute

    private


    contains


    subroutine rh2_diagnostics_setup(all_pools, simulation_clock)

        use mpas_pool_routines, only : mpas_pool_get_subpool

        implicit none

        type (MPAS_pool_type), pointer :: all_pools
        type (MPAS_clock_type), pointer :: simulation_clock

        call mpas_pool_get_subpool(all_pools, 'mesh', mesh)
        call mpas_pool_get_subpool(all_pools, 'diag', diag)
#ifdef DO_PHYSICS
        call mpas_pool_get_subpool(all_pools, 'diag_physics', diag_physics)
#endif

        clock => simulation_clock

    end subroutine rh2_diagnostics_setup


    subroutine rh2_diagnostics_compute()

        use mpas_pool_routines, only : mpas_pool_get_dimension, mpas_pool_get_array

        implicit none

#ifdef DO_PHYSICS
        integer, pointer :: nCellsSolve
        real (kind=RKIND), dimension(:), pointer :: rh2, t2m, q2, psfc

        integer :: iCell
        real (kind=RKIND) :: tc, q, es, e
        real (kind=RKIND), parameter :: fill_value = -999.0_RKIND
        real (kind=RKIND), parameter :: eps = 0.622_RKIND

        nullify(rh2)
        nullify(t2m)
        nullify(q2)
        nullify(psfc)

        call mpas_pool_get_dimension(mesh, 'nCellsSolve', nCellsSolve)
        call mpas_pool_get_array(diag, 'rh2', rh2)
        call mpas_pool_get_array(diag_physics, 't2m', t2m)
        call mpas_pool_get_array(diag_physics, 'q2', q2)
        call mpas_pool_get_array(diag, 'surface_pressure', psfc)

        ! rh2, t2m and q2 only exist when the MYNN or YSU packages are active
        if (.not. associated(rh2)) return
        if (.not. associated(t2m)) return
        if (.not. associated(q2)) return
        if (.not. associated(psfc)) return
        if (size(rh2) < nCellsSolve .or. size(t2m) < nCellsSolve .or. &
            size(q2) < nCellsSolve .or. size(psfc) < nCellsSolve) return

        do iCell = 1, nCellsSolve
            if (t2m(iCell) < 180.0_RKIND .or. psfc(iCell) <= 0.0_RKIND) then
                rh2(iCell) = fill_value
            else
                tc = t2m(iCell) - 273.15_RKIND
                ! q2 can be slightly negative; floor it at zero for this RH diagnostic.
                q  = max(q2(iCell), 0.0_RKIND)
                es = 611.2_RKIND * exp(17.67_RKIND * tc / (tc + 243.5_RKIND))
                e  = q * psfc(iCell) / (eps + q)
                rh2(iCell) = max(0.0_RKIND, min(100.0_RKIND, 100.0_RKIND * e / es))
            end if
        end do
#endif

    end subroutine rh2_diagnostics_compute

end module mpas_rh2_diagnostics
