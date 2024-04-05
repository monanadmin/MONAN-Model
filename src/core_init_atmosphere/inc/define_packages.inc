   function init_atm_define_packages(packagePool) result(iErr)
      use mpas_derived_types
      use mpas_pool_routines
      use mpas_io_units
      implicit none
      type (mpas_pool_type), intent(inout) :: packagePool !< Input: MPAS Pool for containing package logicals.

      integer :: iErr

      iErr = 0
      call mpas_pool_add_package(packagePool, 'initial_condsActive', .false.)
      call mpas_pool_add_package(packagePool, 'sfc_updateActive', .false.)
      call mpas_pool_add_package(packagePool, 'lbcsActive', .false.)
      call mpas_pool_add_package(packagePool, 'gwd_stage_inActive', .false.)
      call mpas_pool_add_package(packagePool, 'gwd_stage_outActive', .false.)
      call mpas_pool_add_package(packagePool, 'vertical_stage_inActive', .false.)
      call mpas_pool_add_package(packagePool, 'vertical_stage_outActive', .false.)
      call mpas_pool_add_package(packagePool, 'met_stage_inActive', .false.)
      call mpas_pool_add_package(packagePool, 'met_stage_outActive', .false.)
      call mpas_pool_add_package(packagePool, 'first_guess_fieldActive', .false.)
   end function init_atm_define_packages
