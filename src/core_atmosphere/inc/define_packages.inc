   function atm_define_packages(packagePool) result(iErr)
      use mpas_derived_types
      use mpas_pool_routines
      use mpas_io_units
      implicit none
      type (mpas_pool_type), intent(inout) :: packagePool !< Input: MPAS Pool for containing package logicals.

      integer :: iErr

      iErr = 0
      call mpas_pool_add_package(packagePool, 'mp_kessler_inActive', .false.)
      call mpas_pool_add_package(packagePool, 'mp_thompson_inActive', .false.)
      call mpas_pool_add_package(packagePool, 'mp_wsm6_inActive', .false.)
      call mpas_pool_add_package(packagePool, 'cu_grell_freitas_inActive', .false.)
      call mpas_pool_add_package(packagePool, 'cu_kain_fritsch_inActive', .false.)
      call mpas_pool_add_package(packagePool, 'cu_tiedtke_inActive', .false.)
      call mpas_pool_add_package(packagePool, 'bl_ysu_inActive', .false.)
      call mpas_pool_add_package(packagePool, 'bl_mynn_inActive', .false.)
      call mpas_pool_add_package(packagePool, 'iauActive', .false.)
      call mpas_pool_add_package(packagePool, 'limited_areaActive', .false.)
      call mpas_pool_add_package(packagePool, 'jedi_daActive', .false.)
   end function atm_define_packages
