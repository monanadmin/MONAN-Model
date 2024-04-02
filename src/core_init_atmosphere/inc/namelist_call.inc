   function init_atm_setup_namelists(configPool, namelistFilename, dminfo) result(iErr)
      use mpas_derived_types
      use mpas_pool_routines
      use mpas_io_units
      use mpas_abort, only : mpas_dmpar_global_abort
      use mpas_log, only : mpas_log_write
      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      character (len=*), intent(in) :: namelistFilename
      type (dm_info), intent(in) :: dminfo
      integer :: iErr

      integer :: unitNumber
      logical :: nmlExists

      iErr = 0
      unitNumber = 21
      call mpas_log_write('Reading namelist from file '//trim(namelistFilename))
      inquire(file=trim(namelistFilename), exist=nmlExists)
      if ( .not. nmlExists ) then
         call mpas_dmpar_global_abort('ERROR: Namelist file '//trim(namelistFilename)//' does not exist.')
      end if
      open(unitNumber,file=trim(namelistFilename),status='old',form='formatted')

      call init_atm_setup_nmlrec_nhyd_model(configPool, unitNumber, dminfo)
      call init_atm_setup_nmlrec_dimensions(configPool, unitNumber, dminfo)
      call init_atm_setup_nmlrec_data_sources(configPool, unitNumber, dminfo)
      call init_atm_setup_nmlrec_vertical_grid(configPool, unitNumber, dminfo)
      call init_atm_setup_nmlrec_interpolation_control(configPool, unitNumber, dminfo)
      call init_atm_setup_nmlrec_preproc_stages(configPool, unitNumber, dminfo)
      call init_atm_setup_nmlrec_physics(configPool, unitNumber, dminfo)
      call init_atm_setup_nmlrec_io(configPool, unitNumber, dminfo)
      call init_atm_setup_nmlrec_decomposition(configPool, unitNumber, dminfo)

      close(unitNumber)
   end function init_atm_setup_namelists
