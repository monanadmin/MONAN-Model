   function init_atm_setup_derived_dimensions(readDimensions, dimensionPool, configPool) result(iErr)

      use mpas_derived_types
      use mpas_pool_routines
      use mpas_io_units
      use mpas_log, only : mpas_log_write

      implicit none

      type (mpas_pool_type), intent(inout) :: readDimensions !< Input: Pool to pull read dimensions from
      type (mpas_pool_type), intent(inout) :: configPool !< Input: Pool containing namelist options with configs
      type (mpas_pool_type), intent(inout) :: dimensionPool !< Input/Output: Pool to add dimensions into

      integer :: iErr, errLevel

      integer, pointer :: nCells
      integer, pointer :: nEdges
      integer, pointer :: maxEdges
      integer, pointer :: maxEdges2
      integer, pointer :: nVertices
      integer, pointer :: TWO
      integer, pointer :: THREE
      integer, pointer :: vertexDegree
      integer, pointer :: FIFTEEN
      integer, pointer :: TWENTYONE
      integer, pointer :: R3
      integer, pointer :: nVertLevels
      integer, pointer :: nSoilLevels
      integer, pointer :: nFGLevels
      integer, pointer :: nFGSoilLevels
      integer, pointer :: nVertLevelsP1
      integer, pointer :: nMonths
      integer, pointer :: config_nvertlevels
      integer, pointer :: config_nsoillevels
      integer, pointer :: config_nfglevels
      integer, pointer :: config_nfgsoillevels
      integer, pointer :: config_months

      iErr = 0
      errLevel = mpas_pool_get_error_level()
      call mpas_pool_set_error_level(MPAS_POOL_SILENT)

      nullify(config_nvertlevels)
      call mpas_pool_get_config(configPool, 'config_nvertlevels', config_nvertlevels)
      nullify(config_nsoillevels)
      call mpas_pool_get_config(configPool, 'config_nsoillevels', config_nsoillevels)
      nullify(config_nfglevels)
      call mpas_pool_get_config(configPool, 'config_nfglevels', config_nfglevels)
      nullify(config_nfgsoillevels)
      call mpas_pool_get_config(configPool, 'config_nfgsoillevels', config_nfgsoillevels)
      nullify(config_months)
      call mpas_pool_get_config(configPool, 'config_months', config_months)

      nullify(nCells)
      call mpas_pool_get_dimension(dimensionPool, 'nCells', nCells)
      nullify(nEdges)
      call mpas_pool_get_dimension(dimensionPool, 'nEdges', nEdges)
      nullify(maxEdges)
      call mpas_pool_get_dimension(dimensionPool, 'maxEdges', maxEdges)
      nullify(maxEdges2)
      call mpas_pool_get_dimension(dimensionPool, 'maxEdges2', maxEdges2)
      nullify(nVertices)
      call mpas_pool_get_dimension(dimensionPool, 'nVertices', nVertices)
      nullify(TWO)
      call mpas_pool_get_dimension(dimensionPool, 'TWO', TWO)
      nullify(THREE)
      call mpas_pool_get_dimension(dimensionPool, 'THREE', THREE)
      nullify(vertexDegree)
      call mpas_pool_get_dimension(dimensionPool, 'vertexDegree', vertexDegree)
      nullify(FIFTEEN)
      call mpas_pool_get_dimension(dimensionPool, 'FIFTEEN', FIFTEEN)
      nullify(TWENTYONE)
      call mpas_pool_get_dimension(dimensionPool, 'TWENTYONE', TWENTYONE)
      nullify(R3)
      call mpas_pool_get_dimension(dimensionPool, 'R3', R3)
      nullify(nVertLevels)
      call mpas_pool_get_dimension(dimensionPool, 'nVertLevels', nVertLevels)
      nullify(nSoilLevels)
      call mpas_pool_get_dimension(dimensionPool, 'nSoilLevels', nSoilLevels)
      nullify(nFGLevels)
      call mpas_pool_get_dimension(dimensionPool, 'nFGLevels', nFGLevels)
      nullify(nFGSoilLevels)
      call mpas_pool_get_dimension(dimensionPool, 'nFGSoilLevels', nFGSoilLevels)
      nullify(nVertLevelsP1)
      call mpas_pool_get_dimension(dimensionPool, 'nVertLevelsP1', nVertLevelsP1)
      nullify(nMonths)
      call mpas_pool_get_dimension(dimensionPool, 'nMonths', nMonths)

call mpas_log_write('Assigning remaining dimensions from definitions in Registry.xml ...')
      if ( .not. associated(nCells) ) then
         allocate(nCells)
         nCells = MPAS_MISSING_DIM
         call mpas_pool_add_dimension(dimensionPool, 'nCells', nCells)
      end if

      if ( .not. associated(nEdges) ) then
         allocate(nEdges)
         nEdges = MPAS_MISSING_DIM
         call mpas_pool_add_dimension(dimensionPool, 'nEdges', nEdges)
      end if

      if ( .not. associated(maxEdges) ) then
         allocate(maxEdges)
         maxEdges = MPAS_MISSING_DIM
         call mpas_pool_add_dimension(dimensionPool, 'maxEdges', maxEdges)
      end if

      if ( .not. associated(maxEdges2) ) then
         allocate(maxEdges2)
         maxEdges2 = MPAS_MISSING_DIM
         call mpas_pool_add_dimension(dimensionPool, 'maxEdges2', maxEdges2)
      end if

      if ( .not. associated(nVertices) ) then
         allocate(nVertices)
         nVertices = MPAS_MISSING_DIM
         call mpas_pool_add_dimension(dimensionPool, 'nVertices', nVertices)
      end if

      call mpas_pool_get_dimension(dimensionPool, 'TWO', TWO)
      if ( .not. associated(TWO) ) then
         allocate(TWO)
         TWO = 2
call mpas_log_write('       TWO = $i', intArgs=(/2/))
         call mpas_pool_add_dimension(dimensionPool, 'TWO', TWO)
          else if ( TWO == MPAS_MISSING_DIM ) then
         TWO = 2
          end if

      call mpas_pool_get_dimension(dimensionPool, 'THREE', THREE)
      if ( .not. associated(THREE) ) then
         allocate(THREE)
         THREE = 3
call mpas_log_write('       THREE = $i', intArgs=(/3/))
         call mpas_pool_add_dimension(dimensionPool, 'THREE', THREE)
          else if ( THREE == MPAS_MISSING_DIM ) then
         THREE = 3
          end if

      if ( .not. associated(vertexDegree) ) then
         allocate(vertexDegree)
         vertexDegree = MPAS_MISSING_DIM
         call mpas_pool_add_dimension(dimensionPool, 'vertexDegree', vertexDegree)
      end if

      call mpas_pool_get_dimension(dimensionPool, 'FIFTEEN', FIFTEEN)
      if ( .not. associated(FIFTEEN) ) then
         allocate(FIFTEEN)
         FIFTEEN = 15
call mpas_log_write('       FIFTEEN = $i', intArgs=(/15/))
         call mpas_pool_add_dimension(dimensionPool, 'FIFTEEN', FIFTEEN)
          else if ( FIFTEEN == MPAS_MISSING_DIM ) then
         FIFTEEN = 15
          end if

      call mpas_pool_get_dimension(dimensionPool, 'TWENTYONE', TWENTYONE)
      if ( .not. associated(TWENTYONE) ) then
         allocate(TWENTYONE)
         TWENTYONE = 21
call mpas_log_write('       TWENTYONE = $i', intArgs=(/21/))
         call mpas_pool_add_dimension(dimensionPool, 'TWENTYONE', TWENTYONE)
          else if ( TWENTYONE == MPAS_MISSING_DIM ) then
         TWENTYONE = 21
          end if

      call mpas_pool_get_dimension(dimensionPool, 'R3', R3)
      if ( .not. associated(R3) ) then
         allocate(R3)
         R3 = 3
call mpas_log_write('       R3 = $i', intArgs=(/3/))
         call mpas_pool_add_dimension(dimensionPool, 'R3', R3)
          else if ( R3 == MPAS_MISSING_DIM ) then
         R3 = 3
          end if

      call mpas_pool_get_dimension(dimensionPool, 'nVertLevels', nVertLevels)
      if ( .not. associated(nVertLevels) ) then
         allocate(nVertLevels)
         nVertLevels = config_nvertlevels
call mpas_log_write('       nVertLevels = $i (config_nvertlevels)', intArgs=(/config_nvertlevels/))
         call mpas_pool_add_dimension(dimensionPool, 'nVertLevels', nVertLevels)
          else if ( nVertLevels == MPAS_MISSING_DIM ) then
         nVertLevels = config_nvertlevels
          end if

      call mpas_pool_get_dimension(dimensionPool, 'nSoilLevels', nSoilLevels)
      if ( .not. associated(nSoilLevels) ) then
         allocate(nSoilLevels)
         nSoilLevels = config_nsoillevels
call mpas_log_write('       nSoilLevels = $i (config_nsoillevels)', intArgs=(/config_nsoillevels/))
         call mpas_pool_add_dimension(dimensionPool, 'nSoilLevels', nSoilLevels)
          else if ( nSoilLevels == MPAS_MISSING_DIM ) then
         nSoilLevels = config_nsoillevels
          end if

      call mpas_pool_get_dimension(dimensionPool, 'nFGLevels', nFGLevels)
      if ( .not. associated(nFGLevels) ) then
         allocate(nFGLevels)
         nFGLevels = config_nfglevels
call mpas_log_write('       nFGLevels = $i (config_nfglevels)', intArgs=(/config_nfglevels/))
         call mpas_pool_add_dimension(dimensionPool, 'nFGLevels', nFGLevels)
          else if ( nFGLevels == MPAS_MISSING_DIM ) then
         nFGLevels = config_nfglevels
          end if

      call mpas_pool_get_dimension(dimensionPool, 'nFGSoilLevels', nFGSoilLevels)
      if ( .not. associated(nFGSoilLevels) ) then
         allocate(nFGSoilLevels)
         nFGSoilLevels = config_nfgsoillevels
call mpas_log_write('       nFGSoilLevels = $i (config_nfgsoillevels)', intArgs=(/config_nfgsoillevels/))
         call mpas_pool_add_dimension(dimensionPool, 'nFGSoilLevels', nFGSoilLevels)
          else if ( nFGSoilLevels == MPAS_MISSING_DIM ) then
         nFGSoilLevels = config_nfgsoillevels
          end if

      call mpas_pool_get_dimension(dimensionPool, 'nVertLevelsP1', nVertLevelsP1)
      if ( .not. associated(nVertLevelsP1) ) then
         allocate(nVertLevelsP1)
         nVertLevelsP1 = nVertLevels+1
call mpas_log_write('       nVertLevelsP1 = $i', intArgs=(/nVertLevels+1/))
         call mpas_pool_add_dimension(dimensionPool, 'nVertLevelsP1', nVertLevelsP1)
          else if ( nVertLevelsP1 == MPAS_MISSING_DIM ) then
         nVertLevelsP1 = nVertLevels+1
          end if

      call mpas_pool_get_dimension(dimensionPool, 'nMonths', nMonths)
      if ( .not. associated(nMonths) ) then
         allocate(nMonths)
         nMonths = config_months
call mpas_log_write('       nMonths = $i (config_months)', intArgs=(/config_months/))
         call mpas_pool_add_dimension(dimensionPool, 'nMonths', nMonths)
          else if ( nMonths == MPAS_MISSING_DIM ) then
         nMonths = config_months
          end if

      call mpas_log_write(' ')
      call mpas_log_write(' ----- done assigning dimensions from Registry.xml -----')
      call mpas_log_write(' ')
      call mpas_log_write(' ')
      call mpas_pool_set_error_level(errLevel)

   end function init_atm_setup_derived_dimensions


   function init_atm_setup_decomposed_dimensions(block, manager, readDimensions, dimensionPool, totalBlocks) result(iErr)

      use mpas_derived_types
      use mpas_decomp
      use mpas_pool_routines
      use mpas_io_units
      use mpas_abort, only : mpas_dmpar_global_abort
      use mpas_log, only : mpas_log_write

      implicit none

      type (block_type), intent(inout) :: block !< Input: Pointer to block
      type (mpas_streamManager_type), intent(inout) :: manager !< Input: Stream manager
      type (mpas_pool_type), intent(inout) :: readDimensions !< Input: Pool to pull read dimensions from
      type (mpas_pool_type), intent(inout) :: dimensionPool !< Input/Output: Pool to add dimensions into
      integer, intent(in) :: totalBlocks !< Input: Number of blocks

      integer :: iErr
      type (field1DInteger), pointer :: ownedIndices
      procedure (mpas_decomp_function), pointer :: decompFunc


      iErr = 0
      call mpas_log_write('Processing decomposed dimensions ...')

      call mpas_log_write(' ')
      call mpas_log_write(' ----- done processing decomposed dimensions -----')
      call mpas_log_write(' ')
      call mpas_log_write(' ')

   end function init_atm_setup_decomposed_dimensions
