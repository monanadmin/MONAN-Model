   subroutine init_atm_generate_pool_mesh(block, structPool, dimensionPool, packagePool)
      use mpas_derived_types
      use mpas_pool_routines
      use mpas_io_units
      use mpas_io, only : MPAS_REAL_FILLVAL, MPAS_INT_FILLVAL, MPAS_CHAR_FILLVAL
      implicit none
      type (block_type), intent(inout), pointer :: block
      type (mpas_pool_type), intent(inout) :: structPool
      type (mpas_pool_type), intent(inout) :: dimensionPool
      type (mpas_pool_type), intent(in) :: packagePool

      type (field0DReal), pointer :: r0Ptr
      type (field1DReal), pointer :: r1Ptr
      type (field2DReal), pointer :: r2Ptr
      type (field3DReal), pointer :: r3Ptr
      type (field4DReal), pointer :: r4Ptr
      type (field5DReal), pointer :: r5Ptr
      type (field0DInteger), pointer :: i0Ptr
      type (field1DInteger), pointer :: i1Ptr
      type (field2DInteger), pointer :: i2Ptr
      type (field3DInteger), pointer :: i3Ptr
      type (field0DChar), pointer :: c0Ptr
      type (field1DChar), pointer :: c1Ptr
      type (field0DReal), dimension(:), pointer :: r0aPtr
      type (field1DReal), dimension(:), pointer :: r1aPtr
      type (field2DReal), dimension(:), pointer :: r2aPtr
      type (field3DReal), dimension(:), pointer :: r3aPtr
      type (field4DReal), dimension(:), pointer :: r4aPtr
      type (field5DReal), dimension(:), pointer :: r5aPtr
      type (field0DInteger), dimension(:), pointer :: i0aPtr
      type (field1DInteger), dimension(:), pointer :: i1aPtr
      type (field2DInteger), dimension(:), pointer :: i2aPtr
      type (field3DInteger), dimension(:), pointer :: i3aPtr
      type (field0DChar), dimension(:), pointer :: c0aPtr
      type (field1DChar), dimension(:), pointer :: c1aPtr

      type (mpas_pool_type), pointer :: newSubPool
      integer :: group_counter
      logical :: group_started
      integer :: group_start
      integer :: index_counter
      integer, pointer :: const_index

      logical, pointer :: initial_condsActive
      logical, pointer :: sfc_updateActive
      logical, pointer :: lbcsActive
      logical, pointer :: gwd_stage_inActive
      logical, pointer :: gwd_stage_outActive
      logical, pointer :: vertical_stage_inActive
      logical, pointer :: vertical_stage_outActive
      logical, pointer :: met_stage_inActive
      logical, pointer :: met_stage_outActive
      logical, pointer :: first_guess_fieldActive


      integer :: numConstituents

      nullify(newSubPool)
      group_counter = -1
      group_started = .false.
      group_start = -1
      call mpas_pool_get_package(packagePool, 'initial_condsActive', initial_condsActive)
      call mpas_pool_get_package(packagePool, 'sfc_updateActive', sfc_updateActive)
      call mpas_pool_get_package(packagePool, 'lbcsActive', lbcsActive)
      call mpas_pool_get_package(packagePool, 'gwd_stage_inActive', gwd_stage_inActive)
      call mpas_pool_get_package(packagePool, 'gwd_stage_outActive', gwd_stage_outActive)
      call mpas_pool_get_package(packagePool, 'vertical_stage_inActive', vertical_stage_inActive)
      call mpas_pool_get_package(packagePool, 'vertical_stage_outActive', vertical_stage_outActive)
      call mpas_pool_get_package(packagePool, 'met_stage_inActive', met_stage_inActive)
      call mpas_pool_get_package(packagePool, 'met_stage_outActive', met_stage_outActive)
      call mpas_pool_get_package(packagePool, 'first_guess_fieldActive', first_guess_fieldActive)

      allocate(newSubPool)
      call mpas_pool_create_pool(newSubPool)
      call mpas_pool_add_subpool(structPool, 'mesh', newSubPool)
      call mpas_pool_add_subpool(block % allStructs, 'mesh', newSubPool)

! Define variable latCell
      allocate(r1Ptr)

! Setting up time level 1
      r1Ptr % fieldName = 'latCell'
      r1Ptr % isVarArray = .false.
      r1Ptr % isDecomposed = .true.
      r1Ptr % hasTimeDimension = .false.
      r1Ptr % isPersistent = .true.
      r1Ptr % isActive = .false.
! Setting up dimensions
      r1Ptr % dimNames(1) = 'nCells'
      r1Ptr % defaultValue = 0.0
      allocate(r1Ptr % attLists(1))
      allocate(r1Ptr % attLists(1) % attList)
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'units', 'rad')
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'long_name', 'Latitude of cells')
      r1Ptr % missingValue = MPAS_REAL_FILLVAL
      r1Ptr % block => block

      r1Ptr % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'latCell', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'latCell', r1Ptr)

! Define variable lonCell
      allocate(r1Ptr)

! Setting up time level 1
      r1Ptr % fieldName = 'lonCell'
      r1Ptr % isVarArray = .false.
      r1Ptr % isDecomposed = .true.
      r1Ptr % hasTimeDimension = .false.
      r1Ptr % isPersistent = .true.
      r1Ptr % isActive = .false.
! Setting up dimensions
      r1Ptr % dimNames(1) = 'nCells'
      r1Ptr % defaultValue = 0.0
      allocate(r1Ptr % attLists(1))
      allocate(r1Ptr % attLists(1) % attList)
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'units', 'rad')
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'long_name', 'Longitude of cells')
      r1Ptr % missingValue = MPAS_REAL_FILLVAL
      r1Ptr % block => block

      r1Ptr % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'lonCell', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'lonCell', r1Ptr)

! Define variable xCell
      allocate(r1Ptr)

! Setting up time level 1
      r1Ptr % fieldName = 'xCell'
      r1Ptr % isVarArray = .false.
      r1Ptr % isDecomposed = .true.
      r1Ptr % hasTimeDimension = .false.
      r1Ptr % isPersistent = .true.
      r1Ptr % isActive = .false.
! Setting up dimensions
      r1Ptr % dimNames(1) = 'nCells'
      r1Ptr % defaultValue = 0.0
      allocate(r1Ptr % attLists(1))
      allocate(r1Ptr % attLists(1) % attList)
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'units', 'm')
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'long_name', 'Cartesian x-coordinate of cells')
      r1Ptr % missingValue = MPAS_REAL_FILLVAL
      r1Ptr % block => block

      r1Ptr % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'xCell', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'xCell', r1Ptr)

! Define variable yCell
      allocate(r1Ptr)

! Setting up time level 1
      r1Ptr % fieldName = 'yCell'
      r1Ptr % isVarArray = .false.
      r1Ptr % isDecomposed = .true.
      r1Ptr % hasTimeDimension = .false.
      r1Ptr % isPersistent = .true.
      r1Ptr % isActive = .false.
! Setting up dimensions
      r1Ptr % dimNames(1) = 'nCells'
      r1Ptr % defaultValue = 0.0
      allocate(r1Ptr % attLists(1))
      allocate(r1Ptr % attLists(1) % attList)
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'units', 'm')
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'long_name', 'Cartesian y-coordinate of cells')
      r1Ptr % missingValue = MPAS_REAL_FILLVAL
      r1Ptr % block => block

      r1Ptr % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'yCell', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'yCell', r1Ptr)

! Define variable zCell
      allocate(r1Ptr)

! Setting up time level 1
      r1Ptr % fieldName = 'zCell'
      r1Ptr % isVarArray = .false.
      r1Ptr % isDecomposed = .true.
      r1Ptr % hasTimeDimension = .false.
      r1Ptr % isPersistent = .true.
      r1Ptr % isActive = .false.
! Setting up dimensions
      r1Ptr % dimNames(1) = 'nCells'
      r1Ptr % defaultValue = 0.0
      allocate(r1Ptr % attLists(1))
      allocate(r1Ptr % attLists(1) % attList)
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'units', 'm')
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'long_name', 'Cartesian z-coordinate of cells')
      r1Ptr % missingValue = MPAS_REAL_FILLVAL
      r1Ptr % block => block

      r1Ptr % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'zCell', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'zCell', r1Ptr)

! Define variable indexToCellID
      allocate(i1Ptr)

! Setting up time level 1
      i1Ptr % fieldName = 'indexToCellID'
      i1Ptr % isVarArray = .false.
      i1Ptr % isDecomposed = .true.
      i1Ptr % hasTimeDimension = .false.
      i1Ptr % isPersistent = .true.
      i1Ptr % isActive = .false.
! Setting up dimensions
      i1Ptr % dimNames(1) = 'nCells'
      i1Ptr % defaultValue = 0
      allocate(i1Ptr % attLists(1))
      allocate(i1Ptr % attLists(1) % attList)
      call mpas_add_att(i1Ptr % attLists(1) % attList, 'units', '-')
      call mpas_add_att(i1Ptr % attLists(1) % attList, 'long_name', 'Mapping from local array index to global cell ID')
      i1Ptr % missingValue = MPAS_INT_FILLVAL
      i1Ptr % block => block

      i1Ptr % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'indexToCellID', i1Ptr)
      call mpas_pool_add_field(block % allFields, 'indexToCellID', i1Ptr)

! Define variable latEdge
      allocate(r1Ptr)

! Setting up time level 1
      r1Ptr % fieldName = 'latEdge'
      r1Ptr % isVarArray = .false.
      r1Ptr % isDecomposed = .true.
      r1Ptr % hasTimeDimension = .false.
      r1Ptr % isPersistent = .true.
      r1Ptr % isActive = .false.
! Setting up dimensions
      r1Ptr % dimNames(1) = 'nEdges'
      r1Ptr % defaultValue = 0.0
      allocate(r1Ptr % attLists(1))
      allocate(r1Ptr % attLists(1) % attList)
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'units', 'rad')
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'long_name', 'Latitude of edges')
      r1Ptr % missingValue = MPAS_REAL_FILLVAL
      r1Ptr % block => block

      r1Ptr % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'latEdge', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'latEdge', r1Ptr)

! Define variable lonEdge
      allocate(r1Ptr)

! Setting up time level 1
      r1Ptr % fieldName = 'lonEdge'
      r1Ptr % isVarArray = .false.
      r1Ptr % isDecomposed = .true.
      r1Ptr % hasTimeDimension = .false.
      r1Ptr % isPersistent = .true.
      r1Ptr % isActive = .false.
! Setting up dimensions
      r1Ptr % dimNames(1) = 'nEdges'
      r1Ptr % defaultValue = 0.0
      allocate(r1Ptr % attLists(1))
      allocate(r1Ptr % attLists(1) % attList)
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'units', 'rad')
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'long_name', 'Longitude of edges')
      r1Ptr % missingValue = MPAS_REAL_FILLVAL
      r1Ptr % block => block

      r1Ptr % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'lonEdge', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'lonEdge', r1Ptr)

! Define variable xEdge
      allocate(r1Ptr)

! Setting up time level 1
      r1Ptr % fieldName = 'xEdge'
      r1Ptr % isVarArray = .false.
      r1Ptr % isDecomposed = .true.
      r1Ptr % hasTimeDimension = .false.
      r1Ptr % isPersistent = .true.
      r1Ptr % isActive = .false.
! Setting up dimensions
      r1Ptr % dimNames(1) = 'nEdges'
      r1Ptr % defaultValue = 0.0
      allocate(r1Ptr % attLists(1))
      allocate(r1Ptr % attLists(1) % attList)
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'units', 'm')
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'long_name', 'Cartesian x-coordinate of edges')
      r1Ptr % missingValue = MPAS_REAL_FILLVAL
      r1Ptr % block => block

      r1Ptr % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'xEdge', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'xEdge', r1Ptr)

! Define variable yEdge
      allocate(r1Ptr)

! Setting up time level 1
      r1Ptr % fieldName = 'yEdge'
      r1Ptr % isVarArray = .false.
      r1Ptr % isDecomposed = .true.
      r1Ptr % hasTimeDimension = .false.
      r1Ptr % isPersistent = .true.
      r1Ptr % isActive = .false.
! Setting up dimensions
      r1Ptr % dimNames(1) = 'nEdges'
      r1Ptr % defaultValue = 0.0
      allocate(r1Ptr % attLists(1))
      allocate(r1Ptr % attLists(1) % attList)
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'units', 'm')
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'long_name', 'Cartesian y-coordinate of edges')
      r1Ptr % missingValue = MPAS_REAL_FILLVAL
      r1Ptr % block => block

      r1Ptr % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'yEdge', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'yEdge', r1Ptr)

! Define variable zEdge
      allocate(r1Ptr)

! Setting up time level 1
      r1Ptr % fieldName = 'zEdge'
      r1Ptr % isVarArray = .false.
      r1Ptr % isDecomposed = .true.
      r1Ptr % hasTimeDimension = .false.
      r1Ptr % isPersistent = .true.
      r1Ptr % isActive = .false.
! Setting up dimensions
      r1Ptr % dimNames(1) = 'nEdges'
      r1Ptr % defaultValue = 0.0
      allocate(r1Ptr % attLists(1))
      allocate(r1Ptr % attLists(1) % attList)
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'units', 'm')
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'long_name', 'Cartesian z-coordinate of edges')
      r1Ptr % missingValue = MPAS_REAL_FILLVAL
      r1Ptr % block => block

      r1Ptr % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'zEdge', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'zEdge', r1Ptr)

! Define variable indexToEdgeID
      allocate(i1Ptr)

! Setting up time level 1
      i1Ptr % fieldName = 'indexToEdgeID'
      i1Ptr % isVarArray = .false.
      i1Ptr % isDecomposed = .true.
      i1Ptr % hasTimeDimension = .false.
      i1Ptr % isPersistent = .true.
      i1Ptr % isActive = .false.
! Setting up dimensions
      i1Ptr % dimNames(1) = 'nEdges'
      i1Ptr % defaultValue = 0
      allocate(i1Ptr % attLists(1))
      allocate(i1Ptr % attLists(1) % attList)
      call mpas_add_att(i1Ptr % attLists(1) % attList, 'units', '-')
      call mpas_add_att(i1Ptr % attLists(1) % attList, 'long_name', 'Mapping from local array index to global edge ID')
      i1Ptr % missingValue = MPAS_INT_FILLVAL
      i1Ptr % block => block

      i1Ptr % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'indexToEdgeID', i1Ptr)
      call mpas_pool_add_field(block % allFields, 'indexToEdgeID', i1Ptr)

! Define variable latVertex
      allocate(r1Ptr)

! Setting up time level 1
      r1Ptr % fieldName = 'latVertex'
      r1Ptr % isVarArray = .false.
      r1Ptr % isDecomposed = .true.
      r1Ptr % hasTimeDimension = .false.
      r1Ptr % isPersistent = .true.
      r1Ptr % isActive = .false.
! Setting up dimensions
      r1Ptr % dimNames(1) = 'nVertices'
      r1Ptr % defaultValue = 0.0
      allocate(r1Ptr % attLists(1))
      allocate(r1Ptr % attLists(1) % attList)
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'units', 'rad')
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'long_name', 'Latitude of vertices')
      r1Ptr % missingValue = MPAS_REAL_FILLVAL
      r1Ptr % block => block

      r1Ptr % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'latVertex', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'latVertex', r1Ptr)

! Define variable lonVertex
      allocate(r1Ptr)

! Setting up time level 1
      r1Ptr % fieldName = 'lonVertex'
      r1Ptr % isVarArray = .false.
      r1Ptr % isDecomposed = .true.
      r1Ptr % hasTimeDimension = .false.
      r1Ptr % isPersistent = .true.
      r1Ptr % isActive = .false.
! Setting up dimensions
      r1Ptr % dimNames(1) = 'nVertices'
      r1Ptr % defaultValue = 0.0
      allocate(r1Ptr % attLists(1))
      allocate(r1Ptr % attLists(1) % attList)
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'units', 'rad')
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'long_name', 'Longitude of vertices')
      r1Ptr % missingValue = MPAS_REAL_FILLVAL
      r1Ptr % block => block

      r1Ptr % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'lonVertex', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'lonVertex', r1Ptr)

! Define variable xVertex
      allocate(r1Ptr)

! Setting up time level 1
      r1Ptr % fieldName = 'xVertex'
      r1Ptr % isVarArray = .false.
      r1Ptr % isDecomposed = .true.
      r1Ptr % hasTimeDimension = .false.
      r1Ptr % isPersistent = .true.
      r1Ptr % isActive = .false.
! Setting up dimensions
      r1Ptr % dimNames(1) = 'nVertices'
      r1Ptr % defaultValue = 0.0
      allocate(r1Ptr % attLists(1))
      allocate(r1Ptr % attLists(1) % attList)
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'units', 'm')
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'long_name', 'Cartesian x-coordinate of vertices')
      r1Ptr % missingValue = MPAS_REAL_FILLVAL
      r1Ptr % block => block

      r1Ptr % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'xVertex', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'xVertex', r1Ptr)

! Define variable yVertex
      allocate(r1Ptr)

! Setting up time level 1
      r1Ptr % fieldName = 'yVertex'
      r1Ptr % isVarArray = .false.
      r1Ptr % isDecomposed = .true.
      r1Ptr % hasTimeDimension = .false.
      r1Ptr % isPersistent = .true.
      r1Ptr % isActive = .false.
! Setting up dimensions
      r1Ptr % dimNames(1) = 'nVertices'
      r1Ptr % defaultValue = 0.0
      allocate(r1Ptr % attLists(1))
      allocate(r1Ptr % attLists(1) % attList)
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'units', 'm')
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'long_name', 'Cartesian y-coordinate of vertices')
      r1Ptr % missingValue = MPAS_REAL_FILLVAL
      r1Ptr % block => block

      r1Ptr % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'yVertex', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'yVertex', r1Ptr)

! Define variable zVertex
      allocate(r1Ptr)

! Setting up time level 1
      r1Ptr % fieldName = 'zVertex'
      r1Ptr % isVarArray = .false.
      r1Ptr % isDecomposed = .true.
      r1Ptr % hasTimeDimension = .false.
      r1Ptr % isPersistent = .true.
      r1Ptr % isActive = .false.
! Setting up dimensions
      r1Ptr % dimNames(1) = 'nVertices'
      r1Ptr % defaultValue = 0.0
      allocate(r1Ptr % attLists(1))
      allocate(r1Ptr % attLists(1) % attList)
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'units', 'm')
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'long_name', 'Cartesian z-coordinate of vertices')
      r1Ptr % missingValue = MPAS_REAL_FILLVAL
      r1Ptr % block => block

      r1Ptr % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'zVertex', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'zVertex', r1Ptr)

! Define variable indexToVertexID
      allocate(i1Ptr)

! Setting up time level 1
      i1Ptr % fieldName = 'indexToVertexID'
      i1Ptr % isVarArray = .false.
      i1Ptr % isDecomposed = .true.
      i1Ptr % hasTimeDimension = .false.
      i1Ptr % isPersistent = .true.
      i1Ptr % isActive = .false.
! Setting up dimensions
      i1Ptr % dimNames(1) = 'nVertices'
      i1Ptr % defaultValue = 0
      allocate(i1Ptr % attLists(1))
      allocate(i1Ptr % attLists(1) % attList)
      call mpas_add_att(i1Ptr % attLists(1) % attList, 'units', '-')
      call mpas_add_att(i1Ptr % attLists(1) % attList, 'long_name', 'Mapping from local array index to global vertex ID')
      i1Ptr % missingValue = MPAS_INT_FILLVAL
      i1Ptr % block => block

      i1Ptr % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'indexToVertexID', i1Ptr)
      call mpas_pool_add_field(block % allFields, 'indexToVertexID', i1Ptr)

! Define variable cellsOnEdge
      allocate(i2Ptr)

! Setting up time level 1
      i2Ptr % fieldName = 'cellsOnEdge'
      i2Ptr % isVarArray = .false.
      i2Ptr % isDecomposed = .true.
      i2Ptr % hasTimeDimension = .false.
      i2Ptr % isPersistent = .true.
      i2Ptr % isActive = .false.
! Setting up dimensions
      i2Ptr % dimNames(1) = 'TWO'
      i2Ptr % dimNames(2) = 'nEdges'
      i2Ptr % defaultValue = 0
      allocate(i2Ptr % attLists(1))
      allocate(i2Ptr % attLists(1) % attList)
      call mpas_add_att(i2Ptr % attLists(1) % attList, 'units', '-')
      call mpas_add_att(i2Ptr % attLists(1) % attList, 'long_name', 'IDs of cells divided by an edge')
      i2Ptr % missingValue = MPAS_INT_FILLVAL
      i2Ptr % block => block

      i2Ptr % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'cellsOnEdge', i2Ptr)
      call mpas_pool_add_field(block % allFields, 'cellsOnEdge', i2Ptr)

! Define variable nEdgesOnCell
      allocate(i1Ptr)

! Setting up time level 1
      i1Ptr % fieldName = 'nEdgesOnCell'
      i1Ptr % isVarArray = .false.
      i1Ptr % isDecomposed = .true.
      i1Ptr % hasTimeDimension = .false.
      i1Ptr % isPersistent = .true.
      i1Ptr % isActive = .false.
! Setting up dimensions
      i1Ptr % dimNames(1) = 'nCells'
      i1Ptr % defaultValue = 0
      allocate(i1Ptr % attLists(1))
      allocate(i1Ptr % attLists(1) % attList)
      call mpas_add_att(i1Ptr % attLists(1) % attList, 'units', '-')
      call mpas_add_att(i1Ptr % attLists(1) % attList, 'long_name', 'Number of edges forming the boundary of a cell')
      i1Ptr % missingValue = MPAS_INT_FILLVAL
      i1Ptr % block => block

      i1Ptr % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'nEdgesOnCell', i1Ptr)
      call mpas_pool_add_field(block % allFields, 'nEdgesOnCell', i1Ptr)

! Define variable nEdgesOnEdge
      allocate(i1Ptr)

! Setting up time level 1
      i1Ptr % fieldName = 'nEdgesOnEdge'
      i1Ptr % isVarArray = .false.
      i1Ptr % isDecomposed = .true.
      i1Ptr % hasTimeDimension = .false.
      i1Ptr % isPersistent = .true.
      i1Ptr % isActive = .false.
! Setting up dimensions
      i1Ptr % dimNames(1) = 'nEdges'
      i1Ptr % defaultValue = 0
      allocate(i1Ptr % attLists(1))
      allocate(i1Ptr % attLists(1) % attList)
      call mpas_add_att(i1Ptr % attLists(1) % attList, 'units', '-')
      call mpas_add_att(i1Ptr % attLists(1) % attList, 'long_name', 'Number of edges involved in reconstruction of tangential '&
//'velocity for an edge')
      i1Ptr % missingValue = MPAS_INT_FILLVAL
      i1Ptr % block => block

      i1Ptr % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'nEdgesOnEdge', i1Ptr)
      call mpas_pool_add_field(block % allFields, 'nEdgesOnEdge', i1Ptr)

! Define variable edgesOnCell
      allocate(i2Ptr)

! Setting up time level 1
      i2Ptr % fieldName = 'edgesOnCell'
      i2Ptr % isVarArray = .false.
      i2Ptr % isDecomposed = .true.
      i2Ptr % hasTimeDimension = .false.
      i2Ptr % isPersistent = .true.
      i2Ptr % isActive = .false.
! Setting up dimensions
      i2Ptr % dimNames(1) = 'maxEdges'
      i2Ptr % dimNames(2) = 'nCells'
      i2Ptr % defaultValue = 0
      allocate(i2Ptr % attLists(1))
      allocate(i2Ptr % attLists(1) % attList)
      call mpas_add_att(i2Ptr % attLists(1) % attList, 'units', '-')
      call mpas_add_att(i2Ptr % attLists(1) % attList, 'long_name', 'IDs of edges forming the boundary of a cell')
      i2Ptr % missingValue = MPAS_INT_FILLVAL
      i2Ptr % block => block

      i2Ptr % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'edgesOnCell', i2Ptr)
      call mpas_pool_add_field(block % allFields, 'edgesOnCell', i2Ptr)

! Define variable edgesOnEdge
      allocate(i2Ptr)

! Setting up time level 1
      i2Ptr % fieldName = 'edgesOnEdge'
      i2Ptr % isVarArray = .false.
      i2Ptr % isDecomposed = .true.
      i2Ptr % hasTimeDimension = .false.
      i2Ptr % isPersistent = .true.
      i2Ptr % isActive = .false.
! Setting up dimensions
      i2Ptr % dimNames(1) = 'maxEdges2'
      i2Ptr % dimNames(2) = 'nEdges'
      i2Ptr % defaultValue = 0
      allocate(i2Ptr % attLists(1))
      allocate(i2Ptr % attLists(1) % attList)
      call mpas_add_att(i2Ptr % attLists(1) % attList, 'units', '-')
      call mpas_add_att(i2Ptr % attLists(1) % attList, 'long_name', 'IDs of edges involved in reconstruction of tangential '&
//'velocity for an edge')
      i2Ptr % missingValue = MPAS_INT_FILLVAL
      i2Ptr % block => block

      i2Ptr % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'edgesOnEdge', i2Ptr)
      call mpas_pool_add_field(block % allFields, 'edgesOnEdge', i2Ptr)

! Define variable weightsOnEdge
      allocate(r2Ptr)

! Setting up time level 1
      r2Ptr % fieldName = 'weightsOnEdge'
      r2Ptr % isVarArray = .false.
      r2Ptr % isDecomposed = .true.
      r2Ptr % hasTimeDimension = .false.
      r2Ptr % isPersistent = .true.
      r2Ptr % isActive = .false.
! Setting up dimensions
      r2Ptr % dimNames(1) = 'maxEdges2'
      r2Ptr % dimNames(2) = 'nEdges'
      r2Ptr % defaultValue = 0.0
      allocate(r2Ptr % attLists(1))
      allocate(r2Ptr % attLists(1) % attList)
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'units', '-')
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'long_name', 'Weights used in reconstruction of tangential velocity for an '&
//'edge')
      r2Ptr % missingValue = MPAS_REAL_FILLVAL
      r2Ptr % block => block

      r2Ptr % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'weightsOnEdge', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'weightsOnEdge', r2Ptr)

! Define variable dvEdge
      allocate(r1Ptr)

! Setting up time level 1
      r1Ptr % fieldName = 'dvEdge'
      r1Ptr % isVarArray = .false.
      r1Ptr % isDecomposed = .true.
      r1Ptr % hasTimeDimension = .false.
      r1Ptr % isPersistent = .true.
      r1Ptr % isActive = .false.
! Setting up dimensions
      r1Ptr % dimNames(1) = 'nEdges'
      r1Ptr % defaultValue = 0.0
      allocate(r1Ptr % attLists(1))
      allocate(r1Ptr % attLists(1) % attList)
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'units', '-')
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'long_name', 'Spherical distance between vertex endpoints of an edge')
      r1Ptr % missingValue = MPAS_REAL_FILLVAL
      r1Ptr % block => block

      r1Ptr % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'dvEdge', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'dvEdge', r1Ptr)

! Define variable dcEdge
      allocate(r1Ptr)

! Setting up time level 1
      r1Ptr % fieldName = 'dcEdge'
      r1Ptr % isVarArray = .false.
      r1Ptr % isDecomposed = .true.
      r1Ptr % hasTimeDimension = .false.
      r1Ptr % isPersistent = .true.
      r1Ptr % isActive = .false.
! Setting up dimensions
      r1Ptr % dimNames(1) = 'nEdges'
      r1Ptr % defaultValue = 0.0
      allocate(r1Ptr % attLists(1))
      allocate(r1Ptr % attLists(1) % attList)
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'units', 'm')
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'long_name', 'Spherical distance between cells separated by an edge')
      r1Ptr % missingValue = MPAS_REAL_FILLVAL
      r1Ptr % block => block

      r1Ptr % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'dcEdge', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'dcEdge', r1Ptr)

! Define variable angleEdge
      allocate(r1Ptr)

! Setting up time level 1
      r1Ptr % fieldName = 'angleEdge'
      r1Ptr % isVarArray = .false.
      r1Ptr % isDecomposed = .true.
      r1Ptr % hasTimeDimension = .false.
      r1Ptr % isPersistent = .true.
      r1Ptr % isActive = .false.
! Setting up dimensions
      r1Ptr % dimNames(1) = 'nEdges'
      r1Ptr % defaultValue = 0.0
      allocate(r1Ptr % attLists(1))
      allocate(r1Ptr % attLists(1) % attList)
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'units', 'rad')
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'long_name', 'Angle between local north and the positive tangential '&
//'direction of an edge')
      r1Ptr % missingValue = MPAS_REAL_FILLVAL
      r1Ptr % block => block

      r1Ptr % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'angleEdge', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'angleEdge', r1Ptr)

! Define variable areaCell
      allocate(r1Ptr)

! Setting up time level 1
      r1Ptr % fieldName = 'areaCell'
      r1Ptr % isVarArray = .false.
      r1Ptr % isDecomposed = .true.
      r1Ptr % hasTimeDimension = .false.
      r1Ptr % isPersistent = .true.
      r1Ptr % isActive = .false.
! Setting up dimensions
      r1Ptr % dimNames(1) = 'nCells'
      r1Ptr % defaultValue = 0.0
      allocate(r1Ptr % attLists(1))
      allocate(r1Ptr % attLists(1) % attList)
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'units', 'm^2')
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'long_name', 'Spherical area of a Voronoi cell')
      r1Ptr % missingValue = MPAS_REAL_FILLVAL
      r1Ptr % block => block

      r1Ptr % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'areaCell', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'areaCell', r1Ptr)

! Define variable areaTriangle
      allocate(r1Ptr)

! Setting up time level 1
      r1Ptr % fieldName = 'areaTriangle'
      r1Ptr % isVarArray = .false.
      r1Ptr % isDecomposed = .true.
      r1Ptr % hasTimeDimension = .false.
      r1Ptr % isPersistent = .true.
      r1Ptr % isActive = .false.
! Setting up dimensions
      r1Ptr % dimNames(1) = 'nVertices'
      r1Ptr % defaultValue = 0.0
      allocate(r1Ptr % attLists(1))
      allocate(r1Ptr % attLists(1) % attList)
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'units', 'm^2')
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'long_name', 'Spherical area of a Delaunay triangle')
      r1Ptr % missingValue = MPAS_REAL_FILLVAL
      r1Ptr % block => block

      r1Ptr % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'areaTriangle', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'areaTriangle', r1Ptr)

! Define variable edgeNormalVectors
      allocate(r2Ptr)

! Setting up time level 1
      r2Ptr % fieldName = 'edgeNormalVectors'
      r2Ptr % isVarArray = .false.
      r2Ptr % isDecomposed = .true.
      r2Ptr % hasTimeDimension = .false.
      r2Ptr % isPersistent = .true.
      r2Ptr % isActive = .false.
! Setting up dimensions
      r2Ptr % dimNames(1) = 'R3'
      r2Ptr % dimNames(2) = 'nEdges'
      r2Ptr % defaultValue = 0.0
      allocate(r2Ptr % attLists(1))
      allocate(r2Ptr % attLists(1) % attList)
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'long_name', 'Cartesian components of the vector normal to an edge and '&
//'tangential to the surface of the sphere')
      r2Ptr % missingValue = MPAS_REAL_FILLVAL
      r2Ptr % block => block

      r2Ptr % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'edgeNormalVectors', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'edgeNormalVectors', r2Ptr)

! Define variable localVerticalUnitVectors
      allocate(r2Ptr)

! Setting up time level 1
      r2Ptr % fieldName = 'localVerticalUnitVectors'
      r2Ptr % isVarArray = .false.
      r2Ptr % isDecomposed = .true.
      r2Ptr % hasTimeDimension = .false.
      r2Ptr % isPersistent = .true.
      r2Ptr % isActive = .false.
! Setting up dimensions
      r2Ptr % dimNames(1) = 'R3'
      r2Ptr % dimNames(2) = 'nCells'
      r2Ptr % defaultValue = 0.0
      allocate(r2Ptr % attLists(1))
      allocate(r2Ptr % attLists(1) % attList)
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'long_name', 'Cartesian components of the vector pointing in the local '&
//'vertical direction for a cell')
      r2Ptr % missingValue = MPAS_REAL_FILLVAL
      r2Ptr % block => block

      r2Ptr % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'localVerticalUnitVectors', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'localVerticalUnitVectors', r2Ptr)

! Define variable cellTangentPlane
      allocate(r3Ptr)

! Setting up time level 1
      r3Ptr % fieldName = 'cellTangentPlane'
      r3Ptr % isVarArray = .false.
      r3Ptr % isDecomposed = .true.
      r3Ptr % hasTimeDimension = .false.
      r3Ptr % isPersistent = .true.
      r3Ptr % isActive = .false.
! Setting up dimensions
      r3Ptr % dimNames(1) = 'R3'
      r3Ptr % dimNames(2) = 'TWO'
      r3Ptr % dimNames(3) = 'nCells'
      r3Ptr % defaultValue = 0.0
      allocate(r3Ptr % attLists(1))
      allocate(r3Ptr % attLists(1) % attList)
      call mpas_add_att(r3Ptr % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r3Ptr % attLists(1) % attList, 'long_name', 'Components of a pair of vectors defining the tangent plane at '&
//'a cell')
      r3Ptr % missingValue = MPAS_REAL_FILLVAL
      r3Ptr % block => block

      r3Ptr % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'cellTangentPlane', r3Ptr)
      call mpas_pool_add_field(block % allFields, 'cellTangentPlane', r3Ptr)

! Define variable cellsOnCell
      allocate(i2Ptr)

! Setting up time level 1
      i2Ptr % fieldName = 'cellsOnCell'
      i2Ptr % isVarArray = .false.
      i2Ptr % isDecomposed = .true.
      i2Ptr % hasTimeDimension = .false.
      i2Ptr % isPersistent = .true.
      i2Ptr % isActive = .false.
! Setting up dimensions
      i2Ptr % dimNames(1) = 'maxEdges'
      i2Ptr % dimNames(2) = 'nCells'
      i2Ptr % defaultValue = 0
      allocate(i2Ptr % attLists(1))
      allocate(i2Ptr % attLists(1) % attList)
      call mpas_add_att(i2Ptr % attLists(1) % attList, 'units', '-')
      call mpas_add_att(i2Ptr % attLists(1) % attList, 'long_name', 'IDs of cells neighboring a cell')
      i2Ptr % missingValue = MPAS_INT_FILLVAL
      i2Ptr % block => block

      i2Ptr % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'cellsOnCell', i2Ptr)
      call mpas_pool_add_field(block % allFields, 'cellsOnCell', i2Ptr)

! Define variable verticesOnCell
      allocate(i2Ptr)

! Setting up time level 1
      i2Ptr % fieldName = 'verticesOnCell'
      i2Ptr % isVarArray = .false.
      i2Ptr % isDecomposed = .true.
      i2Ptr % hasTimeDimension = .false.
      i2Ptr % isPersistent = .true.
      i2Ptr % isActive = .false.
! Setting up dimensions
      i2Ptr % dimNames(1) = 'maxEdges'
      i2Ptr % dimNames(2) = 'nCells'
      i2Ptr % defaultValue = 0
      allocate(i2Ptr % attLists(1))
      allocate(i2Ptr % attLists(1) % attList)
      call mpas_add_att(i2Ptr % attLists(1) % attList, 'units', '-')
      call mpas_add_att(i2Ptr % attLists(1) % attList, 'long_name', 'IDs of vertices (corner points) of a cell')
      i2Ptr % missingValue = MPAS_INT_FILLVAL
      i2Ptr % block => block

      i2Ptr % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'verticesOnCell', i2Ptr)
      call mpas_pool_add_field(block % allFields, 'verticesOnCell', i2Ptr)

! Define variable verticesOnEdge
      allocate(i2Ptr)

! Setting up time level 1
      i2Ptr % fieldName = 'verticesOnEdge'
      i2Ptr % isVarArray = .false.
      i2Ptr % isDecomposed = .true.
      i2Ptr % hasTimeDimension = .false.
      i2Ptr % isPersistent = .true.
      i2Ptr % isActive = .false.
! Setting up dimensions
      i2Ptr % dimNames(1) = 'TWO'
      i2Ptr % dimNames(2) = 'nEdges'
      i2Ptr % defaultValue = 0
      allocate(i2Ptr % attLists(1))
      allocate(i2Ptr % attLists(1) % attList)
      call mpas_add_att(i2Ptr % attLists(1) % attList, 'units', '-')
      call mpas_add_att(i2Ptr % attLists(1) % attList, 'long_name', 'IDs of the two vertex endpoints of an edge')
      i2Ptr % missingValue = MPAS_INT_FILLVAL
      i2Ptr % block => block

      i2Ptr % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'verticesOnEdge', i2Ptr)
      call mpas_pool_add_field(block % allFields, 'verticesOnEdge', i2Ptr)

! Define variable edgesOnVertex
      allocate(i2Ptr)

! Setting up time level 1
      i2Ptr % fieldName = 'edgesOnVertex'
      i2Ptr % isVarArray = .false.
      i2Ptr % isDecomposed = .true.
      i2Ptr % hasTimeDimension = .false.
      i2Ptr % isPersistent = .true.
      i2Ptr % isActive = .false.
! Setting up dimensions
      i2Ptr % dimNames(1) = 'vertexDegree'
      i2Ptr % dimNames(2) = 'nVertices'
      i2Ptr % defaultValue = 0
      allocate(i2Ptr % attLists(1))
      allocate(i2Ptr % attLists(1) % attList)
      call mpas_add_att(i2Ptr % attLists(1) % attList, 'units', '-')
      call mpas_add_att(i2Ptr % attLists(1) % attList, 'long_name', 'IDs of the edges that meet at a vertex')
      i2Ptr % missingValue = MPAS_INT_FILLVAL
      i2Ptr % block => block

      i2Ptr % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'edgesOnVertex', i2Ptr)
      call mpas_pool_add_field(block % allFields, 'edgesOnVertex', i2Ptr)

! Define variable cellsOnVertex
      allocate(i2Ptr)

! Setting up time level 1
      i2Ptr % fieldName = 'cellsOnVertex'
      i2Ptr % isVarArray = .false.
      i2Ptr % isDecomposed = .true.
      i2Ptr % hasTimeDimension = .false.
      i2Ptr % isPersistent = .true.
      i2Ptr % isActive = .false.
! Setting up dimensions
      i2Ptr % dimNames(1) = 'vertexDegree'
      i2Ptr % dimNames(2) = 'nVertices'
      i2Ptr % defaultValue = 0
      allocate(i2Ptr % attLists(1))
      allocate(i2Ptr % attLists(1) % attList)
      call mpas_add_att(i2Ptr % attLists(1) % attList, 'units', '-')
      call mpas_add_att(i2Ptr % attLists(1) % attList, 'long_name', 'IDs of the cells that meet at a vertex')
      i2Ptr % missingValue = MPAS_INT_FILLVAL
      i2Ptr % block => block

      i2Ptr % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'cellsOnVertex', i2Ptr)
      call mpas_pool_add_field(block % allFields, 'cellsOnVertex', i2Ptr)

! Define variable kiteAreasOnVertex
      allocate(r2Ptr)

! Setting up time level 1
      r2Ptr % fieldName = 'kiteAreasOnVertex'
      r2Ptr % isVarArray = .false.
      r2Ptr % isDecomposed = .true.
      r2Ptr % hasTimeDimension = .false.
      r2Ptr % isPersistent = .true.
      r2Ptr % isActive = .false.
! Setting up dimensions
      r2Ptr % dimNames(1) = 'vertexDegree'
      r2Ptr % dimNames(2) = 'nVertices'
      r2Ptr % defaultValue = 0.0
      allocate(r2Ptr % attLists(1))
      allocate(r2Ptr % attLists(1) % attList)
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'units', 'm^2')
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'long_name', 'Intersection areas between primal (Voronoi) and dual '&
//'(triangular) mesh cells')
      r2Ptr % missingValue = MPAS_REAL_FILLVAL
      r2Ptr % block => block

      r2Ptr % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'kiteAreasOnVertex', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'kiteAreasOnVertex', r2Ptr)

! Define variable fEdge
      allocate(r1Ptr)

! Setting up time level 1
      r1Ptr % fieldName = 'fEdge'
      r1Ptr % isVarArray = .false.
      r1Ptr % isDecomposed = .true.
      r1Ptr % hasTimeDimension = .false.
      r1Ptr % isPersistent = .true.
      r1Ptr % isActive = .false.
! Setting up dimensions
      r1Ptr % dimNames(1) = 'nEdges'
      r1Ptr % defaultValue = 0.0
      allocate(r1Ptr % attLists(1))
      allocate(r1Ptr % attLists(1) % attList)
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'long_name', 'Coriolis parameter at an edge')
      r1Ptr % missingValue = MPAS_REAL_FILLVAL
      r1Ptr % block => block

      r1Ptr % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'fEdge', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'fEdge', r1Ptr)

! Define variable fVertex
      allocate(r1Ptr)

! Setting up time level 1
      r1Ptr % fieldName = 'fVertex'
      r1Ptr % isVarArray = .false.
      r1Ptr % isDecomposed = .true.
      r1Ptr % hasTimeDimension = .false.
      r1Ptr % isPersistent = .true.
      r1Ptr % isActive = .false.
! Setting up dimensions
      r1Ptr % dimNames(1) = 'nVertices'
      r1Ptr % defaultValue = 0.0
      allocate(r1Ptr % attLists(1))
      allocate(r1Ptr % attLists(1) % attList)
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'long_name', 'Coriolis parameter at a vertex')
      r1Ptr % missingValue = MPAS_REAL_FILLVAL
      r1Ptr % block => block

      r1Ptr % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'fVertex', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'fVertex', r1Ptr)

! Define variable meshDensity
      allocate(r1Ptr)

! Setting up time level 1
      r1Ptr % fieldName = 'meshDensity'
      r1Ptr % isVarArray = .false.
      r1Ptr % isDecomposed = .true.
      r1Ptr % hasTimeDimension = .false.
      r1Ptr % isPersistent = .true.
      r1Ptr % isActive = .false.
! Setting up dimensions
      r1Ptr % dimNames(1) = 'nCells'
      r1Ptr % defaultValue = 0.0
      allocate(r1Ptr % attLists(1))
      allocate(r1Ptr % attLists(1) % attList)
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'long_name', 'Mesh density function (used when generating the mesh) '&
//'evaluated at a cell')
      r1Ptr % missingValue = MPAS_REAL_FILLVAL
      r1Ptr % block => block

      r1Ptr % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'meshDensity', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'meshDensity', r1Ptr)

! Define variable nominalMinDc
      allocate(r0Ptr)

! Setting up time level 1
      r0Ptr % fieldName = 'nominalMinDc'
      r0Ptr % isVarArray = .false.
      r0Ptr % isDecomposed = .false.
      r0Ptr % hasTimeDimension = .false.
      r0Ptr % defaultValue = 0.0
      r0Ptr % scalar = 0.0
      allocate(r0Ptr % attLists(1))
      allocate(r0Ptr % attLists(1) % attList)
      call mpas_add_att(r0Ptr % attLists(1) % attList, 'units', 'm')
      call mpas_add_att(r0Ptr % attLists(1) % attList, 'long_name', 'Nominal minimum dcEdge value where meshDensity == 1.0')
      r0Ptr % missingValue = MPAS_REAL_FILLVAL
      r0Ptr % block => block

      r0Ptr % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'nominalMinDc', r0Ptr)
      call mpas_pool_add_field(block % allFields, 'nominalMinDc', r0Ptr)

! Define variable cf1
      allocate(r0Ptr)

! Setting up time level 1
      r0Ptr % fieldName = 'cf1'
      r0Ptr % isVarArray = .false.
      r0Ptr % isDecomposed = .false.
      r0Ptr % hasTimeDimension = .false.
      r0Ptr % defaultValue = 0.0
      r0Ptr % scalar = 0.0
      allocate(r0Ptr % attLists(1))
      allocate(r0Ptr % attLists(1) % attList)
      call mpas_add_att(r0Ptr % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r0Ptr % attLists(1) % attList, 'long_name', 'Surface interpolation weight for level k=1 value')
      r0Ptr % missingValue = MPAS_REAL_FILLVAL
      r0Ptr % block => block

      r0Ptr % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'cf1', r0Ptr)
      call mpas_pool_add_field(block % allFields, 'cf1', r0Ptr)

! Define variable cf2
      allocate(r0Ptr)

! Setting up time level 1
      r0Ptr % fieldName = 'cf2'
      r0Ptr % isVarArray = .false.
      r0Ptr % isDecomposed = .false.
      r0Ptr % hasTimeDimension = .false.
      r0Ptr % defaultValue = 0.0
      r0Ptr % scalar = 0.0
      allocate(r0Ptr % attLists(1))
      allocate(r0Ptr % attLists(1) % attList)
      call mpas_add_att(r0Ptr % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r0Ptr % attLists(1) % attList, 'long_name', 'Surface interpolation weight for level k=2 value')
      r0Ptr % missingValue = MPAS_REAL_FILLVAL
      r0Ptr % block => block

      r0Ptr % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'cf2', r0Ptr)
      call mpas_pool_add_field(block % allFields, 'cf2', r0Ptr)

! Define variable cf3
      allocate(r0Ptr)

! Setting up time level 1
      r0Ptr % fieldName = 'cf3'
      r0Ptr % isVarArray = .false.
      r0Ptr % isDecomposed = .false.
      r0Ptr % hasTimeDimension = .false.
      r0Ptr % defaultValue = 0.0
      r0Ptr % scalar = 0.0
      allocate(r0Ptr % attLists(1))
      allocate(r0Ptr % attLists(1) % attList)
      call mpas_add_att(r0Ptr % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r0Ptr % attLists(1) % attList, 'long_name', 'Surface interpolation weight for level k=3 value')
      r0Ptr % missingValue = MPAS_REAL_FILLVAL
      r0Ptr % block => block

      r0Ptr % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'cf3', r0Ptr)
      call mpas_pool_add_field(block % allFields, 'cf3', r0Ptr)

! Define variable ter
      allocate(r1Ptr)

! Setting up time level 1
      r1Ptr % fieldName = 'ter'
      r1Ptr % isVarArray = .false.
      r1Ptr % isDecomposed = .true.
      r1Ptr % hasTimeDimension = .false.
      r1Ptr % isPersistent = .true.
      r1Ptr % isActive = .false.
! Setting up dimensions
      r1Ptr % dimNames(1) = 'nCells'
      r1Ptr % defaultValue = 0.0
      allocate(r1Ptr % attLists(1))
      allocate(r1Ptr % attLists(1) % attList)
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'units', 'm')
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'long_name', 'terrain height')
      r1Ptr % missingValue = MPAS_REAL_FILLVAL
      r1Ptr % block => block

      r1Ptr % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'ter', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'ter', r1Ptr)

! Define variable PHIS
      allocate(r1Ptr)

! Setting up time level 1
      r1Ptr % fieldName = 'PHIS'
      r1Ptr % isVarArray = .false.
      r1Ptr % isDecomposed = .true.
      r1Ptr % hasTimeDimension = .false.
      r1Ptr % isPersistent = .true.
      r1Ptr % isActive = .false.
! Setting up dimensions
      r1Ptr % dimNames(1) = 'nCells'
      r1Ptr % defaultValue = 0.0
      allocate(r1Ptr % attLists(1))
      allocate(r1Ptr % attLists(1) % attList)
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'units', 'm^2 s^{-2}')
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'long_name', 'CAM Surface geopotential')
      r1Ptr % missingValue = MPAS_REAL_FILLVAL
      r1Ptr % block => block

      r1Ptr % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'PHIS', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'PHIS', r1Ptr)

! Define variable landmask
      allocate(i1Ptr)

! Setting up time level 1
      i1Ptr % fieldName = 'landmask'
      i1Ptr % isVarArray = .false.
      i1Ptr % isDecomposed = .true.
      i1Ptr % hasTimeDimension = .false.
      i1Ptr % isPersistent = .true.
      i1Ptr % isActive = .false.
! Setting up dimensions
      i1Ptr % dimNames(1) = 'nCells'
      i1Ptr % defaultValue = 0
      allocate(i1Ptr % attLists(1))
      allocate(i1Ptr % attLists(1) % attList)
      call mpas_add_att(i1Ptr % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(i1Ptr % attLists(1) % attList, 'long_name', 'land-ocean mask (1=land ; 0=ocean)')
      i1Ptr % missingValue = MPAS_INT_FILLVAL
      i1Ptr % block => block

      i1Ptr % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'landmask', i1Ptr)
      call mpas_pool_add_field(block % allFields, 'landmask', i1Ptr)

! Define variable ivgtyp
      allocate(i1Ptr)

! Setting up time level 1
      i1Ptr % fieldName = 'ivgtyp'
      i1Ptr % isVarArray = .false.
      i1Ptr % isDecomposed = .true.
      i1Ptr % hasTimeDimension = .false.
      i1Ptr % isPersistent = .true.
      i1Ptr % isActive = .false.
! Setting up dimensions
      i1Ptr % dimNames(1) = 'nCells'
      i1Ptr % defaultValue = 0
      allocate(i1Ptr % attLists(1))
      allocate(i1Ptr % attLists(1) % attList)
      call mpas_add_att(i1Ptr % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(i1Ptr % attLists(1) % attList, 'long_name', 'dominant vegetation category')
      i1Ptr % missingValue = MPAS_INT_FILLVAL
      i1Ptr % block => block

      i1Ptr % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'lu_index', i1Ptr)
      call mpas_pool_add_field(block % allFields, 'ivgtyp', i1Ptr)

! Define variable mminlu
      allocate(c0Ptr)

! Setting up time level 1
      c0Ptr % fieldName = 'mminlu'
      c0Ptr % isVarArray = .false.
      c0Ptr % isDecomposed = .false.
      c0Ptr % hasTimeDimension = .false.
      c0Ptr % defaultValue = ''
      c0Ptr % scalar = ''
      allocate(c0Ptr % attLists(1))
      allocate(c0Ptr % attLists(1) % attList)
      call mpas_add_att(c0Ptr % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(c0Ptr % attLists(1) % attList, 'long_name', 'land use classification')
      c0Ptr % missingValue = MPAS_CHAR_FILLVAL
      c0Ptr % block => block

      c0Ptr % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'mminlu', c0Ptr)
      call mpas_pool_add_field(block % allFields, 'mminlu', c0Ptr)

! Define variable isltyp
      allocate(i1Ptr)

! Setting up time level 1
      i1Ptr % fieldName = 'isltyp'
      i1Ptr % isVarArray = .false.
      i1Ptr % isDecomposed = .true.
      i1Ptr % hasTimeDimension = .false.
      i1Ptr % isPersistent = .true.
      i1Ptr % isActive = .false.
! Setting up dimensions
      i1Ptr % dimNames(1) = 'nCells'
      i1Ptr % defaultValue = 0
      allocate(i1Ptr % attLists(1))
      allocate(i1Ptr % attLists(1) % attList)
      call mpas_add_att(i1Ptr % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(i1Ptr % attLists(1) % attList, 'long_name', 'dominate soil category')
      i1Ptr % missingValue = MPAS_INT_FILLVAL
      i1Ptr % block => block

      i1Ptr % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'soilcat_top', i1Ptr)
      call mpas_pool_add_field(block % allFields, 'isltyp', i1Ptr)

! Define variable snoalb
      allocate(r1Ptr)

! Setting up time level 1
      r1Ptr % fieldName = 'snoalb'
      r1Ptr % isVarArray = .false.
      r1Ptr % isDecomposed = .true.
      r1Ptr % hasTimeDimension = .false.
      r1Ptr % isPersistent = .true.
      r1Ptr % isActive = .false.
! Setting up dimensions
      r1Ptr % dimNames(1) = 'nCells'
      r1Ptr % defaultValue = 0.0
      allocate(r1Ptr % attLists(1))
      allocate(r1Ptr % attLists(1) % attList)
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'long_name', 'annual maximum snow albedo')
      r1Ptr % missingValue = MPAS_REAL_FILLVAL
      r1Ptr % block => block

      r1Ptr % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'snoalb', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'snoalb', r1Ptr)

! Define variable soiltemp
      allocate(r1Ptr)

! Setting up time level 1
      r1Ptr % fieldName = 'soiltemp'
      r1Ptr % isVarArray = .false.
      r1Ptr % isDecomposed = .true.
      r1Ptr % hasTimeDimension = .false.
      r1Ptr % isPersistent = .true.
      r1Ptr % isActive = .false.
! Setting up dimensions
      r1Ptr % dimNames(1) = 'nCells'
      r1Ptr % defaultValue = 0.0
      allocate(r1Ptr % attLists(1))
      allocate(r1Ptr % attLists(1) % attList)
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'units', 'K')
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'long_name', 'annual mean deep soil temperature')
      r1Ptr % missingValue = MPAS_REAL_FILLVAL
      r1Ptr % block => block

      r1Ptr % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'soiltemp', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'soiltemp', r1Ptr)

! Define variable greenfrac
      allocate(r2Ptr)

! Setting up time level 1
      r2Ptr % fieldName = 'greenfrac'
      r2Ptr % isVarArray = .false.
      r2Ptr % isDecomposed = .true.
      r2Ptr % hasTimeDimension = .false.
      r2Ptr % isPersistent = .true.
      r2Ptr % isActive = .false.
! Setting up dimensions
      r2Ptr % dimNames(1) = 'nMonths'
      r2Ptr % dimNames(2) = 'nCells'
      r2Ptr % defaultValue = 0.0
      allocate(r2Ptr % attLists(1))
      allocate(r2Ptr % attLists(1) % attList)
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'long_name', 'monthly-mean climatological greeness fraction')
      r2Ptr % missingValue = MPAS_REAL_FILLVAL
      r2Ptr % block => block

      r2Ptr % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'greenfrac', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'greenfrac', r2Ptr)

! Define variable shdmin
      allocate(r1Ptr)

! Setting up time level 1
      r1Ptr % fieldName = 'shdmin'
      r1Ptr % isVarArray = .false.
      r1Ptr % isDecomposed = .true.
      r1Ptr % hasTimeDimension = .false.
      r1Ptr % isPersistent = .true.
      r1Ptr % isActive = .false.
! Setting up dimensions
      r1Ptr % dimNames(1) = 'nCells'
      r1Ptr % defaultValue = 0.0
      allocate(r1Ptr % attLists(1))
      allocate(r1Ptr % attLists(1) % attList)
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'long_name', 'minimum fractional coverage of annual green vegetation '&
//'fraction')
      r1Ptr % missingValue = MPAS_REAL_FILLVAL
      r1Ptr % block => block

      r1Ptr % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'shdmin', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'shdmin', r1Ptr)

! Define variable shdmax
      allocate(r1Ptr)

! Setting up time level 1
      r1Ptr % fieldName = 'shdmax'
      r1Ptr % isVarArray = .false.
      r1Ptr % isDecomposed = .true.
      r1Ptr % hasTimeDimension = .false.
      r1Ptr % isPersistent = .true.
      r1Ptr % isActive = .false.
! Setting up dimensions
      r1Ptr % dimNames(1) = 'nCells'
      r1Ptr % defaultValue = 0.0
      allocate(r1Ptr % attLists(1))
      allocate(r1Ptr % attLists(1) % attList)
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'long_name', 'maximum fractional coverage of annual green vegetation '&
//'fraction')
      r1Ptr % missingValue = MPAS_REAL_FILLVAL
      r1Ptr % block => block

      r1Ptr % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'shdmax', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'shdmax', r1Ptr)

! Define variable albedo12m
      allocate(r2Ptr)

! Setting up time level 1
      r2Ptr % fieldName = 'albedo12m'
      r2Ptr % isVarArray = .false.
      r2Ptr % isDecomposed = .true.
      r2Ptr % hasTimeDimension = .false.
      r2Ptr % isPersistent = .true.
      r2Ptr % isActive = .false.
! Setting up dimensions
      r2Ptr % dimNames(1) = 'nMonths'
      r2Ptr % dimNames(2) = 'nCells'
      r2Ptr % defaultValue = 0.0
      allocate(r2Ptr % attLists(1))
      allocate(r2Ptr % attLists(1) % attList)
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'long_name', 'monthly-mean climatological surface albedo')
      r2Ptr % missingValue = MPAS_REAL_FILLVAL
      r2Ptr % block => block

      r2Ptr % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'albedo12m', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'albedo12m', r2Ptr)

! Define variable isice_lu
      allocate(i0Ptr)

! Setting up time level 1
      i0Ptr % fieldName = 'isice_lu'
      i0Ptr % isVarArray = .false.
      i0Ptr % isDecomposed = .false.
      i0Ptr % hasTimeDimension = .false.
      i0Ptr % defaultValue = 24
      i0Ptr % scalar = 24
      allocate(i0Ptr % attLists(1))
      allocate(i0Ptr % attLists(1) % attList)
      call mpas_add_att(i0Ptr % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(i0Ptr % attLists(1) % attList, 'long_name', 'Index category for snow/ice')
      i0Ptr % missingValue = MPAS_INT_FILLVAL
      i0Ptr % block => block

      i0Ptr % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'isice_lu', i0Ptr)
      call mpas_pool_add_field(block % allFields, 'isice_lu', i0Ptr)

! Define variable iswater_lu
      allocate(i0Ptr)

! Setting up time level 1
      i0Ptr % fieldName = 'iswater_lu'
      i0Ptr % isVarArray = .false.
      i0Ptr % isDecomposed = .false.
      i0Ptr % hasTimeDimension = .false.
      i0Ptr % defaultValue = 16
      i0Ptr % scalar = 16
      allocate(i0Ptr % attLists(1))
      allocate(i0Ptr % attLists(1) % attList)
      call mpas_add_att(i0Ptr % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(i0Ptr % attLists(1) % attList, 'long_name', 'Index category for water')
      i0Ptr % missingValue = MPAS_INT_FILLVAL
      i0Ptr % block => block

      i0Ptr % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'iswater_lu', i0Ptr)
      call mpas_pool_add_field(block % allFields, 'iswater_lu', i0Ptr)

! Define variable var2d
      allocate(r1Ptr)

! Setting up time level 1
      r1Ptr % fieldName = 'var2d'
      r1Ptr % isVarArray = .false.
      r1Ptr % isDecomposed = .true.
      r1Ptr % hasTimeDimension = .false.
      r1Ptr % isPersistent = .true.
      r1Ptr % isActive = .false.
! Setting up dimensions
      r1Ptr % dimNames(1) = 'nCells'
      r1Ptr % defaultValue = 0.0
      allocate(r1Ptr % attLists(1))
      allocate(r1Ptr % attLists(1) % attList)
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'units', 'm')
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'long_name', 'standard deviation of subgrid-scale orography')
      r1Ptr % missingValue = MPAS_REAL_FILLVAL
      r1Ptr % block => block

      if (gwd_stage_outActive .or. vertical_stage_outActive .or. met_stage_outActive) then
         r1Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'var2d', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'var2d', r1Ptr)

! Define variable con
      allocate(r1Ptr)

! Setting up time level 1
      r1Ptr % fieldName = 'con'
      r1Ptr % isVarArray = .false.
      r1Ptr % isDecomposed = .true.
      r1Ptr % hasTimeDimension = .false.
      r1Ptr % isPersistent = .true.
      r1Ptr % isActive = .false.
! Setting up dimensions
      r1Ptr % dimNames(1) = 'nCells'
      r1Ptr % defaultValue = 0.0
      allocate(r1Ptr % attLists(1))
      allocate(r1Ptr % attLists(1) % attList)
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'long_name', 'orographic convexity')
      r1Ptr % missingValue = MPAS_REAL_FILLVAL
      r1Ptr % block => block

      if (gwd_stage_outActive .or. vertical_stage_outActive .or. met_stage_outActive) then
         r1Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'con', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'con', r1Ptr)

! Define variable oa1
      allocate(r1Ptr)

! Setting up time level 1
      r1Ptr % fieldName = 'oa1'
      r1Ptr % isVarArray = .false.
      r1Ptr % isDecomposed = .true.
      r1Ptr % hasTimeDimension = .false.
      r1Ptr % isPersistent = .true.
      r1Ptr % isActive = .false.
! Setting up dimensions
      r1Ptr % dimNames(1) = 'nCells'
      r1Ptr % defaultValue = 0.0
      allocate(r1Ptr % attLists(1))
      allocate(r1Ptr % attLists(1) % attList)
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'long_name', 'asymmetry of subgrid-scale orography for westerly flow')
      r1Ptr % missingValue = MPAS_REAL_FILLVAL
      r1Ptr % block => block

      if (gwd_stage_outActive .or. vertical_stage_outActive .or. met_stage_outActive) then
         r1Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'oa1', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'oa1', r1Ptr)

! Define variable oa2
      allocate(r1Ptr)

! Setting up time level 1
      r1Ptr % fieldName = 'oa2'
      r1Ptr % isVarArray = .false.
      r1Ptr % isDecomposed = .true.
      r1Ptr % hasTimeDimension = .false.
      r1Ptr % isPersistent = .true.
      r1Ptr % isActive = .false.
! Setting up dimensions
      r1Ptr % dimNames(1) = 'nCells'
      r1Ptr % defaultValue = 0.0
      allocate(r1Ptr % attLists(1))
      allocate(r1Ptr % attLists(1) % attList)
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'long_name', 'asymmetry of subgrid-scale orography for southerly flow')
      r1Ptr % missingValue = MPAS_REAL_FILLVAL
      r1Ptr % block => block

      if (gwd_stage_outActive .or. vertical_stage_outActive .or. met_stage_outActive) then
         r1Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'oa2', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'oa2', r1Ptr)

! Define variable oa3
      allocate(r1Ptr)

! Setting up time level 1
      r1Ptr % fieldName = 'oa3'
      r1Ptr % isVarArray = .false.
      r1Ptr % isDecomposed = .true.
      r1Ptr % hasTimeDimension = .false.
      r1Ptr % isPersistent = .true.
      r1Ptr % isActive = .false.
! Setting up dimensions
      r1Ptr % dimNames(1) = 'nCells'
      r1Ptr % defaultValue = 0.0
      allocate(r1Ptr % attLists(1))
      allocate(r1Ptr % attLists(1) % attList)
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'long_name', 'asymmetry of subgrid-scale orography for south-westerly flow')
      r1Ptr % missingValue = MPAS_REAL_FILLVAL
      r1Ptr % block => block

      if (gwd_stage_outActive .or. vertical_stage_outActive .or. met_stage_outActive) then
         r1Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'oa3', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'oa3', r1Ptr)

! Define variable oa4
      allocate(r1Ptr)

! Setting up time level 1
      r1Ptr % fieldName = 'oa4'
      r1Ptr % isVarArray = .false.
      r1Ptr % isDecomposed = .true.
      r1Ptr % hasTimeDimension = .false.
      r1Ptr % isPersistent = .true.
      r1Ptr % isActive = .false.
! Setting up dimensions
      r1Ptr % dimNames(1) = 'nCells'
      r1Ptr % defaultValue = 0.0
      allocate(r1Ptr % attLists(1))
      allocate(r1Ptr % attLists(1) % attList)
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'long_name', 'asymmetry of subgrid-scale orography for north-westerly flow')
      r1Ptr % missingValue = MPAS_REAL_FILLVAL
      r1Ptr % block => block

      if (gwd_stage_outActive .or. vertical_stage_outActive .or. met_stage_outActive) then
         r1Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'oa4', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'oa4', r1Ptr)

! Define variable ol1
      allocate(r1Ptr)

! Setting up time level 1
      r1Ptr % fieldName = 'ol1'
      r1Ptr % isVarArray = .false.
      r1Ptr % isDecomposed = .true.
      r1Ptr % hasTimeDimension = .false.
      r1Ptr % isPersistent = .true.
      r1Ptr % isActive = .false.
! Setting up dimensions
      r1Ptr % dimNames(1) = 'nCells'
      r1Ptr % defaultValue = 0.0
      allocate(r1Ptr % attLists(1))
      allocate(r1Ptr % attLists(1) % attList)
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'long_name', 'effective orographic length for westerly flow')
      r1Ptr % missingValue = MPAS_REAL_FILLVAL
      r1Ptr % block => block

      if (gwd_stage_outActive .or. vertical_stage_outActive .or. met_stage_outActive) then
         r1Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'ol1', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'ol1', r1Ptr)

! Define variable ol2
      allocate(r1Ptr)

! Setting up time level 1
      r1Ptr % fieldName = 'ol2'
      r1Ptr % isVarArray = .false.
      r1Ptr % isDecomposed = .true.
      r1Ptr % hasTimeDimension = .false.
      r1Ptr % isPersistent = .true.
      r1Ptr % isActive = .false.
! Setting up dimensions
      r1Ptr % dimNames(1) = 'nCells'
      r1Ptr % defaultValue = 0.0
      allocate(r1Ptr % attLists(1))
      allocate(r1Ptr % attLists(1) % attList)
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'long_name', 'effective orographic length for southerly flow')
      r1Ptr % missingValue = MPAS_REAL_FILLVAL
      r1Ptr % block => block

      if (gwd_stage_outActive .or. vertical_stage_outActive .or. met_stage_outActive) then
         r1Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'ol2', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'ol2', r1Ptr)

! Define variable ol3
      allocate(r1Ptr)

! Setting up time level 1
      r1Ptr % fieldName = 'ol3'
      r1Ptr % isVarArray = .false.
      r1Ptr % isDecomposed = .true.
      r1Ptr % hasTimeDimension = .false.
      r1Ptr % isPersistent = .true.
      r1Ptr % isActive = .false.
! Setting up dimensions
      r1Ptr % dimNames(1) = 'nCells'
      r1Ptr % defaultValue = 0.0
      allocate(r1Ptr % attLists(1))
      allocate(r1Ptr % attLists(1) % attList)
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'units', 'unitles')
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'long_name', 'effective orographic length for south-westerly flow')
      r1Ptr % missingValue = MPAS_REAL_FILLVAL
      r1Ptr % block => block

      if (gwd_stage_outActive .or. vertical_stage_outActive .or. met_stage_outActive) then
         r1Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'ol3', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'ol3', r1Ptr)

! Define variable ol4
      allocate(r1Ptr)

! Setting up time level 1
      r1Ptr % fieldName = 'ol4'
      r1Ptr % isVarArray = .false.
      r1Ptr % isDecomposed = .true.
      r1Ptr % hasTimeDimension = .false.
      r1Ptr % isPersistent = .true.
      r1Ptr % isActive = .false.
! Setting up dimensions
      r1Ptr % dimNames(1) = 'nCells'
      r1Ptr % defaultValue = 0.0
      allocate(r1Ptr % attLists(1))
      allocate(r1Ptr % attLists(1) % attList)
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'long_name', 'effective orographic length for north-westerly flow')
      r1Ptr % missingValue = MPAS_REAL_FILLVAL
      r1Ptr % block => block

      if (gwd_stage_outActive .or. vertical_stage_outActive .or. met_stage_outActive) then
         r1Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'ol4', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'ol4', r1Ptr)

! Define variable hx
      allocate(r2Ptr)

! Setting up time level 1
      r2Ptr % fieldName = 'hx'
      r2Ptr % isVarArray = .false.
      r2Ptr % isDecomposed = .true.
      r2Ptr % hasTimeDimension = .false.
      r2Ptr % isPersistent = .true.
      r2Ptr % isActive = .false.
! Setting up dimensions
      r2Ptr % dimNames(1) = 'nVertLevelsP1'
      r2Ptr % dimNames(2) = 'nCells'
      r2Ptr % defaultValue = 0.0
      allocate(r2Ptr % attLists(1))
      allocate(r2Ptr % attLists(1) % attList)
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'units', 'm')
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'long_name', 'terrain influence in vertical coordinate, $h_s(x,y,\zeta)$ in '&
//'Klemp (MWR 2011)')
      r2Ptr % missingValue = MPAS_REAL_FILLVAL
      r2Ptr % block => block

      if (vertical_stage_outActive .or. met_stage_outActive) then
         r2Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'hx', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'hx', r2Ptr)

! Define variable zgrid
      allocate(r2Ptr)

! Setting up time level 1
      r2Ptr % fieldName = 'zgrid'
      r2Ptr % isVarArray = .false.
      r2Ptr % isDecomposed = .true.
      r2Ptr % hasTimeDimension = .false.
      r2Ptr % isPersistent = .true.
      r2Ptr % isActive = .false.
! Setting up dimensions
      r2Ptr % dimNames(1) = 'nVertLevelsP1'
      r2Ptr % dimNames(2) = 'nCells'
      r2Ptr % defaultValue = 0.0
      allocate(r2Ptr % attLists(1))
      allocate(r2Ptr % attLists(1) % attList)
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'units', 'm MSL')
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'long_name', 'Geometric height of layer interfaces')
      r2Ptr % missingValue = MPAS_REAL_FILLVAL
      r2Ptr % block => block

      if (vertical_stage_outActive .or. met_stage_outActive) then
         r2Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'zgrid', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'zgrid', r2Ptr)

! Define variable rdzw
      allocate(r1Ptr)

! Setting up time level 1
      r1Ptr % fieldName = 'rdzw'
      r1Ptr % isVarArray = .false.
      r1Ptr % isDecomposed = .false.
      r1Ptr % hasTimeDimension = .false.
      r1Ptr % isPersistent = .true.
      r1Ptr % isActive = .false.
! Setting up dimensions
      r1Ptr % dimNames(1) = 'nVertLevels'
      r1Ptr % defaultValue = 0.0
      allocate(r1Ptr % attLists(1))
      allocate(r1Ptr % attLists(1) % attList)
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'long_name', 'Reciprocal dzw')
      r1Ptr % missingValue = MPAS_REAL_FILLVAL
      r1Ptr % block => block

      if (vertical_stage_outActive .or. met_stage_outActive) then
         r1Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'rdzw', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'rdzw', r1Ptr)

! Define variable dzu
      allocate(r1Ptr)

! Setting up time level 1
      r1Ptr % fieldName = 'dzu'
      r1Ptr % isVarArray = .false.
      r1Ptr % isDecomposed = .false.
      r1Ptr % hasTimeDimension = .false.
      r1Ptr % isPersistent = .true.
      r1Ptr % isActive = .false.
! Setting up dimensions
      r1Ptr % dimNames(1) = 'nVertLevels'
      r1Ptr % defaultValue = 0.0
      allocate(r1Ptr % attLists(1))
      allocate(r1Ptr % attLists(1) % attList)
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'long_name', 'd(zeta) at w levels')
      r1Ptr % missingValue = MPAS_REAL_FILLVAL
      r1Ptr % block => block

      if (vertical_stage_outActive .or. met_stage_outActive) then
         r1Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'dzu', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'dzu', r1Ptr)

! Define variable rdzu
      allocate(r1Ptr)

! Setting up time level 1
      r1Ptr % fieldName = 'rdzu'
      r1Ptr % isVarArray = .false.
      r1Ptr % isDecomposed = .false.
      r1Ptr % hasTimeDimension = .false.
      r1Ptr % isPersistent = .true.
      r1Ptr % isActive = .false.
! Setting up dimensions
      r1Ptr % dimNames(1) = 'nVertLevels'
      r1Ptr % defaultValue = 0.0
      allocate(r1Ptr % attLists(1))
      allocate(r1Ptr % attLists(1) % attList)
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'long_name', 'Reciprocal dzu')
      r1Ptr % missingValue = MPAS_REAL_FILLVAL
      r1Ptr % block => block

      if (vertical_stage_outActive .or. met_stage_outActive) then
         r1Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'rdzu', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'rdzu', r1Ptr)

! Define variable fzm
      allocate(r1Ptr)

! Setting up time level 1
      r1Ptr % fieldName = 'fzm'
      r1Ptr % isVarArray = .false.
      r1Ptr % isDecomposed = .false.
      r1Ptr % hasTimeDimension = .false.
      r1Ptr % isPersistent = .true.
      r1Ptr % isActive = .false.
! Setting up dimensions
      r1Ptr % dimNames(1) = 'nVertLevels'
      r1Ptr % defaultValue = 0.0
      allocate(r1Ptr % attLists(1))
      allocate(r1Ptr % attLists(1) % attList)
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'long_name', 'Weight for linear interpolation to w(k) point for u(k) level '&
//'variable')
      r1Ptr % missingValue = MPAS_REAL_FILLVAL
      r1Ptr % block => block

      if (vertical_stage_outActive .or. met_stage_outActive) then
         r1Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'fzm', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'fzm', r1Ptr)

! Define variable fzp
      allocate(r1Ptr)

! Setting up time level 1
      r1Ptr % fieldName = 'fzp'
      r1Ptr % isVarArray = .false.
      r1Ptr % isDecomposed = .false.
      r1Ptr % hasTimeDimension = .false.
      r1Ptr % isPersistent = .true.
      r1Ptr % isActive = .false.
! Setting up dimensions
      r1Ptr % dimNames(1) = 'nVertLevels'
      r1Ptr % defaultValue = 0.0
      allocate(r1Ptr % attLists(1))
      allocate(r1Ptr % attLists(1) % attList)
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'long_name', 'Weight for linear interpolation to w(k) point for u(k-1) '&
//'level variable')
      r1Ptr % missingValue = MPAS_REAL_FILLVAL
      r1Ptr % block => block

      if (vertical_stage_outActive .or. met_stage_outActive) then
         r1Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'fzp', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'fzp', r1Ptr)

! Define variable zxu
      allocate(r2Ptr)

! Setting up time level 1
      r2Ptr % fieldName = 'zxu'
      r2Ptr % isVarArray = .false.
      r2Ptr % isDecomposed = .true.
      r2Ptr % hasTimeDimension = .false.
      r2Ptr % isPersistent = .true.
      r2Ptr % isActive = .false.
! Setting up dimensions
      r2Ptr % dimNames(1) = 'nVertLevels'
      r2Ptr % dimNames(2) = 'nEdges'
      r2Ptr % defaultValue = 0.0
      allocate(r2Ptr % attLists(1))
      allocate(r2Ptr % attLists(1) % attList)
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'long_name', 'dz/dx on horizontal coordinate surfaces at u levels')
      r2Ptr % missingValue = MPAS_REAL_FILLVAL
      r2Ptr % block => block

      if (vertical_stage_outActive .or. met_stage_outActive) then
         r2Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'zxu', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'zxu', r2Ptr)

! Define variable zz
      allocate(r2Ptr)

! Setting up time level 1
      r2Ptr % fieldName = 'zz'
      r2Ptr % isVarArray = .false.
      r2Ptr % isDecomposed = .true.
      r2Ptr % hasTimeDimension = .false.
      r2Ptr % isPersistent = .true.
      r2Ptr % isActive = .false.
! Setting up dimensions
      r2Ptr % dimNames(1) = 'nVertLevels'
      r2Ptr % dimNames(2) = 'nCells'
      r2Ptr % defaultValue = 0.0
      allocate(r2Ptr % attLists(1))
      allocate(r2Ptr % attLists(1) % attList)
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'long_name', 'd(zeta)/dz, vertical metric term')
      r2Ptr % missingValue = MPAS_REAL_FILLVAL
      r2Ptr % block => block

      if (vertical_stage_outActive .or. met_stage_outActive) then
         r2Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'zz', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'zz', r2Ptr)

! Define variable zb
      allocate(r3Ptr)

! Setting up time level 1
      r3Ptr % fieldName = 'zb'
      r3Ptr % isVarArray = .false.
      r3Ptr % isDecomposed = .true.
      r3Ptr % hasTimeDimension = .false.
      r3Ptr % isPersistent = .true.
      r3Ptr % isActive = .false.
! Setting up dimensions
      r3Ptr % dimNames(1) = 'nVertLevelsP1'
      r3Ptr % dimNames(2) = 'TWO'
      r3Ptr % dimNames(3) = 'nEdges'
      r3Ptr % defaultValue = 0.0
      allocate(r3Ptr % attLists(1))
      allocate(r3Ptr % attLists(1) % attList)
      call mpas_add_att(r3Ptr % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r3Ptr % attLists(1) % attList, 'long_name', 'Coefficients for contribution from u to omega diagnosis, '&
//'edge-oriented')
      r3Ptr % missingValue = MPAS_REAL_FILLVAL
      r3Ptr % block => block

      if (vertical_stage_outActive .or. met_stage_outActive) then
         r3Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'zb', r3Ptr)
      call mpas_pool_add_field(block % allFields, 'zb', r3Ptr)

! Define variable zb3
      allocate(r3Ptr)

! Setting up time level 1
      r3Ptr % fieldName = 'zb3'
      r3Ptr % isVarArray = .false.
      r3Ptr % isDecomposed = .true.
      r3Ptr % hasTimeDimension = .false.
      r3Ptr % isPersistent = .true.
      r3Ptr % isActive = .false.
! Setting up dimensions
      r3Ptr % dimNames(1) = 'nVertLevelsP1'
      r3Ptr % dimNames(2) = 'TWO'
      r3Ptr % dimNames(3) = 'nEdges'
      r3Ptr % defaultValue = 0.0
      allocate(r3Ptr % attLists(1))
      allocate(r3Ptr % attLists(1) % attList)
      call mpas_add_att(r3Ptr % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r3Ptr % attLists(1) % attList, 'long_name', 'Coefficients for 3rd-order correction to contribution from u '&
//'to omega diagnosis, edge-oriented')
      r3Ptr % missingValue = MPAS_REAL_FILLVAL
      r3Ptr % block => block

      if (vertical_stage_outActive .or. met_stage_outActive) then
         r3Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'zb3', r3Ptr)
      call mpas_pool_add_field(block % allFields, 'zb3', r3Ptr)

! Define variable dss
      allocate(r2Ptr)

! Setting up time level 1
      r2Ptr % fieldName = 'dss'
      r2Ptr % isVarArray = .false.
      r2Ptr % isDecomposed = .true.
      r2Ptr % hasTimeDimension = .false.
      r2Ptr % isPersistent = .true.
      r2Ptr % isActive = .false.
! Setting up dimensions
      r2Ptr % dimNames(1) = 'nVertLevels'
      r2Ptr % dimNames(2) = 'nCells'
      r2Ptr % defaultValue = 0.0
      allocate(r2Ptr % attLists(1))
      allocate(r2Ptr % attLists(1) % attList)
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'long_name', 'w-damping coefficient')
      r2Ptr % missingValue = MPAS_REAL_FILLVAL
      r2Ptr % block => block

      if (vertical_stage_outActive .or. met_stage_outActive) then
         r2Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'dss', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'dss', r2Ptr)

! Define variable u_init
      allocate(r1Ptr)

! Setting up time level 1
      r1Ptr % fieldName = 'u_init'
      r1Ptr % isVarArray = .false.
      r1Ptr % isDecomposed = .false.
      r1Ptr % hasTimeDimension = .false.
      r1Ptr % isPersistent = .true.
      r1Ptr % isActive = .false.
! Setting up dimensions
      r1Ptr % dimNames(1) = 'nVertLevels'
      r1Ptr % defaultValue = 0.0
      allocate(r1Ptr % attLists(1))
      allocate(r1Ptr % attLists(1) % attList)
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'units', 'm s^{-1}')
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'long_name', 'u reference profile')
      r1Ptr % missingValue = MPAS_REAL_FILLVAL
      r1Ptr % block => block

      if (met_stage_outActive) then
         r1Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'u_init', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'u_init', r1Ptr)

! Define variable v_init
      allocate(r1Ptr)

! Setting up time level 1
      r1Ptr % fieldName = 'v_init'
      r1Ptr % isVarArray = .false.
      r1Ptr % isDecomposed = .false.
      r1Ptr % hasTimeDimension = .false.
      r1Ptr % isPersistent = .true.
      r1Ptr % isActive = .false.
! Setting up dimensions
      r1Ptr % dimNames(1) = 'nVertLevels'
      r1Ptr % defaultValue = 0.0
      allocate(r1Ptr % attLists(1))
      allocate(r1Ptr % attLists(1) % attList)
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'units', 'm s^{-1}')
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'long_name', 'v reference profile')
      r1Ptr % missingValue = MPAS_REAL_FILLVAL
      r1Ptr % block => block

      if (met_stage_outActive) then
         r1Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'v_init', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'v_init', r1Ptr)

! Define variable t_init
      allocate(r2Ptr)

! Setting up time level 1
      r2Ptr % fieldName = 't_init'
      r2Ptr % isVarArray = .false.
      r2Ptr % isDecomposed = .true.
      r2Ptr % hasTimeDimension = .false.
      r2Ptr % isPersistent = .true.
      r2Ptr % isActive = .false.
! Setting up dimensions
      r2Ptr % dimNames(1) = 'nVertLevels'
      r2Ptr % dimNames(2) = 'nCells'
      r2Ptr % defaultValue = 0.0
      allocate(r2Ptr % attLists(1))
      allocate(r2Ptr % attLists(1) % attList)
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'units', 'K')
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'long_name', 'theta reference profile')
      r2Ptr % missingValue = MPAS_REAL_FILLVAL
      r2Ptr % block => block

      if (met_stage_outActive) then
         r2Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 't_init', r2Ptr)
      call mpas_pool_add_field(block % allFields, 't_init', r2Ptr)

! Define variable qv_init
      allocate(r1Ptr)

! Setting up time level 1
      r1Ptr % fieldName = 'qv_init'
      r1Ptr % isVarArray = .false.
      r1Ptr % isDecomposed = .false.
      r1Ptr % hasTimeDimension = .false.
      r1Ptr % isPersistent = .true.
      r1Ptr % isActive = .false.
! Setting up dimensions
      r1Ptr % dimNames(1) = 'nVertLevels'
      r1Ptr % defaultValue = 0.0
      allocate(r1Ptr % attLists(1))
      allocate(r1Ptr % attLists(1) % attList)
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'units', 'kg kg^{-1}')
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'long_name', 'qv reference profile')
      r1Ptr % missingValue = MPAS_REAL_FILLVAL
      r1Ptr % block => block

      if (met_stage_outActive) then
         r1Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'qv_init', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'qv_init', r1Ptr)

! Define variable deriv_two
      allocate(r3Ptr)

! Setting up time level 1
      r3Ptr % fieldName = 'deriv_two'
      r3Ptr % isVarArray = .false.
      r3Ptr % isDecomposed = .true.
      r3Ptr % hasTimeDimension = .false.
      r3Ptr % isPersistent = .true.
      r3Ptr % isActive = .false.
! Setting up dimensions
      r3Ptr % dimNames(1) = 'FIFTEEN'
      r3Ptr % dimNames(2) = 'TWO'
      r3Ptr % dimNames(3) = 'nEdges'
      r3Ptr % defaultValue = 0.0
      allocate(r3Ptr % attLists(1))
      allocate(r3Ptr % attLists(1) % attList)
      call mpas_add_att(r3Ptr % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r3Ptr % attLists(1) % attList, 'long_name', 'weights for cell-centered second derivative, normal to edge, '&
//'for transport scheme')
      r3Ptr % missingValue = MPAS_REAL_FILLVAL
      r3Ptr % block => block

      r3Ptr % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'deriv_two', r3Ptr)
      call mpas_pool_add_field(block % allFields, 'deriv_two', r3Ptr)

! Define variable advCells
      allocate(i2Ptr)

! Setting up time level 1
      i2Ptr % fieldName = 'advCells'
      i2Ptr % isVarArray = .false.
      i2Ptr % isDecomposed = .true.
      i2Ptr % hasTimeDimension = .false.
      i2Ptr % isPersistent = .true.
      i2Ptr % isActive = .false.
! Setting up dimensions
      i2Ptr % dimNames(1) = 'TWENTYONE'
      i2Ptr % dimNames(2) = 'nCells'
      i2Ptr % defaultValue = 0
      allocate(i2Ptr % attLists(1))
      allocate(i2Ptr % attLists(1) % attList)
      call mpas_add_att(i2Ptr % attLists(1) % attList, 'units', '-')
      call mpas_add_att(i2Ptr % attLists(1) % attList, 'long_name', 'cells used in least squares fit of polynomial for advection')
      i2Ptr % missingValue = MPAS_INT_FILLVAL
      i2Ptr % block => block

      i2Ptr % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'advCells', i2Ptr)
      call mpas_pool_add_field(block % allFields, 'advCells', i2Ptr)

! Define variable defc_a
      allocate(r2Ptr)

! Setting up time level 1
      r2Ptr % fieldName = 'defc_a'
      r2Ptr % isVarArray = .false.
      r2Ptr % isDecomposed = .true.
      r2Ptr % hasTimeDimension = .false.
      r2Ptr % isPersistent = .true.
      r2Ptr % isActive = .false.
! Setting up dimensions
      r2Ptr % dimNames(1) = 'maxEdges'
      r2Ptr % dimNames(2) = 'nCells'
      r2Ptr % defaultValue = 0.0
      allocate(r2Ptr % attLists(1))
      allocate(r2Ptr % attLists(1) % attList)
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'long_name', 'Coefficients for computing the off-diagonal components of the '&
//'horizontal deformation')
      r2Ptr % missingValue = MPAS_REAL_FILLVAL
      r2Ptr % block => block

      r2Ptr % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'defc_a', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'defc_a', r2Ptr)

! Define variable defc_b
      allocate(r2Ptr)

! Setting up time level 1
      r2Ptr % fieldName = 'defc_b'
      r2Ptr % isVarArray = .false.
      r2Ptr % isDecomposed = .true.
      r2Ptr % hasTimeDimension = .false.
      r2Ptr % isPersistent = .true.
      r2Ptr % isActive = .false.
! Setting up dimensions
      r2Ptr % dimNames(1) = 'maxEdges'
      r2Ptr % dimNames(2) = 'nCells'
      r2Ptr % defaultValue = 0.0
      allocate(r2Ptr % attLists(1))
      allocate(r2Ptr % attLists(1) % attList)
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'long_name', 'Coefficients for computing the diagonal components of the '&
//'horizontal deformation')
      r2Ptr % missingValue = MPAS_REAL_FILLVAL
      r2Ptr % block => block

      r2Ptr % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'defc_b', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'defc_b', r2Ptr)

! Define variable cell_gradient_coef_x
      allocate(r2Ptr)

! Setting up time level 1
      r2Ptr % fieldName = 'cell_gradient_coef_x'
      r2Ptr % isVarArray = .false.
      r2Ptr % isDecomposed = .true.
      r2Ptr % hasTimeDimension = .false.
      r2Ptr % isPersistent = .true.
      r2Ptr % isActive = .false.
! Setting up dimensions
      r2Ptr % dimNames(1) = 'maxEdges'
      r2Ptr % dimNames(2) = 'nCells'
      r2Ptr % defaultValue = 0.0
      allocate(r2Ptr % attLists(1))
      allocate(r2Ptr % attLists(1) % attList)
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'units', 'm^-1')
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'long_name', 'Coefficients for computing the x (zonal) derivative of a '&
//'cell-centered variable')
      r2Ptr % missingValue = MPAS_REAL_FILLVAL
      r2Ptr % block => block

      r2Ptr % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'cell_gradient_coef_x', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'cell_gradient_coef_x', r2Ptr)

! Define variable cell_gradient_coef_y
      allocate(r2Ptr)

! Setting up time level 1
      r2Ptr % fieldName = 'cell_gradient_coef_y'
      r2Ptr % isVarArray = .false.
      r2Ptr % isDecomposed = .true.
      r2Ptr % hasTimeDimension = .false.
      r2Ptr % isPersistent = .true.
      r2Ptr % isActive = .false.
! Setting up dimensions
      r2Ptr % dimNames(1) = 'maxEdges'
      r2Ptr % dimNames(2) = 'nCells'
      r2Ptr % defaultValue = 0.0
      allocate(r2Ptr % attLists(1))
      allocate(r2Ptr % attLists(1) % attList)
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'units', 'm^-1')
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'long_name', 'Coefficients for computing the y (meridional) derivative of a '&
//'cell-centered variable')
      r2Ptr % missingValue = MPAS_REAL_FILLVAL
      r2Ptr % block => block

      r2Ptr % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'cell_gradient_coef_y', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'cell_gradient_coef_y', r2Ptr)

! Define variable coeffs_reconstruct
      allocate(r3Ptr)

! Setting up time level 1
      r3Ptr % fieldName = 'coeffs_reconstruct'
      r3Ptr % isVarArray = .false.
      r3Ptr % isDecomposed = .true.
      r3Ptr % hasTimeDimension = .false.
      r3Ptr % isPersistent = .true.
      r3Ptr % isActive = .false.
! Setting up dimensions
      r3Ptr % dimNames(1) = 'R3'
      r3Ptr % dimNames(2) = 'maxEdges'
      r3Ptr % dimNames(3) = 'nCells'
      r3Ptr % defaultValue = 0.0
      allocate(r3Ptr % attLists(1))
      allocate(r3Ptr % attLists(1) % attList)
      call mpas_add_att(r3Ptr % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r3Ptr % attLists(1) % attList, 'long_name', 'Coefficients to reconstruct velocity vectors at cell centers')
      r3Ptr % missingValue = MPAS_REAL_FILLVAL
      r3Ptr % block => block

      r3Ptr % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'coeffs_reconstruct', r3Ptr)
      call mpas_pool_add_field(block % allFields, 'coeffs_reconstruct', r3Ptr)

! Define variable bdyMaskCell
      allocate(i1Ptr)

! Setting up time level 1
      i1Ptr % fieldName = 'bdyMaskCell'
      i1Ptr % isVarArray = .false.
      i1Ptr % isDecomposed = .true.
      i1Ptr % hasTimeDimension = .false.
      i1Ptr % isPersistent = .true.
      i1Ptr % isActive = .false.
! Setting up dimensions
      i1Ptr % dimNames(1) = 'nCells'
      i1Ptr % defaultValue = 0
      allocate(i1Ptr % attLists(1))
      allocate(i1Ptr % attLists(1) % attList)
      call mpas_add_att(i1Ptr % attLists(1) % attList, 'units', '-')
      call mpas_add_att(i1Ptr % attLists(1) % attList, 'long_name', 'Indicator of whether a cell is an interior cell, a '&
//'relaxation-zone cell, or a specified-zone cell')
      i1Ptr % missingValue = MPAS_INT_FILLVAL
      i1Ptr % block => block

      i1Ptr % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'bdyMaskCell', i1Ptr)
      call mpas_pool_add_field(block % allFields, 'bdyMaskCell', i1Ptr)

! Define variable bdyMaskEdge
      allocate(i1Ptr)

! Setting up time level 1
      i1Ptr % fieldName = 'bdyMaskEdge'
      i1Ptr % isVarArray = .false.
      i1Ptr % isDecomposed = .true.
      i1Ptr % hasTimeDimension = .false.
      i1Ptr % isPersistent = .true.
      i1Ptr % isActive = .false.
! Setting up dimensions
      i1Ptr % dimNames(1) = 'nEdges'
      i1Ptr % defaultValue = 0
      allocate(i1Ptr % attLists(1))
      allocate(i1Ptr % attLists(1) % attList)
      call mpas_add_att(i1Ptr % attLists(1) % attList, 'units', '-')
      call mpas_add_att(i1Ptr % attLists(1) % attList, 'long_name', 'Indicator of whether an edge is an interior edge, a '&
//'relaxation-zone edge, or a specified-zone edge')
      i1Ptr % missingValue = MPAS_INT_FILLVAL
      i1Ptr % block => block

      i1Ptr % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'bdyMaskEdge', i1Ptr)
      call mpas_pool_add_field(block % allFields, 'bdyMaskEdge', i1Ptr)

! Define variable bdyMaskVertex
      allocate(i1Ptr)

! Setting up time level 1
      i1Ptr % fieldName = 'bdyMaskVertex'
      i1Ptr % isVarArray = .false.
      i1Ptr % isDecomposed = .true.
      i1Ptr % hasTimeDimension = .false.
      i1Ptr % isPersistent = .true.
      i1Ptr % isActive = .false.
! Setting up dimensions
      i1Ptr % dimNames(1) = 'nVertices'
      i1Ptr % defaultValue = 0
      allocate(i1Ptr % attLists(1))
      allocate(i1Ptr % attLists(1) % attList)
      call mpas_add_att(i1Ptr % attLists(1) % attList, 'units', '-')
      call mpas_add_att(i1Ptr % attLists(1) % attList, 'long_name', 'Indicator of whether an vertex is an interior vertex, a '&
//'relaxation-zone vertex, or a specified-zone vertex')
      i1Ptr % missingValue = MPAS_INT_FILLVAL
      i1Ptr % block => block

      i1Ptr % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'bdyMaskVertex', i1Ptr)
      call mpas_pool_add_field(block % allFields, 'bdyMaskVertex', i1Ptr)



      if (associated(newSubPool)) then
         call mpas_pool_add_config(newSubPool, 'on_a_sphere', block % domain % on_a_sphere)
         call mpas_pool_add_config(newSubPool, 'sphere_radius', block % domain % sphere_radius)
         call mpas_pool_add_config(newSubPool, 'is_periodic', block % domain % is_periodic)
         call mpas_pool_add_config(newSubPool, 'x_period', block % domain % x_period)
         call mpas_pool_add_config(newSubPool, 'y_period', block % domain % y_period)
      end if

   end subroutine init_atm_generate_pool_mesh


   subroutine init_atm_generate_pool_state(block, structPool, dimensionPool, packagePool)
      use mpas_derived_types
      use mpas_pool_routines
      use mpas_io_units
      use mpas_io, only : MPAS_REAL_FILLVAL, MPAS_INT_FILLVAL, MPAS_CHAR_FILLVAL
      implicit none
      type (block_type), intent(inout), pointer :: block
      type (mpas_pool_type), intent(inout) :: structPool
      type (mpas_pool_type), intent(inout) :: dimensionPool
      type (mpas_pool_type), intent(in) :: packagePool

      type (field0DReal), pointer :: r0Ptr
      type (field1DReal), pointer :: r1Ptr
      type (field2DReal), pointer :: r2Ptr
      type (field3DReal), pointer :: r3Ptr
      type (field4DReal), pointer :: r4Ptr
      type (field5DReal), pointer :: r5Ptr
      type (field0DInteger), pointer :: i0Ptr
      type (field1DInteger), pointer :: i1Ptr
      type (field2DInteger), pointer :: i2Ptr
      type (field3DInteger), pointer :: i3Ptr
      type (field0DChar), pointer :: c0Ptr
      type (field1DChar), pointer :: c1Ptr
      type (field0DReal), dimension(:), pointer :: r0aPtr
      type (field1DReal), dimension(:), pointer :: r1aPtr
      type (field2DReal), dimension(:), pointer :: r2aPtr
      type (field3DReal), dimension(:), pointer :: r3aPtr
      type (field4DReal), dimension(:), pointer :: r4aPtr
      type (field5DReal), dimension(:), pointer :: r5aPtr
      type (field0DInteger), dimension(:), pointer :: i0aPtr
      type (field1DInteger), dimension(:), pointer :: i1aPtr
      type (field2DInteger), dimension(:), pointer :: i2aPtr
      type (field3DInteger), dimension(:), pointer :: i3aPtr
      type (field0DChar), dimension(:), pointer :: c0aPtr
      type (field1DChar), dimension(:), pointer :: c1aPtr

      type (mpas_pool_type), pointer :: newSubPool
      integer :: group_counter
      logical :: group_started
      integer :: group_start
      integer :: index_counter
      integer, pointer :: const_index

      logical, pointer :: initial_condsActive
      logical, pointer :: sfc_updateActive
      logical, pointer :: lbcsActive
      logical, pointer :: gwd_stage_inActive
      logical, pointer :: gwd_stage_outActive
      logical, pointer :: vertical_stage_inActive
      logical, pointer :: vertical_stage_outActive
      logical, pointer :: met_stage_inActive
      logical, pointer :: met_stage_outActive
      logical, pointer :: first_guess_fieldActive


      integer :: numConstituents

      nullify(newSubPool)
      group_counter = -1
      group_started = .false.
      group_start = -1
      call mpas_pool_get_package(packagePool, 'initial_condsActive', initial_condsActive)
      call mpas_pool_get_package(packagePool, 'sfc_updateActive', sfc_updateActive)
      call mpas_pool_get_package(packagePool, 'lbcsActive', lbcsActive)
      call mpas_pool_get_package(packagePool, 'gwd_stage_inActive', gwd_stage_inActive)
      call mpas_pool_get_package(packagePool, 'gwd_stage_outActive', gwd_stage_outActive)
      call mpas_pool_get_package(packagePool, 'vertical_stage_inActive', vertical_stage_inActive)
      call mpas_pool_get_package(packagePool, 'vertical_stage_outActive', vertical_stage_outActive)
      call mpas_pool_get_package(packagePool, 'met_stage_inActive', met_stage_inActive)
      call mpas_pool_get_package(packagePool, 'met_stage_outActive', met_stage_outActive)
      call mpas_pool_get_package(packagePool, 'first_guess_fieldActive', first_guess_fieldActive)

      allocate(newSubPool)
      call mpas_pool_create_pool(newSubPool)
      call mpas_pool_add_subpool(structPool, 'state', newSubPool)
      call mpas_pool_add_subpool(block % allStructs, 'state', newSubPool)

! Define var array scalars
      allocate(r3Ptr)
      index_counter = 0
      group_counter = -1
      group_start = -1
      group_started = .false.

! Starting group moist
! Define constituent var qv
! My Packages are met_stage_out
      if (met_stage_outActive) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_qv', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'moist_start', group_start)
            end if
            group_started = .true.
         end if
      else
           if (associated(newSubPool)) then
              call mpas_pool_add_dimension(newSubPool, 'index_qv', -1)
           end if
      end if
! Define constituent var qc
! My packages are met_stage_out
      if (met_stage_outActive) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_qc', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'moist_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_qc', -1)
         end if
      end if
! Define constituent var qr
! My packages are met_stage_out
      if (met_stage_outActive) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_qr', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'moist_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_qr', -1)
         end if
      end if
         if (.not. group_started) then
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'moist_start', -1)
               call mpas_pool_add_dimension(newSubPool, 'moist_end', -1)
            end if
         else
            group_started = .false.
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'moist_end', index_counter)
            end if
         end if
! End of group       

      numConstituents = index_counter
      if (associated(newSubPool)) then
         call mpas_pool_add_dimension(newSubPool, 'num_scalars', numConstituents)
      end if
! Defining time level 1
      allocate( r3Ptr % constituentNames(numConstituents) )
      r3Ptr % fieldName = 'scalars'
      r3Ptr % isDecomposed = .true.
      r3Ptr % hasTimeDimension = .true.
      r3Ptr % isVarArray = .true.
      r3Ptr % isPersistent = .true.
      r3Ptr % isActive = .false.

      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_qv', const_index)
      end if
      if (const_index > 0) then
         r3Ptr % constituentNames(const_index) = 'qv'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_qc', const_index)
      end if
      if (const_index > 0) then
         r3Ptr % constituentNames(const_index) = 'qc'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_qr', const_index)
      end if
      if (const_index > 0) then
         r3Ptr % constituentNames(const_index) = 'qr'
      end if

! Setup dimensions for       
      r3Ptr % dimNames(1) = 'num_scalars'
      r3Ptr % dimNames(2) = 'nVertLevels'
      r3Ptr % dimNames(3) = 'nCells'

      r3Ptr % defaultValue = 0.0
      allocate(r3Ptr % attLists(size(r3Ptr % constituentNames, dim=1)))
      do index_counter = 1, size(r3Ptr % constituentNames, dim=1)
         allocate(r3Ptr % attLists(index_counter) % attList)
      end do
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_qv', const_index)
      end if
      if (const_index > 0) then
         call mpas_add_att(r3Ptr % attLists(const_index) % attList, 'long_name', 'Water vapor mixing ratio')
         call mpas_add_att(r3Ptr % attLists(const_index) % attList, 'units', 'kg kg^{-1}')
         r3Ptr % missingValue = MPAS_REAL_FILLVAL
         r3Ptr % constituentNames(const_index) = 'qv'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_qc', const_index)
      end if
      if (const_index > 0) then
         call mpas_add_att(r3Ptr % attLists(const_index) % attList, 'long_name', 'Cloud water mixing ratio')
         call mpas_add_att(r3Ptr % attLists(const_index) % attList, 'units', 'kg kg^{-1}')
         r3Ptr % missingValue = MPAS_REAL_FILLVAL
         r3Ptr % constituentNames(const_index) = 'qc'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_qr', const_index)
      end if
      if (const_index > 0) then
         call mpas_add_att(r3Ptr % attLists(const_index) % attList, 'long_name', 'Rain water mixing ratio')
         call mpas_add_att(r3Ptr % attLists(const_index) % attList, 'units', 'kg kg^{-1}')
         r3Ptr % missingValue = MPAS_REAL_FILLVAL
         r3Ptr % constituentNames(const_index) = 'qr'
      end if
      r3Ptr % block => block

      if (met_stage_outActive) then
         r3Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'scalars', r3Ptr)
      call mpas_pool_add_field(block % allFields, 'scalars', r3Ptr)

! Define variable initial_time
      allocate(c0Ptr)

! Setting up time level 1
      c0Ptr % fieldName = 'initial_time'
      c0Ptr % isVarArray = .false.
      c0Ptr % isDecomposed = .false.
      c0Ptr % hasTimeDimension = .false.
      c0Ptr % defaultValue = ''
      c0Ptr % scalar = ''
      allocate(c0Ptr % attLists(1))
      allocate(c0Ptr % attLists(1) % attList)
      call mpas_add_att(c0Ptr % attLists(1) % attList, 'units', 'YYYY-MM-DD_hh:mm:ss')
      call mpas_add_att(c0Ptr % attLists(1) % attList, 'long_name', 'Model initialization time')
      c0Ptr % missingValue = MPAS_CHAR_FILLVAL
      c0Ptr % block => block

      c0Ptr % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'initial_time', c0Ptr)
      call mpas_pool_add_field(block % allFields, 'initial_time', c0Ptr)

! Define variable xtime
      allocate(c0Ptr)

! Setting up time level 1
      c0Ptr % fieldName = 'xtime'
      c0Ptr % isVarArray = .false.
      c0Ptr % isDecomposed = .false.
      c0Ptr % hasTimeDimension = .true.
      c0Ptr % defaultValue = ''
      c0Ptr % scalar = ''
      allocate(c0Ptr % attLists(1))
      allocate(c0Ptr % attLists(1) % attList)
      call mpas_add_att(c0Ptr % attLists(1) % attList, 'units', 'YYYY-MM-DD_hh:mm:ss')
      call mpas_add_att(c0Ptr % attLists(1) % attList, 'long_name', 'Model valid time')
      c0Ptr % missingValue = MPAS_CHAR_FILLVAL
      c0Ptr % block => block

      c0Ptr % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'xtime', c0Ptr)
      call mpas_pool_add_field(block % allFields, 'xtime', c0Ptr)

! Define variable u
      allocate(r2Ptr)

! Setting up time level 1
      r2Ptr % fieldName = 'u'
      r2Ptr % isVarArray = .false.
      r2Ptr % isDecomposed = .true.
      r2Ptr % hasTimeDimension = .true.
      r2Ptr % isPersistent = .true.
      r2Ptr % isActive = .false.
! Setting up dimensions
      r2Ptr % dimNames(1) = 'nVertLevels'
      r2Ptr % dimNames(2) = 'nEdges'
      r2Ptr % defaultValue = 0.0
      allocate(r2Ptr % attLists(1))
      allocate(r2Ptr % attLists(1) % attList)
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'units', 'm s^{-1}')
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'long_name', 'Horizontal normal velocity at edges')
      r2Ptr % missingValue = MPAS_REAL_FILLVAL
      r2Ptr % block => block

      if (met_stage_outActive) then
         r2Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'u', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'u', r2Ptr)

! Define variable w
      allocate(r2Ptr)

! Setting up time level 1
      r2Ptr % fieldName = 'w'
      r2Ptr % isVarArray = .false.
      r2Ptr % isDecomposed = .true.
      r2Ptr % hasTimeDimension = .true.
      r2Ptr % isPersistent = .true.
      r2Ptr % isActive = .false.
! Setting up dimensions
      r2Ptr % dimNames(1) = 'nVertLevelsP1'
      r2Ptr % dimNames(2) = 'nCells'
      r2Ptr % defaultValue = 0.0
      allocate(r2Ptr % attLists(1))
      allocate(r2Ptr % attLists(1) % attList)
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'units', 'm s^{-1}')
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'long_name', 'Vertical velocity at vertical cell faces')
      r2Ptr % missingValue = MPAS_REAL_FILLVAL
      r2Ptr % block => block

      if (met_stage_outActive) then
         r2Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'w', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'w', r2Ptr)

! Define variable rho_zz
      allocate(r2Ptr)

! Setting up time level 1
      r2Ptr % fieldName = 'rho_zz'
      r2Ptr % isVarArray = .false.
      r2Ptr % isDecomposed = .true.
      r2Ptr % hasTimeDimension = .true.
      r2Ptr % isPersistent = .true.
      r2Ptr % isActive = .false.
! Setting up dimensions
      r2Ptr % dimNames(1) = 'nVertLevels'
      r2Ptr % dimNames(2) = 'nCells'
      r2Ptr % defaultValue = 0.0
      allocate(r2Ptr % attLists(1))
      allocate(r2Ptr % attLists(1) % attList)
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'units', 'kg m^{-3}')
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'long_name', 'Dry air density divided by d(zeta)/dz')
      r2Ptr % missingValue = MPAS_REAL_FILLVAL
      r2Ptr % block => block

      if (met_stage_outActive) then
         r2Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'rho_zz', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'rho_zz', r2Ptr)

! Define variable theta_m
      allocate(r2Ptr)

! Setting up time level 1
      r2Ptr % fieldName = 'theta_m'
      r2Ptr % isVarArray = .false.
      r2Ptr % isDecomposed = .true.
      r2Ptr % hasTimeDimension = .true.
      r2Ptr % isPersistent = .true.
      r2Ptr % isActive = .false.
! Setting up dimensions
      r2Ptr % dimNames(1) = 'nVertLevels'
      r2Ptr % dimNames(2) = 'nCells'
      r2Ptr % defaultValue = 0.0
      allocate(r2Ptr % attLists(1))
      allocate(r2Ptr % attLists(1) % attList)
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'units', 'K')
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'long_name', 'Moist potential temperature: theta*(1+q_v*R_v/R_d)')
      r2Ptr % missingValue = MPAS_REAL_FILLVAL
      r2Ptr % block => block

      if (met_stage_outActive) then
         r2Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'theta_m', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'theta_m', r2Ptr)

! Define variable h_oml_initial
      allocate(r1Ptr)

! Setting up time level 1
      r1Ptr % fieldName = 'h_oml_initial'
      r1Ptr % isVarArray = .false.
      r1Ptr % isDecomposed = .true.
      r1Ptr % hasTimeDimension = .true.
      r1Ptr % isPersistent = .true.
      r1Ptr % isActive = .false.
! Setting up dimensions
      r1Ptr % dimNames(1) = 'nCells'
      r1Ptr % defaultValue = 0.0
      allocate(r1Ptr % attLists(1))
      allocate(r1Ptr % attLists(1) % attList)
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'units', 'm')
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'long_name', 'Initial depth of ocean mix layer')
      r1Ptr % missingValue = MPAS_REAL_FILLVAL
      r1Ptr % block => block

      if (met_stage_outActive) then
         r1Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'h_oml_initial', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'h_oml_initial', r1Ptr)



      if (associated(newSubPool)) then
         call mpas_pool_add_config(newSubPool, 'on_a_sphere', block % domain % on_a_sphere)
         call mpas_pool_add_config(newSubPool, 'sphere_radius', block % domain % sphere_radius)
         call mpas_pool_add_config(newSubPool, 'is_periodic', block % domain % is_periodic)
         call mpas_pool_add_config(newSubPool, 'x_period', block % domain % x_period)
         call mpas_pool_add_config(newSubPool, 'y_period', block % domain % y_period)
      end if

   end subroutine init_atm_generate_pool_state


   subroutine init_atm_generate_pool_lbc_state(block, structPool, dimensionPool, packagePool)
      use mpas_derived_types
      use mpas_pool_routines
      use mpas_io_units
      use mpas_io, only : MPAS_REAL_FILLVAL, MPAS_INT_FILLVAL, MPAS_CHAR_FILLVAL
      implicit none
      type (block_type), intent(inout), pointer :: block
      type (mpas_pool_type), intent(inout) :: structPool
      type (mpas_pool_type), intent(inout) :: dimensionPool
      type (mpas_pool_type), intent(in) :: packagePool

      type (field0DReal), pointer :: r0Ptr
      type (field1DReal), pointer :: r1Ptr
      type (field2DReal), pointer :: r2Ptr
      type (field3DReal), pointer :: r3Ptr
      type (field4DReal), pointer :: r4Ptr
      type (field5DReal), pointer :: r5Ptr
      type (field0DInteger), pointer :: i0Ptr
      type (field1DInteger), pointer :: i1Ptr
      type (field2DInteger), pointer :: i2Ptr
      type (field3DInteger), pointer :: i3Ptr
      type (field0DChar), pointer :: c0Ptr
      type (field1DChar), pointer :: c1Ptr
      type (field0DReal), dimension(:), pointer :: r0aPtr
      type (field1DReal), dimension(:), pointer :: r1aPtr
      type (field2DReal), dimension(:), pointer :: r2aPtr
      type (field3DReal), dimension(:), pointer :: r3aPtr
      type (field4DReal), dimension(:), pointer :: r4aPtr
      type (field5DReal), dimension(:), pointer :: r5aPtr
      type (field0DInteger), dimension(:), pointer :: i0aPtr
      type (field1DInteger), dimension(:), pointer :: i1aPtr
      type (field2DInteger), dimension(:), pointer :: i2aPtr
      type (field3DInteger), dimension(:), pointer :: i3aPtr
      type (field0DChar), dimension(:), pointer :: c0aPtr
      type (field1DChar), dimension(:), pointer :: c1aPtr

      type (mpas_pool_type), pointer :: newSubPool
      integer :: group_counter
      logical :: group_started
      integer :: group_start
      integer :: index_counter
      integer, pointer :: const_index

      logical, pointer :: initial_condsActive
      logical, pointer :: sfc_updateActive
      logical, pointer :: lbcsActive
      logical, pointer :: gwd_stage_inActive
      logical, pointer :: gwd_stage_outActive
      logical, pointer :: vertical_stage_inActive
      logical, pointer :: vertical_stage_outActive
      logical, pointer :: met_stage_inActive
      logical, pointer :: met_stage_outActive
      logical, pointer :: first_guess_fieldActive


      integer :: numConstituents

      nullify(newSubPool)
      group_counter = -1
      group_started = .false.
      group_start = -1
      call mpas_pool_get_package(packagePool, 'initial_condsActive', initial_condsActive)
      call mpas_pool_get_package(packagePool, 'sfc_updateActive', sfc_updateActive)
      call mpas_pool_get_package(packagePool, 'lbcsActive', lbcsActive)
      call mpas_pool_get_package(packagePool, 'gwd_stage_inActive', gwd_stage_inActive)
      call mpas_pool_get_package(packagePool, 'gwd_stage_outActive', gwd_stage_outActive)
      call mpas_pool_get_package(packagePool, 'vertical_stage_inActive', vertical_stage_inActive)
      call mpas_pool_get_package(packagePool, 'vertical_stage_outActive', vertical_stage_outActive)
      call mpas_pool_get_package(packagePool, 'met_stage_inActive', met_stage_inActive)
      call mpas_pool_get_package(packagePool, 'met_stage_outActive', met_stage_outActive)
      call mpas_pool_get_package(packagePool, 'first_guess_fieldActive', first_guess_fieldActive)

      allocate(newSubPool)
      call mpas_pool_create_pool(newSubPool)
      call mpas_pool_add_subpool(structPool, 'lbc_state', newSubPool)
      call mpas_pool_add_subpool(block % allStructs, 'lbc_state', newSubPool)

! Define var array lbc_scalars
      allocate(r3Ptr)
      index_counter = 0
      group_counter = -1
      group_start = -1
      group_started = .false.

! Starting group moist
! Define constituent var lbc_qv
! My Packages are lbcs
      if (lbcsActive) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_lbc_qv', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'moist_start', group_start)
            end if
            group_started = .true.
         end if
      else
           if (associated(newSubPool)) then
              call mpas_pool_add_dimension(newSubPool, 'index_lbc_qv', -1)
           end if
      end if
! Define constituent var lbc_qc
! My packages are lbcs
      if (lbcsActive) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_lbc_qc', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'moist_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_lbc_qc', -1)
         end if
      end if
! Define constituent var lbc_qr
! My packages are lbcs
      if (lbcsActive) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_lbc_qr', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'moist_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_lbc_qr', -1)
         end if
      end if
         if (.not. group_started) then
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'moist_start', -1)
               call mpas_pool_add_dimension(newSubPool, 'moist_end', -1)
            end if
         else
            group_started = .false.
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'moist_end', index_counter)
            end if
         end if
! End of group       

      numConstituents = index_counter
      if (associated(newSubPool)) then
         call mpas_pool_add_dimension(newSubPool, 'num_lbc_scalars', numConstituents)
      end if
! Defining time level 1
      allocate( r3Ptr % constituentNames(numConstituents) )
      r3Ptr % fieldName = 'lbc_scalars'
      r3Ptr % isDecomposed = .true.
      r3Ptr % hasTimeDimension = .true.
      r3Ptr % isVarArray = .true.
      r3Ptr % isPersistent = .true.
      r3Ptr % isActive = .false.

      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_lbc_qv', const_index)
      end if
      if (const_index > 0) then
         r3Ptr % constituentNames(const_index) = 'lbc_qv'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_lbc_qc', const_index)
      end if
      if (const_index > 0) then
         r3Ptr % constituentNames(const_index) = 'lbc_qc'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_lbc_qr', const_index)
      end if
      if (const_index > 0) then
         r3Ptr % constituentNames(const_index) = 'lbc_qr'
      end if

! Setup dimensions for       
      r3Ptr % dimNames(1) = 'num_lbc_scalars'
      r3Ptr % dimNames(2) = 'nVertLevels'
      r3Ptr % dimNames(3) = 'nCells'

      r3Ptr % defaultValue = 0.0
      allocate(r3Ptr % attLists(size(r3Ptr % constituentNames, dim=1)))
      do index_counter = 1, size(r3Ptr % constituentNames, dim=1)
         allocate(r3Ptr % attLists(index_counter) % attList)
      end do
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_lbc_qv', const_index)
      end if
      if (const_index > 0) then
         call mpas_add_att(r3Ptr % attLists(const_index) % attList, 'long_name', 'Water vapor mixing ratio on lateral boundary '&
//'cells')
         call mpas_add_att(r3Ptr % attLists(const_index) % attList, 'units', 'kg kg^{-1}')
         r3Ptr % missingValue = MPAS_REAL_FILLVAL
         r3Ptr % constituentNames(const_index) = 'lbc_qv'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_lbc_qc', const_index)
      end if
      if (const_index > 0) then
         call mpas_add_att(r3Ptr % attLists(const_index) % attList, 'long_name', 'Cloud water mixing ratio on lateral boundary '&
//'cells')
         call mpas_add_att(r3Ptr % attLists(const_index) % attList, 'units', 'kg kg^{-1}')
         r3Ptr % missingValue = MPAS_REAL_FILLVAL
         r3Ptr % constituentNames(const_index) = 'lbc_qc'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_lbc_qr', const_index)
      end if
      if (const_index > 0) then
         call mpas_add_att(r3Ptr % attLists(const_index) % attList, 'long_name', 'Rain water mixing ratio on lateral boundary '&
//'cells')
         call mpas_add_att(r3Ptr % attLists(const_index) % attList, 'units', 'kg kg^{-1}')
         r3Ptr % missingValue = MPAS_REAL_FILLVAL
         r3Ptr % constituentNames(const_index) = 'lbc_qr'
      end if
      r3Ptr % block => block

      if (lbcsActive) then
         r3Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'lbc_scalars', r3Ptr)
      call mpas_pool_add_field(block % allFields, 'lbc_scalars', r3Ptr)

! Define variable lbc_u
      allocate(r2Ptr)

! Setting up time level 1
      r2Ptr % fieldName = 'lbc_u'
      r2Ptr % isVarArray = .false.
      r2Ptr % isDecomposed = .true.
      r2Ptr % hasTimeDimension = .true.
      r2Ptr % isPersistent = .true.
      r2Ptr % isActive = .false.
! Setting up dimensions
      r2Ptr % dimNames(1) = 'nVertLevels'
      r2Ptr % dimNames(2) = 'nEdges'
      r2Ptr % defaultValue = 0.0
      allocate(r2Ptr % attLists(1))
      allocate(r2Ptr % attLists(1) % attList)
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'units', 'm s^{-1}')
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'long_name', 'Horizontal normal velocity on domain lateral boundary edges')
      r2Ptr % missingValue = MPAS_REAL_FILLVAL
      r2Ptr % block => block

      if (lbcsActive) then
         r2Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'lbc_u', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'lbc_u', r2Ptr)

! Define variable lbc_w
      allocate(r2Ptr)

! Setting up time level 1
      r2Ptr % fieldName = 'lbc_w'
      r2Ptr % isVarArray = .false.
      r2Ptr % isDecomposed = .true.
      r2Ptr % hasTimeDimension = .true.
      r2Ptr % isPersistent = .true.
      r2Ptr % isActive = .false.
! Setting up dimensions
      r2Ptr % dimNames(1) = 'nVertLevelsP1'
      r2Ptr % dimNames(2) = 'nCells'
      r2Ptr % defaultValue = 0.0
      allocate(r2Ptr % attLists(1))
      allocate(r2Ptr % attLists(1) % attList)
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'units', 'm s^{-1}')
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'long_name', 'Vertical velocity on domain lateral boundary vertical cell '&
//'faces')
      r2Ptr % missingValue = MPAS_REAL_FILLVAL
      r2Ptr % block => block

      if (lbcsActive) then
         r2Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'lbc_w', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'lbc_w', r2Ptr)

! Define variable lbc_rho
      allocate(r2Ptr)

! Setting up time level 1
      r2Ptr % fieldName = 'lbc_rho'
      r2Ptr % isVarArray = .false.
      r2Ptr % isDecomposed = .true.
      r2Ptr % hasTimeDimension = .true.
      r2Ptr % isPersistent = .true.
      r2Ptr % isActive = .false.
! Setting up dimensions
      r2Ptr % dimNames(1) = 'nVertLevels'
      r2Ptr % dimNames(2) = 'nCells'
      r2Ptr % defaultValue = 0.0
      allocate(r2Ptr % attLists(1))
      allocate(r2Ptr % attLists(1) % attList)
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'units', 'kg m^{-3}')
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'long_name', 'Dry air density on lateral boundary cells')
      r2Ptr % missingValue = MPAS_REAL_FILLVAL
      r2Ptr % block => block

      if (lbcsActive) then
         r2Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'lbc_rho', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'lbc_rho', r2Ptr)

! Define variable lbc_theta
      allocate(r2Ptr)

! Setting up time level 1
      r2Ptr % fieldName = 'lbc_theta'
      r2Ptr % isVarArray = .false.
      r2Ptr % isDecomposed = .true.
      r2Ptr % hasTimeDimension = .true.
      r2Ptr % isPersistent = .true.
      r2Ptr % isActive = .false.
! Setting up dimensions
      r2Ptr % dimNames(1) = 'nVertLevels'
      r2Ptr % dimNames(2) = 'nCells'
      r2Ptr % defaultValue = 0.0
      allocate(r2Ptr % attLists(1))
      allocate(r2Ptr % attLists(1) % attList)
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'units', 'K')
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'long_name', 'Potential temperature on lateral boundary cells')
      r2Ptr % missingValue = MPAS_REAL_FILLVAL
      r2Ptr % block => block

      if (lbcsActive) then
         r2Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'lbc_theta', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'lbc_theta', r2Ptr)



      if (associated(newSubPool)) then
         call mpas_pool_add_config(newSubPool, 'on_a_sphere', block % domain % on_a_sphere)
         call mpas_pool_add_config(newSubPool, 'sphere_radius', block % domain % sphere_radius)
         call mpas_pool_add_config(newSubPool, 'is_periodic', block % domain % is_periodic)
         call mpas_pool_add_config(newSubPool, 'x_period', block % domain % x_period)
         call mpas_pool_add_config(newSubPool, 'y_period', block % domain % y_period)
      end if

   end subroutine init_atm_generate_pool_lbc_state


   subroutine init_atm_generate_pool_fg(block, structPool, dimensionPool, packagePool)
      use mpas_derived_types
      use mpas_pool_routines
      use mpas_io_units
      use mpas_io, only : MPAS_REAL_FILLVAL, MPAS_INT_FILLVAL, MPAS_CHAR_FILLVAL
      implicit none
      type (block_type), intent(inout), pointer :: block
      type (mpas_pool_type), intent(inout) :: structPool
      type (mpas_pool_type), intent(inout) :: dimensionPool
      type (mpas_pool_type), intent(in) :: packagePool

      type (field0DReal), pointer :: r0Ptr
      type (field1DReal), pointer :: r1Ptr
      type (field2DReal), pointer :: r2Ptr
      type (field3DReal), pointer :: r3Ptr
      type (field4DReal), pointer :: r4Ptr
      type (field5DReal), pointer :: r5Ptr
      type (field0DInteger), pointer :: i0Ptr
      type (field1DInteger), pointer :: i1Ptr
      type (field2DInteger), pointer :: i2Ptr
      type (field3DInteger), pointer :: i3Ptr
      type (field0DChar), pointer :: c0Ptr
      type (field1DChar), pointer :: c1Ptr
      type (field0DReal), dimension(:), pointer :: r0aPtr
      type (field1DReal), dimension(:), pointer :: r1aPtr
      type (field2DReal), dimension(:), pointer :: r2aPtr
      type (field3DReal), dimension(:), pointer :: r3aPtr
      type (field4DReal), dimension(:), pointer :: r4aPtr
      type (field5DReal), dimension(:), pointer :: r5aPtr
      type (field0DInteger), dimension(:), pointer :: i0aPtr
      type (field1DInteger), dimension(:), pointer :: i1aPtr
      type (field2DInteger), dimension(:), pointer :: i2aPtr
      type (field3DInteger), dimension(:), pointer :: i3aPtr
      type (field0DChar), dimension(:), pointer :: c0aPtr
      type (field1DChar), dimension(:), pointer :: c1aPtr

      type (mpas_pool_type), pointer :: newSubPool
      integer :: group_counter
      logical :: group_started
      integer :: group_start
      integer :: index_counter
      integer, pointer :: const_index

      logical, pointer :: initial_condsActive
      logical, pointer :: sfc_updateActive
      logical, pointer :: lbcsActive
      logical, pointer :: gwd_stage_inActive
      logical, pointer :: gwd_stage_outActive
      logical, pointer :: vertical_stage_inActive
      logical, pointer :: vertical_stage_outActive
      logical, pointer :: met_stage_inActive
      logical, pointer :: met_stage_outActive
      logical, pointer :: first_guess_fieldActive


      integer :: numConstituents

      nullify(newSubPool)
      group_counter = -1
      group_started = .false.
      group_start = -1
      call mpas_pool_get_package(packagePool, 'initial_condsActive', initial_condsActive)
      call mpas_pool_get_package(packagePool, 'sfc_updateActive', sfc_updateActive)
      call mpas_pool_get_package(packagePool, 'lbcsActive', lbcsActive)
      call mpas_pool_get_package(packagePool, 'gwd_stage_inActive', gwd_stage_inActive)
      call mpas_pool_get_package(packagePool, 'gwd_stage_outActive', gwd_stage_outActive)
      call mpas_pool_get_package(packagePool, 'vertical_stage_inActive', vertical_stage_inActive)
      call mpas_pool_get_package(packagePool, 'vertical_stage_outActive', vertical_stage_outActive)
      call mpas_pool_get_package(packagePool, 'met_stage_inActive', met_stage_inActive)
      call mpas_pool_get_package(packagePool, 'met_stage_outActive', met_stage_outActive)
      call mpas_pool_get_package(packagePool, 'first_guess_fieldActive', first_guess_fieldActive)

      allocate(newSubPool)
      call mpas_pool_create_pool(newSubPool)
      call mpas_pool_add_subpool(structPool, 'fg', newSubPool)
      call mpas_pool_add_subpool(block % allStructs, 'fg', newSubPool)

! Define variable u_fg
      allocate(r2Ptr)

! Setting up time level 1
      r2Ptr % fieldName = 'u_fg'
      r2Ptr % isVarArray = .false.
      r2Ptr % isDecomposed = .true.
      r2Ptr % hasTimeDimension = .true.
      r2Ptr % isPersistent = .true.
      r2Ptr % isActive = .false.
! Setting up dimensions
      r2Ptr % dimNames(1) = 'nFGLevels'
      r2Ptr % dimNames(2) = 'nEdges'
      r2Ptr % defaultValue = 0.0
      allocate(r2Ptr % attLists(1))
      allocate(r2Ptr % attLists(1) % attList)
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'units', 'm s^{-1}')
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'long_name', 'First guess zonal wind component')
      r2Ptr % missingValue = MPAS_REAL_FILLVAL
      r2Ptr % block => block

      if (first_guess_fieldActive) then
         r2Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'u', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'u_fg', r2Ptr)

! Define variable v_fg
      allocate(r2Ptr)

! Setting up time level 1
      r2Ptr % fieldName = 'v_fg'
      r2Ptr % isVarArray = .false.
      r2Ptr % isDecomposed = .true.
      r2Ptr % hasTimeDimension = .true.
      r2Ptr % isPersistent = .true.
      r2Ptr % isActive = .false.
! Setting up dimensions
      r2Ptr % dimNames(1) = 'nFGLevels'
      r2Ptr % dimNames(2) = 'nEdges'
      r2Ptr % defaultValue = 0.0
      allocate(r2Ptr % attLists(1))
      allocate(r2Ptr % attLists(1) % attList)
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'units', 'm s^{-1}')
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'long_name', 'First guess merdian wind component')
      r2Ptr % missingValue = MPAS_REAL_FILLVAL
      r2Ptr % block => block

      if (first_guess_fieldActive) then
         r2Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'v', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'v_fg', r2Ptr)

! Define variable t_fg
      allocate(r2Ptr)

! Setting up time level 1
      r2Ptr % fieldName = 't_fg'
      r2Ptr % isVarArray = .false.
      r2Ptr % isDecomposed = .true.
      r2Ptr % hasTimeDimension = .true.
      r2Ptr % isPersistent = .true.
      r2Ptr % isActive = .false.
! Setting up dimensions
      r2Ptr % dimNames(1) = 'nFGLevels'
      r2Ptr % dimNames(2) = 'nCells'
      r2Ptr % defaultValue = 0.0
      allocate(r2Ptr % attLists(1))
      allocate(r2Ptr % attLists(1) % attList)
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'units', 'K')
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'long_name', 'First guess temperature')
      r2Ptr % missingValue = MPAS_REAL_FILLVAL
      r2Ptr % block => block

      if (first_guess_fieldActive) then
         r2Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 't', r2Ptr)
      call mpas_pool_add_field(block % allFields, 't_fg', r2Ptr)

! Define variable p_fg
      allocate(r2Ptr)

! Setting up time level 1
      r2Ptr % fieldName = 'p_fg'
      r2Ptr % isVarArray = .false.
      r2Ptr % isDecomposed = .true.
      r2Ptr % hasTimeDimension = .true.
      r2Ptr % isPersistent = .true.
      r2Ptr % isActive = .false.
! Setting up dimensions
      r2Ptr % dimNames(1) = 'nFGLevels'
      r2Ptr % dimNames(2) = 'nCells'
      r2Ptr % defaultValue = 0.0
      allocate(r2Ptr % attLists(1))
      allocate(r2Ptr % attLists(1) % attList)
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'units', 'Pa')
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'long_name', 'First guess pressure')
      r2Ptr % missingValue = MPAS_REAL_FILLVAL
      r2Ptr % block => block

      if (first_guess_fieldActive) then
         r2Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'p', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'p_fg', r2Ptr)

! Define variable z_fg
      allocate(r2Ptr)

! Setting up time level 1
      r2Ptr % fieldName = 'z_fg'
      r2Ptr % isVarArray = .false.
      r2Ptr % isDecomposed = .true.
      r2Ptr % hasTimeDimension = .true.
      r2Ptr % isPersistent = .true.
      r2Ptr % isActive = .false.
! Setting up dimensions
      r2Ptr % dimNames(1) = 'nFGLevels'
      r2Ptr % dimNames(2) = 'nCells'
      r2Ptr % defaultValue = 0.0
      allocate(r2Ptr % attLists(1))
      allocate(r2Ptr % attLists(1) % attList)
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'units', 'm')
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'long_name', 'First guess geopotential height')
      r2Ptr % missingValue = MPAS_REAL_FILLVAL
      r2Ptr % block => block

      if (first_guess_fieldActive) then
         r2Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'z', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'z_fg', r2Ptr)

! Define variable rh_fg
      allocate(r2Ptr)

! Setting up time level 1
      r2Ptr % fieldName = 'rh_fg'
      r2Ptr % isVarArray = .false.
      r2Ptr % isDecomposed = .true.
      r2Ptr % hasTimeDimension = .true.
      r2Ptr % isPersistent = .true.
      r2Ptr % isActive = .false.
! Setting up dimensions
      r2Ptr % dimNames(1) = 'nFGLevels'
      r2Ptr % dimNames(2) = 'nCells'
      r2Ptr % defaultValue = 0.0
      allocate(r2Ptr % attLists(1))
      allocate(r2Ptr % attLists(1) % attList)
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'units', 'percent')
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'long_name', 'First guess relative humidity with respect to liquid water')
      r2Ptr % missingValue = MPAS_REAL_FILLVAL
      r2Ptr % block => block

      if (first_guess_fieldActive) then
         r2Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'rh', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'rh_fg', r2Ptr)

! Define variable sh_fg
      allocate(r2Ptr)

! Setting up time level 1
      r2Ptr % fieldName = 'sh_fg'
      r2Ptr % isVarArray = .false.
      r2Ptr % isDecomposed = .true.
      r2Ptr % hasTimeDimension = .true.
      r2Ptr % isPersistent = .true.
      r2Ptr % isActive = .false.
! Setting up dimensions
      r2Ptr % dimNames(1) = 'nFGLevels'
      r2Ptr % dimNames(2) = 'nCells'
      r2Ptr % defaultValue = 0.0
      allocate(r2Ptr % attLists(1))
      allocate(r2Ptr % attLists(1) % attList)
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'units', 'kg kg^{-1}')
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'long_name', 'First guess specific humidity')
      r2Ptr % missingValue = MPAS_REAL_FILLVAL
      r2Ptr % block => block

      if (first_guess_fieldActive) then
         r2Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'sh', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'sh_fg', r2Ptr)

! Define variable soilz_fg
      allocate(r1Ptr)

! Setting up time level 1
      r1Ptr % fieldName = 'soilz_fg'
      r1Ptr % isVarArray = .false.
      r1Ptr % isDecomposed = .true.
      r1Ptr % hasTimeDimension = .true.
      r1Ptr % isPersistent = .true.
      r1Ptr % isActive = .false.
! Setting up dimensions
      r1Ptr % dimNames(1) = 'nCells'
      r1Ptr % defaultValue = 0.0
      allocate(r1Ptr % attLists(1))
      allocate(r1Ptr % attLists(1) % attList)
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'units', 'm')
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'long_name', 'First guess soil height')
      r1Ptr % missingValue = MPAS_REAL_FILLVAL
      r1Ptr % block => block

      if (first_guess_fieldActive) then
         r1Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'soilz', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'soilz_fg', r1Ptr)

! Define variable psfc_fg
      allocate(r1Ptr)

! Setting up time level 1
      r1Ptr % fieldName = 'psfc_fg'
      r1Ptr % isVarArray = .false.
      r1Ptr % isDecomposed = .true.
      r1Ptr % hasTimeDimension = .true.
      r1Ptr % isPersistent = .true.
      r1Ptr % isActive = .false.
! Setting up dimensions
      r1Ptr % dimNames(1) = 'nCells'
      r1Ptr % defaultValue = 0.0
      allocate(r1Ptr % attLists(1))
      allocate(r1Ptr % attLists(1) % attList)
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'units', 'Pa')
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'long_name', 'First guess surface pressure')
      r1Ptr % missingValue = MPAS_REAL_FILLVAL
      r1Ptr % block => block

      if (first_guess_fieldActive) then
         r1Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'psfc', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'psfc_fg', r1Ptr)

! Define variable pmsl_fg
      allocate(r1Ptr)

! Setting up time level 1
      r1Ptr % fieldName = 'pmsl_fg'
      r1Ptr % isVarArray = .false.
      r1Ptr % isDecomposed = .true.
      r1Ptr % hasTimeDimension = .true.
      r1Ptr % isPersistent = .true.
      r1Ptr % isActive = .false.
! Setting up dimensions
      r1Ptr % dimNames(1) = 'nCells'
      r1Ptr % defaultValue = 0.0
      allocate(r1Ptr % attLists(1))
      allocate(r1Ptr % attLists(1) % attList)
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'units', 'Pa')
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'long_name', 'First guess mean sea level pressure')
      r1Ptr % missingValue = MPAS_REAL_FILLVAL
      r1Ptr % block => block

      if (first_guess_fieldActive) then
         r1Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'pmsl', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'pmsl_fg', r1Ptr)

! Define variable dz_fg
      allocate(r2Ptr)

! Setting up time level 1
      r2Ptr % fieldName = 'dz_fg'
      r2Ptr % isVarArray = .false.
      r2Ptr % isDecomposed = .true.
      r2Ptr % hasTimeDimension = .true.
      r2Ptr % isPersistent = .true.
      r2Ptr % isActive = .false.
! Setting up dimensions
      r2Ptr % dimNames(1) = 'nFGSoilLevels'
      r2Ptr % dimNames(2) = 'nCells'
      r2Ptr % defaultValue = 0.0
      allocate(r2Ptr % attLists(1))
      allocate(r2Ptr % attLists(1) % attList)
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'units', 'cm')
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'long_name', 'First guess depth of soil layer bottom')
      r2Ptr % missingValue = MPAS_REAL_FILLVAL
      r2Ptr % block => block

      if (first_guess_fieldActive) then
         r2Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'dz_fg', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'dz_fg', r2Ptr)

! Define variable dzs_fg
      allocate(r2Ptr)

! Setting up time level 1
      r2Ptr % fieldName = 'dzs_fg'
      r2Ptr % isVarArray = .false.
      r2Ptr % isDecomposed = .true.
      r2Ptr % hasTimeDimension = .true.
      r2Ptr % isPersistent = .true.
      r2Ptr % isActive = .false.
! Setting up dimensions
      r2Ptr % dimNames(1) = 'nFGSoilLevels'
      r2Ptr % dimNames(2) = 'nCells'
      r2Ptr % defaultValue = 0.0
      allocate(r2Ptr % attLists(1))
      allocate(r2Ptr % attLists(1) % attList)
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'units', 'cm')
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'long_name', 'First guess soil layer thickness')
      r2Ptr % missingValue = MPAS_REAL_FILLVAL
      r2Ptr % block => block

      if (first_guess_fieldActive) then
         r2Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'dzs_fg', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'dzs_fg', r2Ptr)

! Define variable zs_fg
      allocate(r2Ptr)

! Setting up time level 1
      r2Ptr % fieldName = 'zs_fg'
      r2Ptr % isVarArray = .false.
      r2Ptr % isDecomposed = .true.
      r2Ptr % hasTimeDimension = .true.
      r2Ptr % isPersistent = .true.
      r2Ptr % isActive = .false.
! Setting up dimensions
      r2Ptr % dimNames(1) = 'nFGSoilLevels'
      r2Ptr % dimNames(2) = 'nCells'
      r2Ptr % defaultValue = 0.0
      allocate(r2Ptr % attLists(1))
      allocate(r2Ptr % attLists(1) % attList)
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'units', 'cm')
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'long_name', 'First guess depth of centers of soil levels')
      r2Ptr % missingValue = MPAS_REAL_FILLVAL
      r2Ptr % block => block

      if (first_guess_fieldActive) then
         r2Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'zs_fg', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'zs_fg', r2Ptr)

! Define variable st_fg
      allocate(r2Ptr)

! Setting up time level 1
      r2Ptr % fieldName = 'st_fg'
      r2Ptr % isVarArray = .false.
      r2Ptr % isDecomposed = .true.
      r2Ptr % hasTimeDimension = .true.
      r2Ptr % isPersistent = .true.
      r2Ptr % isActive = .false.
! Setting up dimensions
      r2Ptr % dimNames(1) = 'nFGSoilLevels'
      r2Ptr % dimNames(2) = 'nCells'
      r2Ptr % defaultValue = 0.0
      allocate(r2Ptr % attLists(1))
      allocate(r2Ptr % attLists(1) % attList)
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'units', 'k')
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'long_name', 'First guess soil temperature')
      r2Ptr % missingValue = MPAS_REAL_FILLVAL
      r2Ptr % block => block

      if (first_guess_fieldActive) then
         r2Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'st_fg', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'st_fg', r2Ptr)

! Define variable sm_fg
      allocate(r2Ptr)

! Setting up time level 1
      r2Ptr % fieldName = 'sm_fg'
      r2Ptr % isVarArray = .false.
      r2Ptr % isDecomposed = .true.
      r2Ptr % hasTimeDimension = .true.
      r2Ptr % isPersistent = .true.
      r2Ptr % isActive = .false.
! Setting up dimensions
      r2Ptr % dimNames(1) = 'nFGSoilLevels'
      r2Ptr % dimNames(2) = 'nCells'
      r2Ptr % defaultValue = 0.0
      allocate(r2Ptr % attLists(1))
      allocate(r2Ptr % attLists(1) % attList)
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'units', 'm^{3} m^{-3}')
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'long_name', 'First guess soil moisture')
      r2Ptr % missingValue = MPAS_REAL_FILLVAL
      r2Ptr % block => block

      if (first_guess_fieldActive) then
         r2Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'sm_fg', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'sm_fg', r2Ptr)

! Define variable dz
      allocate(r2Ptr)

! Setting up time level 1
      r2Ptr % fieldName = 'dz'
      r2Ptr % isVarArray = .false.
      r2Ptr % isDecomposed = .true.
      r2Ptr % hasTimeDimension = .true.
      r2Ptr % isPersistent = .true.
      r2Ptr % isActive = .false.
! Setting up dimensions
      r2Ptr % dimNames(1) = 'nSoilLevels'
      r2Ptr % dimNames(2) = 'nCells'
      r2Ptr % defaultValue = 0.0
      allocate(r2Ptr % attLists(1))
      allocate(r2Ptr % attLists(1) % attList)
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'units', 'm')
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'long_name', 'depth of soil layer bottom')
      r2Ptr % missingValue = MPAS_REAL_FILLVAL
      r2Ptr % block => block

      if (met_stage_outActive) then
         r2Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'dz', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'dz', r2Ptr)

! Define variable dzs
      allocate(r2Ptr)

! Setting up time level 1
      r2Ptr % fieldName = 'dzs'
      r2Ptr % isVarArray = .false.
      r2Ptr % isDecomposed = .true.
      r2Ptr % hasTimeDimension = .true.
      r2Ptr % isPersistent = .true.
      r2Ptr % isActive = .false.
! Setting up dimensions
      r2Ptr % dimNames(1) = 'nSoilLevels'
      r2Ptr % dimNames(2) = 'nCells'
      r2Ptr % defaultValue = 0.0
      allocate(r2Ptr % attLists(1))
      allocate(r2Ptr % attLists(1) % attList)
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'units', 'm')
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'long_name', 'soil layer thickness')
      r2Ptr % missingValue = MPAS_REAL_FILLVAL
      r2Ptr % block => block

      if (met_stage_outActive) then
         r2Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'dzs', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'dzs', r2Ptr)

! Define variable zs
      allocate(r2Ptr)

! Setting up time level 1
      r2Ptr % fieldName = 'zs'
      r2Ptr % isVarArray = .false.
      r2Ptr % isDecomposed = .true.
      r2Ptr % hasTimeDimension = .true.
      r2Ptr % isPersistent = .true.
      r2Ptr % isActive = .false.
! Setting up dimensions
      r2Ptr % dimNames(1) = 'nSoilLevels'
      r2Ptr % dimNames(2) = 'nCells'
      r2Ptr % defaultValue = 0.0
      allocate(r2Ptr % attLists(1))
      allocate(r2Ptr % attLists(1) % attList)
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'units', 'm')
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'long_name', 'depth of centers of soil layers')
      r2Ptr % missingValue = MPAS_REAL_FILLVAL
      r2Ptr % block => block

      if (met_stage_outActive) then
         r2Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'zs', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'zs', r2Ptr)

! Define variable sh2o
      allocate(r2Ptr)

! Setting up time level 1
      r2Ptr % fieldName = 'sh2o'
      r2Ptr % isVarArray = .false.
      r2Ptr % isDecomposed = .true.
      r2Ptr % hasTimeDimension = .true.
      r2Ptr % isPersistent = .true.
      r2Ptr % isActive = .false.
! Setting up dimensions
      r2Ptr % dimNames(1) = 'nSoilLevels'
      r2Ptr % dimNames(2) = 'nCells'
      r2Ptr % defaultValue = 0.0
      allocate(r2Ptr % attLists(1))
      allocate(r2Ptr % attLists(1) % attList)
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'units', 'm3 m^{-3}')
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'long_name', 'soil equivalent liquid water ')
      r2Ptr % missingValue = MPAS_REAL_FILLVAL
      r2Ptr % block => block

      if (met_stage_outActive) then
         r2Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'sh2o', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'sh2o', r2Ptr)

! Define variable smois
      allocate(r2Ptr)

! Setting up time level 1
      r2Ptr % fieldName = 'smois'
      r2Ptr % isVarArray = .false.
      r2Ptr % isDecomposed = .true.
      r2Ptr % hasTimeDimension = .true.
      r2Ptr % isPersistent = .true.
      r2Ptr % isActive = .false.
! Setting up dimensions
      r2Ptr % dimNames(1) = 'nSoilLevels'
      r2Ptr % dimNames(2) = 'nCells'
      r2Ptr % defaultValue = 0.0
      allocate(r2Ptr % attLists(1))
      allocate(r2Ptr % attLists(1) % attList)
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'units', 'm3 m^{-3}')
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'long_name', 'soil moisture')
      r2Ptr % missingValue = MPAS_REAL_FILLVAL
      r2Ptr % block => block

      if (met_stage_outActive) then
         r2Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'smois', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'smois', r2Ptr)

! Define variable tslb
      allocate(r2Ptr)

! Setting up time level 1
      r2Ptr % fieldName = 'tslb'
      r2Ptr % isVarArray = .false.
      r2Ptr % isDecomposed = .true.
      r2Ptr % hasTimeDimension = .true.
      r2Ptr % isPersistent = .true.
      r2Ptr % isActive = .false.
! Setting up dimensions
      r2Ptr % dimNames(1) = 'nSoilLevels'
      r2Ptr % dimNames(2) = 'nCells'
      r2Ptr % defaultValue = 0.0
      allocate(r2Ptr % attLists(1))
      allocate(r2Ptr % attLists(1) % attList)
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'units', 'K')
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'long_name', 'soil layer temperature')
      r2Ptr % missingValue = MPAS_REAL_FILLVAL
      r2Ptr % block => block

      if (met_stage_outActive) then
         r2Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'tslb', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'tslb', r2Ptr)

! Define variable smcrel
      allocate(r2Ptr)

! Setting up time level 1
      r2Ptr % fieldName = 'smcrel'
      r2Ptr % isVarArray = .false.
      r2Ptr % isDecomposed = .true.
      r2Ptr % hasTimeDimension = .true.
      r2Ptr % isPersistent = .true.
      r2Ptr % isActive = .false.
! Setting up dimensions
      r2Ptr % dimNames(1) = 'nSoilLevels'
      r2Ptr % dimNames(2) = 'nCells'
      r2Ptr % defaultValue = 0.0
      allocate(r2Ptr % attLists(1))
      allocate(r2Ptr % attLists(1) % attList)
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'units', 'm3 m^{-3}')
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'long_name', 'soil moisture threshold below which transpiration begins to '&
//'stress')
      r2Ptr % missingValue = MPAS_REAL_FILLVAL
      r2Ptr % block => block

      if (met_stage_outActive) then
         r2Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'smcrel', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'smcrel', r2Ptr)

! Define variable tmn
      allocate(r1Ptr)

! Setting up time level 1
      r1Ptr % fieldName = 'tmn'
      r1Ptr % isVarArray = .false.
      r1Ptr % isDecomposed = .true.
      r1Ptr % hasTimeDimension = .true.
      r1Ptr % isPersistent = .true.
      r1Ptr % isActive = .false.
! Setting up dimensions
      r1Ptr % dimNames(1) = 'nCells'
      r1Ptr % defaultValue = 0.0
      allocate(r1Ptr % attLists(1))
      allocate(r1Ptr % attLists(1) % attList)
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'units', 'K')
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'long_name', 'deep soil temperature')
      r1Ptr % missingValue = MPAS_REAL_FILLVAL
      r1Ptr % block => block

      if (met_stage_outActive) then
         r1Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'tmn', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'tmn', r1Ptr)

! Define variable skintemp
      allocate(r1Ptr)

! Setting up time level 1
      r1Ptr % fieldName = 'skintemp'
      r1Ptr % isVarArray = .false.
      r1Ptr % isDecomposed = .true.
      r1Ptr % hasTimeDimension = .true.
      r1Ptr % isPersistent = .true.
      r1Ptr % isActive = .false.
! Setting up dimensions
      r1Ptr % dimNames(1) = 'nCells'
      r1Ptr % defaultValue = 0.0
      allocate(r1Ptr % attLists(1))
      allocate(r1Ptr % attLists(1) % attList)
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'units', 'K')
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'long_name', 'ground or water surface temperature')
      r1Ptr % missingValue = MPAS_REAL_FILLVAL
      r1Ptr % block => block

      if (met_stage_outActive) then
         r1Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'skintemp', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'skintemp', r1Ptr)

! Define variable sst
      allocate(r1Ptr)

! Setting up time level 1
      r1Ptr % fieldName = 'sst'
      r1Ptr % isVarArray = .false.
      r1Ptr % isDecomposed = .true.
      r1Ptr % hasTimeDimension = .true.
      r1Ptr % isPersistent = .true.
      r1Ptr % isActive = .false.
! Setting up dimensions
      r1Ptr % dimNames(1) = 'nCells'
      r1Ptr % defaultValue = 0.0
      allocate(r1Ptr % attLists(1))
      allocate(r1Ptr % attLists(1) % attList)
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'units', 'K')
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'long_name', 'sea-surface temperature')
      r1Ptr % missingValue = MPAS_REAL_FILLVAL
      r1Ptr % block => block

      if (met_stage_outActive .or. sfc_updateActive) then
         r1Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'sst', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'sst', r1Ptr)

! Define variable snow
      allocate(r1Ptr)

! Setting up time level 1
      r1Ptr % fieldName = 'snow'
      r1Ptr % isVarArray = .false.
      r1Ptr % isDecomposed = .true.
      r1Ptr % hasTimeDimension = .true.
      r1Ptr % isPersistent = .true.
      r1Ptr % isActive = .false.
! Setting up dimensions
      r1Ptr % dimNames(1) = 'nCells'
      r1Ptr % defaultValue = 0.0
      allocate(r1Ptr % attLists(1))
      allocate(r1Ptr % attLists(1) % attList)
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'units', 'kg m^{-2}')
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'long_name', 'snow water equivalent')
      r1Ptr % missingValue = MPAS_REAL_FILLVAL
      r1Ptr % block => block

      if (met_stage_outActive) then
         r1Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'snow', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'snow', r1Ptr)

! Define variable snowc
      allocate(r1Ptr)

! Setting up time level 1
      r1Ptr % fieldName = 'snowc'
      r1Ptr % isVarArray = .false.
      r1Ptr % isDecomposed = .true.
      r1Ptr % hasTimeDimension = .true.
      r1Ptr % isPersistent = .true.
      r1Ptr % isActive = .false.
! Setting up dimensions
      r1Ptr % dimNames(1) = 'nCells'
      r1Ptr % defaultValue = 0.0
      allocate(r1Ptr % attLists(1))
      allocate(r1Ptr % attLists(1) % attList)
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'long_name', 'flag for snow on ground (=0 no snow; =1,otherwise')
      r1Ptr % missingValue = MPAS_REAL_FILLVAL
      r1Ptr % block => block

      if (met_stage_outActive) then
         r1Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'snowc', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'snowc', r1Ptr)

! Define variable snowh
      allocate(r1Ptr)

! Setting up time level 1
      r1Ptr % fieldName = 'snowh'
      r1Ptr % isVarArray = .false.
      r1Ptr % isDecomposed = .true.
      r1Ptr % hasTimeDimension = .true.
      r1Ptr % isPersistent = .true.
      r1Ptr % isActive = .false.
! Setting up dimensions
      r1Ptr % dimNames(1) = 'nCells'
      r1Ptr % defaultValue = 0.0
      allocate(r1Ptr % attLists(1))
      allocate(r1Ptr % attLists(1) % attList)
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'units', 'm')
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'long_name', 'physical snow depth')
      r1Ptr % missingValue = MPAS_REAL_FILLVAL
      r1Ptr % block => block

      if (met_stage_outActive) then
         r1Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'snowh', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'snowh', r1Ptr)

! Define variable xice
      allocate(r1Ptr)

! Setting up time level 1
      r1Ptr % fieldName = 'xice'
      r1Ptr % isVarArray = .false.
      r1Ptr % isDecomposed = .true.
      r1Ptr % hasTimeDimension = .true.
      r1Ptr % isPersistent = .true.
      r1Ptr % isActive = .false.
! Setting up dimensions
      r1Ptr % dimNames(1) = 'nCells'
      r1Ptr % defaultValue = 0.0
      allocate(r1Ptr % attLists(1))
      allocate(r1Ptr % attLists(1) % attList)
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'long_name', 'fractional area coverage of sea-ice')
      r1Ptr % missingValue = MPAS_REAL_FILLVAL
      r1Ptr % block => block

      if (met_stage_outActive .or. sfc_updateActive) then
         r1Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'xice', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'xice', r1Ptr)

! Define variable seaice
      allocate(r1Ptr)

! Setting up time level 1
      r1Ptr % fieldName = 'seaice'
      r1Ptr % isVarArray = .false.
      r1Ptr % isDecomposed = .true.
      r1Ptr % hasTimeDimension = .true.
      r1Ptr % isPersistent = .true.
      r1Ptr % isActive = .false.
! Setting up dimensions
      r1Ptr % dimNames(1) = 'nCells'
      r1Ptr % defaultValue = 0.0
      allocate(r1Ptr % attLists(1))
      allocate(r1Ptr % attLists(1) % attList)
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'long_name', 'sea-ice flag (0=no seaice; =1 otherwise)')
      r1Ptr % missingValue = MPAS_REAL_FILLVAL
      r1Ptr % block => block

      if (met_stage_outActive) then
         r1Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'seaice', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'seaice', r1Ptr)

! Define variable gfs_z
      allocate(r2Ptr)

! Setting up time level 1
      r2Ptr % fieldName = 'gfs_z'
      r2Ptr % isVarArray = .false.
      r2Ptr % isDecomposed = .true.
      r2Ptr % hasTimeDimension = .true.
      r2Ptr % isPersistent = .true.
      r2Ptr % isActive = .false.
! Setting up dimensions
      r2Ptr % dimNames(1) = 'nVertLevels'
      r2Ptr % dimNames(2) = 'nCells'
      r2Ptr % defaultValue = 0.0
      allocate(r2Ptr % attLists(1))
      allocate(r2Ptr % attLists(1) % attList)
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'units', 'm')
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'long_name', 'geopotential height vertically interpolated from first guess')
      r2Ptr % missingValue = MPAS_REAL_FILLVAL
      r2Ptr % block => block

      if (met_stage_outActive) then
         r2Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'gfs_z', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'gfs_z', r2Ptr)

! Define variable vegfra
      allocate(r1Ptr)

! Setting up time level 1
      r1Ptr % fieldName = 'vegfra'
      r1Ptr % isVarArray = .false.
      r1Ptr % isDecomposed = .true.
      r1Ptr % hasTimeDimension = .true.
      r1Ptr % isPersistent = .true.
      r1Ptr % isActive = .false.
! Setting up dimensions
      r1Ptr % dimNames(1) = 'nCells'
      r1Ptr % defaultValue = 0.0
      allocate(r1Ptr % attLists(1))
      allocate(r1Ptr % attLists(1) % attList)
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'long_name', 'vegetation fraction')
      r1Ptr % missingValue = MPAS_REAL_FILLVAL
      r1Ptr % block => block

      if (met_stage_outActive) then
         r1Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'vegfra', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'vegfra', r1Ptr)

! Define variable sfc_albbck
      allocate(r1Ptr)

! Setting up time level 1
      r1Ptr % fieldName = 'sfc_albbck'
      r1Ptr % isVarArray = .false.
      r1Ptr % isDecomposed = .true.
      r1Ptr % hasTimeDimension = .true.
      r1Ptr % isPersistent = .true.
      r1Ptr % isActive = .false.
! Setting up dimensions
      r1Ptr % dimNames(1) = 'nCells'
      r1Ptr % defaultValue = 0.0
      allocate(r1Ptr % attLists(1))
      allocate(r1Ptr % attLists(1) % attList)
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'long_name', 'background surface albedo')
      r1Ptr % missingValue = MPAS_REAL_FILLVAL
      r1Ptr % block => block

      if (met_stage_outActive) then
         r1Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'sfc_albbck', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'sfc_albbck', r1Ptr)

! Define variable xland
      allocate(r1Ptr)

! Setting up time level 1
      r1Ptr % fieldName = 'xland'
      r1Ptr % isVarArray = .false.
      r1Ptr % isDecomposed = .true.
      r1Ptr % hasTimeDimension = .true.
      r1Ptr % isPersistent = .true.
      r1Ptr % isActive = .false.
! Setting up dimensions
      r1Ptr % dimNames(1) = 'nCells'
      r1Ptr % defaultValue = 0.0
      allocate(r1Ptr % attLists(1))
      allocate(r1Ptr % attLists(1) % attList)
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'long_name', 'land-ocean mask (1=land including sea-ice ; 2=ocean)')
      r1Ptr % missingValue = MPAS_REAL_FILLVAL
      r1Ptr % block => block

      if (met_stage_outActive) then
         r1Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'xland', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'xland', r1Ptr)

! Define variable u10
      allocate(r1Ptr)

! Setting up time level 1
      r1Ptr % fieldName = 'u10'
      r1Ptr % isVarArray = .false.
      r1Ptr % isDecomposed = .true.
      r1Ptr % hasTimeDimension = .true.
      r1Ptr % isPersistent = .true.
      r1Ptr % isActive = .false.
! Setting up dimensions
      r1Ptr % dimNames(1) = 'nCells'
      r1Ptr % defaultValue = 0.0
      allocate(r1Ptr % attLists(1))
      allocate(r1Ptr % attLists(1) % attList)
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'units', 'm s^{-1}')
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'long_name', '10-meter zonal wind')
      r1Ptr % missingValue = MPAS_REAL_FILLVAL
      r1Ptr % block => block

      if (met_stage_outActive) then
         r1Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'u10', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'u10', r1Ptr)

! Define variable v10
      allocate(r1Ptr)

! Setting up time level 1
      r1Ptr % fieldName = 'v10'
      r1Ptr % isVarArray = .false.
      r1Ptr % isDecomposed = .true.
      r1Ptr % hasTimeDimension = .true.
      r1Ptr % isPersistent = .true.
      r1Ptr % isActive = .false.
! Setting up dimensions
      r1Ptr % dimNames(1) = 'nCells'
      r1Ptr % defaultValue = 0.0
      allocate(r1Ptr % attLists(1))
      allocate(r1Ptr % attLists(1) % attList)
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'units', 'm s^{-1}')
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'long_name', '10-meter meridional wind')
      r1Ptr % missingValue = MPAS_REAL_FILLVAL
      r1Ptr % block => block

      if (met_stage_outActive) then
         r1Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'v10', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'v10', r1Ptr)

! Define variable q2
      allocate(r1Ptr)

! Setting up time level 1
      r1Ptr % fieldName = 'q2'
      r1Ptr % isVarArray = .false.
      r1Ptr % isDecomposed = .true.
      r1Ptr % hasTimeDimension = .true.
      r1Ptr % isPersistent = .true.
      r1Ptr % isActive = .false.
! Setting up dimensions
      r1Ptr % dimNames(1) = 'nCells'
      r1Ptr % defaultValue = 0.0
      allocate(r1Ptr % attLists(1))
      allocate(r1Ptr % attLists(1) % attList)
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'units', 'kg kg^{-1}')
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'long_name', '2-meter specific humidity')
      r1Ptr % missingValue = MPAS_REAL_FILLVAL
      r1Ptr % block => block

      if (met_stage_outActive) then
         r1Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'q2', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'q2', r1Ptr)

! Define variable rh2
      allocate(r1Ptr)

! Setting up time level 1
      r1Ptr % fieldName = 'rh2'
      r1Ptr % isVarArray = .false.
      r1Ptr % isDecomposed = .true.
      r1Ptr % hasTimeDimension = .true.
      r1Ptr % isPersistent = .true.
      r1Ptr % isActive = .false.
! Setting up dimensions
      r1Ptr % dimNames(1) = 'nCells'
      r1Ptr % defaultValue = 0.0
      allocate(r1Ptr % attLists(1))
      allocate(r1Ptr % attLists(1) % attList)
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'units', 'percent')
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'long_name', '2-meter relative humidity')
      r1Ptr % missingValue = MPAS_REAL_FILLVAL
      r1Ptr % block => block

      if (met_stage_outActive) then
         r1Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'rh2', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'rh2', r1Ptr)

! Define variable t2m
      allocate(r1Ptr)

! Setting up time level 1
      r1Ptr % fieldName = 't2m'
      r1Ptr % isVarArray = .false.
      r1Ptr % isDecomposed = .true.
      r1Ptr % hasTimeDimension = .true.
      r1Ptr % isPersistent = .true.
      r1Ptr % isActive = .false.
! Setting up dimensions
      r1Ptr % dimNames(1) = 'nCells'
      r1Ptr % defaultValue = 0.0
      allocate(r1Ptr % attLists(1))
      allocate(r1Ptr % attLists(1) % attList)
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'units', 'K')
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'long_name', '2-meter temperature')
      r1Ptr % missingValue = MPAS_REAL_FILLVAL
      r1Ptr % block => block

      if (met_stage_outActive) then
         r1Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 't2m', r1Ptr)
      call mpas_pool_add_field(block % allFields, 't2m', r1Ptr)



      if (associated(newSubPool)) then
         call mpas_pool_add_config(newSubPool, 'on_a_sphere', block % domain % on_a_sphere)
         call mpas_pool_add_config(newSubPool, 'sphere_radius', block % domain % sphere_radius)
         call mpas_pool_add_config(newSubPool, 'is_periodic', block % domain % is_periodic)
         call mpas_pool_add_config(newSubPool, 'x_period', block % domain % x_period)
         call mpas_pool_add_config(newSubPool, 'y_period', block % domain % y_period)
      end if

   end subroutine init_atm_generate_pool_fg


   subroutine init_atm_generate_pool_diag(block, structPool, dimensionPool, packagePool)
      use mpas_derived_types
      use mpas_pool_routines
      use mpas_io_units
      use mpas_io, only : MPAS_REAL_FILLVAL, MPAS_INT_FILLVAL, MPAS_CHAR_FILLVAL
      implicit none
      type (block_type), intent(inout), pointer :: block
      type (mpas_pool_type), intent(inout) :: structPool
      type (mpas_pool_type), intent(inout) :: dimensionPool
      type (mpas_pool_type), intent(in) :: packagePool

      type (field0DReal), pointer :: r0Ptr
      type (field1DReal), pointer :: r1Ptr
      type (field2DReal), pointer :: r2Ptr
      type (field3DReal), pointer :: r3Ptr
      type (field4DReal), pointer :: r4Ptr
      type (field5DReal), pointer :: r5Ptr
      type (field0DInteger), pointer :: i0Ptr
      type (field1DInteger), pointer :: i1Ptr
      type (field2DInteger), pointer :: i2Ptr
      type (field3DInteger), pointer :: i3Ptr
      type (field0DChar), pointer :: c0Ptr
      type (field1DChar), pointer :: c1Ptr
      type (field0DReal), dimension(:), pointer :: r0aPtr
      type (field1DReal), dimension(:), pointer :: r1aPtr
      type (field2DReal), dimension(:), pointer :: r2aPtr
      type (field3DReal), dimension(:), pointer :: r3aPtr
      type (field4DReal), dimension(:), pointer :: r4aPtr
      type (field5DReal), dimension(:), pointer :: r5aPtr
      type (field0DInteger), dimension(:), pointer :: i0aPtr
      type (field1DInteger), dimension(:), pointer :: i1aPtr
      type (field2DInteger), dimension(:), pointer :: i2aPtr
      type (field3DInteger), dimension(:), pointer :: i3aPtr
      type (field0DChar), dimension(:), pointer :: c0aPtr
      type (field1DChar), dimension(:), pointer :: c1aPtr

      type (mpas_pool_type), pointer :: newSubPool
      integer :: group_counter
      logical :: group_started
      integer :: group_start
      integer :: index_counter
      integer, pointer :: const_index

      logical, pointer :: initial_condsActive
      logical, pointer :: sfc_updateActive
      logical, pointer :: lbcsActive
      logical, pointer :: gwd_stage_inActive
      logical, pointer :: gwd_stage_outActive
      logical, pointer :: vertical_stage_inActive
      logical, pointer :: vertical_stage_outActive
      logical, pointer :: met_stage_inActive
      logical, pointer :: met_stage_outActive
      logical, pointer :: first_guess_fieldActive


      integer :: numConstituents

      nullify(newSubPool)
      group_counter = -1
      group_started = .false.
      group_start = -1
      call mpas_pool_get_package(packagePool, 'initial_condsActive', initial_condsActive)
      call mpas_pool_get_package(packagePool, 'sfc_updateActive', sfc_updateActive)
      call mpas_pool_get_package(packagePool, 'lbcsActive', lbcsActive)
      call mpas_pool_get_package(packagePool, 'gwd_stage_inActive', gwd_stage_inActive)
      call mpas_pool_get_package(packagePool, 'gwd_stage_outActive', gwd_stage_outActive)
      call mpas_pool_get_package(packagePool, 'vertical_stage_inActive', vertical_stage_inActive)
      call mpas_pool_get_package(packagePool, 'vertical_stage_outActive', vertical_stage_outActive)
      call mpas_pool_get_package(packagePool, 'met_stage_inActive', met_stage_inActive)
      call mpas_pool_get_package(packagePool, 'met_stage_outActive', met_stage_outActive)
      call mpas_pool_get_package(packagePool, 'first_guess_fieldActive', first_guess_fieldActive)

      allocate(newSubPool)
      call mpas_pool_create_pool(newSubPool)
      call mpas_pool_add_subpool(structPool, 'diag', newSubPool)
      call mpas_pool_add_subpool(block % allStructs, 'diag', newSubPool)

! Define variable pressure_p
      allocate(r2Ptr)

! Setting up time level 1
      r2Ptr % fieldName = 'pressure_p'
      r2Ptr % isVarArray = .false.
      r2Ptr % isDecomposed = .true.
      r2Ptr % hasTimeDimension = .true.
      r2Ptr % isPersistent = .true.
      r2Ptr % isActive = .false.
! Setting up dimensions
      r2Ptr % dimNames(1) = 'nVertLevels'
      r2Ptr % dimNames(2) = 'nCells'
      r2Ptr % defaultValue = 0.0
      allocate(r2Ptr % attLists(1))
      allocate(r2Ptr % attLists(1) % attList)
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'units', 'Pa')
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'long_name', 'Perturbation pressure')
      r2Ptr % missingValue = MPAS_REAL_FILLVAL
      r2Ptr % block => block

      if (met_stage_outActive) then
         r2Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'pressure_p', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'pressure_p', r2Ptr)

! Define variable rho
      allocate(r2Ptr)

! Setting up time level 1
      r2Ptr % fieldName = 'rho'
      r2Ptr % isVarArray = .false.
      r2Ptr % isDecomposed = .true.
      r2Ptr % hasTimeDimension = .true.
      r2Ptr % isPersistent = .true.
      r2Ptr % isActive = .false.
! Setting up dimensions
      r2Ptr % dimNames(1) = 'nVertLevels'
      r2Ptr % dimNames(2) = 'nCells'
      r2Ptr % defaultValue = 0.0
      allocate(r2Ptr % attLists(1))
      allocate(r2Ptr % attLists(1) % attList)
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'units', 'kg m^{-3}')
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'long_name', 'Dry air density')
      r2Ptr % missingValue = MPAS_REAL_FILLVAL
      r2Ptr % block => block

      if (met_stage_outActive) then
         r2Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'rho', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'rho', r2Ptr)

! Define variable theta
      allocate(r2Ptr)

! Setting up time level 1
      r2Ptr % fieldName = 'theta'
      r2Ptr % isVarArray = .false.
      r2Ptr % isDecomposed = .true.
      r2Ptr % hasTimeDimension = .true.
      r2Ptr % isPersistent = .true.
      r2Ptr % isActive = .false.
! Setting up dimensions
      r2Ptr % dimNames(1) = 'nVertLevels'
      r2Ptr % dimNames(2) = 'nCells'
      r2Ptr % defaultValue = 0.0
      allocate(r2Ptr % attLists(1))
      allocate(r2Ptr % attLists(1) % attList)
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'units', 'K')
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'long_name', 'Potential temperature')
      r2Ptr % missingValue = MPAS_REAL_FILLVAL
      r2Ptr % block => block

      if (met_stage_outActive) then
         r2Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'theta', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'theta', r2Ptr)

! Define variable v
      allocate(r2Ptr)

! Setting up time level 1
      r2Ptr % fieldName = 'v'
      r2Ptr % isVarArray = .false.
      r2Ptr % isDecomposed = .true.
      r2Ptr % hasTimeDimension = .true.
      r2Ptr % isPersistent = .true.
      r2Ptr % isActive = .false.
! Setting up dimensions
      r2Ptr % dimNames(1) = 'nVertLevels'
      r2Ptr % dimNames(2) = 'nEdges'
      r2Ptr % defaultValue = 0.0
      allocate(r2Ptr % attLists(1))
      allocate(r2Ptr % attLists(1) % attList)
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'units', 'm s^{-1}')
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'long_name', 'Horizontal tangential velocity at edges')
      r2Ptr % missingValue = MPAS_REAL_FILLVAL
      r2Ptr % block => block

      if (met_stage_outActive) then
         r2Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'v', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'v', r2Ptr)

! Define variable relhum
      allocate(r2Ptr)

! Setting up time level 1
      r2Ptr % fieldName = 'relhum'
      r2Ptr % isVarArray = .false.
      r2Ptr % isDecomposed = .true.
      r2Ptr % hasTimeDimension = .true.
      r2Ptr % isPersistent = .true.
      r2Ptr % isActive = .false.
! Setting up dimensions
      r2Ptr % dimNames(1) = 'nVertLevels'
      r2Ptr % dimNames(2) = 'nCells'
      r2Ptr % defaultValue = 0.0
      allocate(r2Ptr % attLists(1))
      allocate(r2Ptr % attLists(1) % attList)
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'units', 'percent')
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'long_name', 'Relative humidity')
      r2Ptr % missingValue = MPAS_REAL_FILLVAL
      r2Ptr % block => block

      if (met_stage_outActive) then
         r2Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'relhum', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'relhum', r2Ptr)

! Define variable spechum
      allocate(r2Ptr)

! Setting up time level 1
      r2Ptr % fieldName = 'spechum'
      r2Ptr % isVarArray = .false.
      r2Ptr % isDecomposed = .true.
      r2Ptr % hasTimeDimension = .true.
      r2Ptr % isPersistent = .true.
      r2Ptr % isActive = .false.
! Setting up dimensions
      r2Ptr % dimNames(1) = 'nVertLevels'
      r2Ptr % dimNames(2) = 'nCells'
      r2Ptr % defaultValue = 0.0
      allocate(r2Ptr % attLists(1))
      allocate(r2Ptr % attLists(1) % attList)
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'units', 'kg kg^{-1}')
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'long_name', 'Specific humidity')
      r2Ptr % missingValue = MPAS_REAL_FILLVAL
      r2Ptr % block => block

      if (met_stage_outActive) then
         r2Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'spechum', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'spechum', r2Ptr)

! Define variable uReconstructX
      allocate(r2Ptr)

! Setting up time level 1
      r2Ptr % fieldName = 'uReconstructX'
      r2Ptr % isVarArray = .false.
      r2Ptr % isDecomposed = .true.
      r2Ptr % hasTimeDimension = .true.
      r2Ptr % isPersistent = .true.
      r2Ptr % isActive = .false.
! Setting up dimensions
      r2Ptr % dimNames(1) = 'nVertLevels'
      r2Ptr % dimNames(2) = 'nCells'
      r2Ptr % defaultValue = 0.0
      allocate(r2Ptr % attLists(1))
      allocate(r2Ptr % attLists(1) % attList)
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'units', 'm s^{-1}')
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'long_name', 'Cartesian x-component of reconstructed horizontal velocity at '&
//'cell centers')
      r2Ptr % missingValue = MPAS_REAL_FILLVAL
      r2Ptr % block => block

      if (met_stage_outActive) then
         r2Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'uReconstructX', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'uReconstructX', r2Ptr)

! Define variable uReconstructY
      allocate(r2Ptr)

! Setting up time level 1
      r2Ptr % fieldName = 'uReconstructY'
      r2Ptr % isVarArray = .false.
      r2Ptr % isDecomposed = .true.
      r2Ptr % hasTimeDimension = .true.
      r2Ptr % isPersistent = .true.
      r2Ptr % isActive = .false.
! Setting up dimensions
      r2Ptr % dimNames(1) = 'nVertLevels'
      r2Ptr % dimNames(2) = 'nCells'
      r2Ptr % defaultValue = 0.0
      allocate(r2Ptr % attLists(1))
      allocate(r2Ptr % attLists(1) % attList)
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'units', 'm s^{-1}')
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'long_name', 'Cartesian y-component of reconstructed horizontal velocity at '&
//'cell centers')
      r2Ptr % missingValue = MPAS_REAL_FILLVAL
      r2Ptr % block => block

      if (met_stage_outActive) then
         r2Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'uReconstructY', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'uReconstructY', r2Ptr)

! Define variable uReconstructZ
      allocate(r2Ptr)

! Setting up time level 1
      r2Ptr % fieldName = 'uReconstructZ'
      r2Ptr % isVarArray = .false.
      r2Ptr % isDecomposed = .true.
      r2Ptr % hasTimeDimension = .true.
      r2Ptr % isPersistent = .true.
      r2Ptr % isActive = .false.
! Setting up dimensions
      r2Ptr % dimNames(1) = 'nVertLevels'
      r2Ptr % dimNames(2) = 'nCells'
      r2Ptr % defaultValue = 0.0
      allocate(r2Ptr % attLists(1))
      allocate(r2Ptr % attLists(1) % attList)
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'units', 'm s^{-1}')
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'long_name', 'Cartesian z-component of reconstructed horizontal velocity at '&
//'cell centers')
      r2Ptr % missingValue = MPAS_REAL_FILLVAL
      r2Ptr % block => block

      if (met_stage_outActive) then
         r2Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'uReconstructZ', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'uReconstructZ', r2Ptr)

! Define variable uReconstructZonal
      allocate(r2Ptr)

! Setting up time level 1
      r2Ptr % fieldName = 'uReconstructZonal'
      r2Ptr % isVarArray = .false.
      r2Ptr % isDecomposed = .true.
      r2Ptr % hasTimeDimension = .true.
      r2Ptr % isPersistent = .true.
      r2Ptr % isActive = .false.
! Setting up dimensions
      r2Ptr % dimNames(1) = 'nVertLevels'
      r2Ptr % dimNames(2) = 'nCells'
      r2Ptr % defaultValue = 0.0
      allocate(r2Ptr % attLists(1))
      allocate(r2Ptr % attLists(1) % attList)
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'units', 'm s^{-1}')
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'long_name', 'Zonal component of reconstructed horizontal velocity at cell '&
//'centers')
      r2Ptr % missingValue = MPAS_REAL_FILLVAL
      r2Ptr % block => block

      if (met_stage_outActive) then
         r2Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'uReconstructZonal', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'uReconstructZonal', r2Ptr)

! Define variable uReconstructMeridional
      allocate(r2Ptr)

! Setting up time level 1
      r2Ptr % fieldName = 'uReconstructMeridional'
      r2Ptr % isVarArray = .false.
      r2Ptr % isDecomposed = .true.
      r2Ptr % hasTimeDimension = .true.
      r2Ptr % isPersistent = .true.
      r2Ptr % isActive = .false.
! Setting up dimensions
      r2Ptr % dimNames(1) = 'nVertLevels'
      r2Ptr % dimNames(2) = 'nCells'
      r2Ptr % defaultValue = 0.0
      allocate(r2Ptr % attLists(1))
      allocate(r2Ptr % attLists(1) % attList)
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'units', 'm s^{-1}')
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'long_name', 'Meridional component of reconstructed horizontal velocity at '&
//'cell centers')
      r2Ptr % missingValue = MPAS_REAL_FILLVAL
      r2Ptr % block => block

      if (met_stage_outActive) then
         r2Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'uReconstructMeridional', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'uReconstructMeridional', r2Ptr)

! Define variable exner
      allocate(r2Ptr)

! Setting up time level 1
      r2Ptr % fieldName = 'exner'
      r2Ptr % isVarArray = .false.
      r2Ptr % isDecomposed = .true.
      r2Ptr % hasTimeDimension = .true.
      r2Ptr % isPersistent = .true.
      r2Ptr % isActive = .false.
! Setting up dimensions
      r2Ptr % dimNames(1) = 'nVertLevels'
      r2Ptr % dimNames(2) = 'nCells'
      r2Ptr % defaultValue = 0.0
      allocate(r2Ptr % attLists(1))
      allocate(r2Ptr % attLists(1) % attList)
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'long_name', 'Exner function')
      r2Ptr % missingValue = MPAS_REAL_FILLVAL
      r2Ptr % block => block

      if (met_stage_outActive) then
         r2Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'exner', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'exner', r2Ptr)

! Define variable exner_base
      allocate(r2Ptr)

! Setting up time level 1
      r2Ptr % fieldName = 'exner_base'
      r2Ptr % isVarArray = .false.
      r2Ptr % isDecomposed = .true.
      r2Ptr % hasTimeDimension = .true.
      r2Ptr % isPersistent = .true.
      r2Ptr % isActive = .false.
! Setting up dimensions
      r2Ptr % dimNames(1) = 'nVertLevels'
      r2Ptr % dimNames(2) = 'nCells'
      r2Ptr % defaultValue = 0.0
      allocate(r2Ptr % attLists(1))
      allocate(r2Ptr % attLists(1) % attList)
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'long_name', 'Base-state Exner function')
      r2Ptr % missingValue = MPAS_REAL_FILLVAL
      r2Ptr % block => block

      if (met_stage_outActive) then
         r2Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'exner_base', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'exner_base', r2Ptr)

! Define variable rtheta_base
      allocate(r2Ptr)

! Setting up time level 1
      r2Ptr % fieldName = 'rtheta_base'
      r2Ptr % isVarArray = .false.
      r2Ptr % isDecomposed = .true.
      r2Ptr % hasTimeDimension = .true.
      r2Ptr % isPersistent = .true.
      r2Ptr % isActive = .false.
! Setting up dimensions
      r2Ptr % dimNames(1) = 'nVertLevels'
      r2Ptr % dimNames(2) = 'nCells'
      r2Ptr % defaultValue = 0.0
      allocate(r2Ptr % attLists(1))
      allocate(r2Ptr % attLists(1) % attList)
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'units', 'kg K m^{-3}')
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'long_name', 'reference state rho*theta/zz')
      r2Ptr % missingValue = MPAS_REAL_FILLVAL
      r2Ptr % block => block

      if (met_stage_outActive) then
         r2Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'rtheta_base', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'rtheta_base', r2Ptr)

! Define variable pressure
      allocate(r2Ptr)

! Setting up time level 1
      r2Ptr % fieldName = 'pressure'
      r2Ptr % isVarArray = .false.
      r2Ptr % isDecomposed = .true.
      r2Ptr % hasTimeDimension = .true.
      r2Ptr % isPersistent = .true.
      r2Ptr % isActive = .false.
! Setting up dimensions
      r2Ptr % dimNames(1) = 'nVertLevels'
      r2Ptr % dimNames(2) = 'nCells'
      r2Ptr % defaultValue = 0.0
      allocate(r2Ptr % attLists(1))
      allocate(r2Ptr % attLists(1) % attList)
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'units', 'Pa')
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'long_name', 'Pressure')
      r2Ptr % missingValue = MPAS_REAL_FILLVAL
      r2Ptr % block => block

      if (met_stage_outActive) then
         r2Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'pressure', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'pressure', r2Ptr)

! Define variable pressure_base
      allocate(r2Ptr)

! Setting up time level 1
      r2Ptr % fieldName = 'pressure_base'
      r2Ptr % isVarArray = .false.
      r2Ptr % isDecomposed = .true.
      r2Ptr % hasTimeDimension = .true.
      r2Ptr % isPersistent = .true.
      r2Ptr % isActive = .false.
! Setting up dimensions
      r2Ptr % dimNames(1) = 'nVertLevels'
      r2Ptr % dimNames(2) = 'nCells'
      r2Ptr % defaultValue = 0.0
      allocate(r2Ptr % attLists(1))
      allocate(r2Ptr % attLists(1) % attList)
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'units', 'Pa')
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'long_name', 'Base state pressure')
      r2Ptr % missingValue = MPAS_REAL_FILLVAL
      r2Ptr % block => block

      if (met_stage_outActive) then
         r2Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'pressure_base', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'pressure_base', r2Ptr)

! Define variable rho_base
      allocate(r2Ptr)

! Setting up time level 1
      r2Ptr % fieldName = 'rho_base'
      r2Ptr % isVarArray = .false.
      r2Ptr % isDecomposed = .true.
      r2Ptr % hasTimeDimension = .true.
      r2Ptr % isPersistent = .true.
      r2Ptr % isActive = .false.
! Setting up dimensions
      r2Ptr % dimNames(1) = 'nVertLevels'
      r2Ptr % dimNames(2) = 'nCells'
      r2Ptr % defaultValue = 0.0
      allocate(r2Ptr % attLists(1))
      allocate(r2Ptr % attLists(1) % attList)
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'units', 'kg m^{-3}')
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'long_name', 'Base state dry air density')
      r2Ptr % missingValue = MPAS_REAL_FILLVAL
      r2Ptr % block => block

      if (met_stage_outActive) then
         r2Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'rho_base', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'rho_base', r2Ptr)

! Define variable theta_base
      allocate(r2Ptr)

! Setting up time level 1
      r2Ptr % fieldName = 'theta_base'
      r2Ptr % isVarArray = .false.
      r2Ptr % isDecomposed = .true.
      r2Ptr % hasTimeDimension = .true.
      r2Ptr % isPersistent = .true.
      r2Ptr % isActive = .false.
! Setting up dimensions
      r2Ptr % dimNames(1) = 'nVertLevels'
      r2Ptr % dimNames(2) = 'nCells'
      r2Ptr % defaultValue = 0.0
      allocate(r2Ptr % attLists(1))
      allocate(r2Ptr % attLists(1) % attList)
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'units', 'K')
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'long_name', 'Base state potential temperature')
      r2Ptr % missingValue = MPAS_REAL_FILLVAL
      r2Ptr % block => block

      if (met_stage_outActive) then
         r2Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'theta_base', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'theta_base', r2Ptr)

! Define variable cqw
      allocate(r2Ptr)

! Setting up time level 1
      r2Ptr % fieldName = 'cqw'
      r2Ptr % isVarArray = .false.
      r2Ptr % isDecomposed = .true.
      r2Ptr % hasTimeDimension = .true.
      r2Ptr % isPersistent = .true.
      r2Ptr % isActive = .false.
! Setting up dimensions
      r2Ptr % dimNames(1) = 'nVertLevels'
      r2Ptr % dimNames(2) = 'nCells'
      r2Ptr % defaultValue = 0.0
      allocate(r2Ptr % attLists(1))
      allocate(r2Ptr % attLists(1) % attList)
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'long_name', 'rho_d/rho_m at w points')
      r2Ptr % missingValue = MPAS_REAL_FILLVAL
      r2Ptr % block => block

      if (met_stage_outActive) then
         r2Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'cqw', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'cqw', r2Ptr)

! Define variable surface_pressure
      allocate(r1Ptr)

! Setting up time level 1
      r1Ptr % fieldName = 'surface_pressure'
      r1Ptr % isVarArray = .false.
      r1Ptr % isDecomposed = .true.
      r1Ptr % hasTimeDimension = .true.
      r1Ptr % isPersistent = .true.
      r1Ptr % isActive = .false.
! Setting up dimensions
      r1Ptr % dimNames(1) = 'nCells'
      r1Ptr % defaultValue = 0.0
      allocate(r1Ptr % attLists(1))
      allocate(r1Ptr % attLists(1) % attList)
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'units', 'Pa')
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'long_name', 'Diagnosed surface pressure')
      r1Ptr % missingValue = MPAS_REAL_FILLVAL
      r1Ptr % block => block

      if (met_stage_outActive) then
         r1Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'surface_pressure', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'surface_pressure', r1Ptr)

! Define variable ru
      allocate(r2Ptr)

! Setting up time level 1
      r2Ptr % fieldName = 'ru'
      r2Ptr % isVarArray = .false.
      r2Ptr % isDecomposed = .true.
      r2Ptr % hasTimeDimension = .true.
      r2Ptr % isPersistent = .true.
      r2Ptr % isActive = .false.
! Setting up dimensions
      r2Ptr % dimNames(1) = 'nVertLevels'
      r2Ptr % dimNames(2) = 'nEdges'
      r2Ptr % defaultValue = 0.0
      allocate(r2Ptr % attLists(1))
      allocate(r2Ptr % attLists(1) % attList)
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'units', 'kg m^{-2} s^{-1}')
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'long_name', 'horizontal momentum at cell edge (rho*u/zz)')
      r2Ptr % missingValue = MPAS_REAL_FILLVAL
      r2Ptr % block => block

      if (met_stage_outActive) then
         r2Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'ru', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'ru', r2Ptr)

! Define variable rw
      allocate(r2Ptr)

! Setting up time level 1
      r2Ptr % fieldName = 'rw'
      r2Ptr % isVarArray = .false.
      r2Ptr % isDecomposed = .true.
      r2Ptr % hasTimeDimension = .true.
      r2Ptr % isPersistent = .true.
      r2Ptr % isActive = .false.
! Setting up dimensions
      r2Ptr % dimNames(1) = 'nVertLevelsP1'
      r2Ptr % dimNames(2) = 'nCells'
      r2Ptr % defaultValue = 0.0
      allocate(r2Ptr % attLists(1))
      allocate(r2Ptr % attLists(1) % attList)
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'units', 'kg m^{-2} s^{-1}')
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'long_name', 'rho*omega/zz carried at w points')
      r2Ptr % missingValue = MPAS_REAL_FILLVAL
      r2Ptr % block => block

      if (met_stage_outActive) then
         r2Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'rw', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'rw', r2Ptr)

! Define variable rtheta_p
      allocate(r2Ptr)

! Setting up time level 1
      r2Ptr % fieldName = 'rtheta_p'
      r2Ptr % isVarArray = .false.
      r2Ptr % isDecomposed = .true.
      r2Ptr % hasTimeDimension = .true.
      r2Ptr % isPersistent = .true.
      r2Ptr % isActive = .false.
! Setting up dimensions
      r2Ptr % dimNames(1) = 'nVertLevels'
      r2Ptr % dimNames(2) = 'nCells'
      r2Ptr % defaultValue = 0.0
      allocate(r2Ptr % attLists(1))
      allocate(r2Ptr % attLists(1) % attList)
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'units', 'kg K m^{-3}')
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'long_name', 'rho*theta_m/zz perturbation from the reference state value')
      r2Ptr % missingValue = MPAS_REAL_FILLVAL
      r2Ptr % block => block

      if (met_stage_outActive) then
         r2Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'rtheta_p', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'rtheta_p', r2Ptr)

! Define variable rho_p
      allocate(r2Ptr)

! Setting up time level 1
      r2Ptr % fieldName = 'rho_p'
      r2Ptr % isVarArray = .false.
      r2Ptr % isDecomposed = .true.
      r2Ptr % hasTimeDimension = .true.
      r2Ptr % isPersistent = .true.
      r2Ptr % isActive = .false.
! Setting up dimensions
      r2Ptr % dimNames(1) = 'nVertLevels'
      r2Ptr % dimNames(2) = 'nCells'
      r2Ptr % defaultValue = 0.0
      allocate(r2Ptr % attLists(1))
      allocate(r2Ptr % attLists(1) % attList)
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'units', 'kg m^{-3}')
      call mpas_add_att(r2Ptr % attLists(1) % attList, 'long_name', 'rho/zz perturbation from the reference state value, advanced '&
//'over acoustic steps')
      r2Ptr % missingValue = MPAS_REAL_FILLVAL
      r2Ptr % block => block

      if (met_stage_outActive) then
         r2Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'rho_p', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'rho_p', r2Ptr)



      if (associated(newSubPool)) then
         call mpas_pool_add_config(newSubPool, 'on_a_sphere', block % domain % on_a_sphere)
         call mpas_pool_add_config(newSubPool, 'sphere_radius', block % domain % sphere_radius)
         call mpas_pool_add_config(newSubPool, 'is_periodic', block % domain % is_periodic)
         call mpas_pool_add_config(newSubPool, 'x_period', block % domain % x_period)
         call mpas_pool_add_config(newSubPool, 'y_period', block % domain % y_period)
      end if

   end subroutine init_atm_generate_pool_diag


   subroutine init_atm_generate_pool_diag_physics(block, structPool, dimensionPool, packagePool)
      use mpas_derived_types
      use mpas_pool_routines
      use mpas_io_units
      use mpas_io, only : MPAS_REAL_FILLVAL, MPAS_INT_FILLVAL, MPAS_CHAR_FILLVAL
      implicit none
      type (block_type), intent(inout), pointer :: block
      type (mpas_pool_type), intent(inout) :: structPool
      type (mpas_pool_type), intent(inout) :: dimensionPool
      type (mpas_pool_type), intent(in) :: packagePool

      type (field0DReal), pointer :: r0Ptr
      type (field1DReal), pointer :: r1Ptr
      type (field2DReal), pointer :: r2Ptr
      type (field3DReal), pointer :: r3Ptr
      type (field4DReal), pointer :: r4Ptr
      type (field5DReal), pointer :: r5Ptr
      type (field0DInteger), pointer :: i0Ptr
      type (field1DInteger), pointer :: i1Ptr
      type (field2DInteger), pointer :: i2Ptr
      type (field3DInteger), pointer :: i3Ptr
      type (field0DChar), pointer :: c0Ptr
      type (field1DChar), pointer :: c1Ptr
      type (field0DReal), dimension(:), pointer :: r0aPtr
      type (field1DReal), dimension(:), pointer :: r1aPtr
      type (field2DReal), dimension(:), pointer :: r2aPtr
      type (field3DReal), dimension(:), pointer :: r3aPtr
      type (field4DReal), dimension(:), pointer :: r4aPtr
      type (field5DReal), dimension(:), pointer :: r5aPtr
      type (field0DInteger), dimension(:), pointer :: i0aPtr
      type (field1DInteger), dimension(:), pointer :: i1aPtr
      type (field2DInteger), dimension(:), pointer :: i2aPtr
      type (field3DInteger), dimension(:), pointer :: i3aPtr
      type (field0DChar), dimension(:), pointer :: c0aPtr
      type (field1DChar), dimension(:), pointer :: c1aPtr

      type (mpas_pool_type), pointer :: newSubPool
      integer :: group_counter
      logical :: group_started
      integer :: group_start
      integer :: index_counter
      integer, pointer :: const_index

      logical, pointer :: initial_condsActive
      logical, pointer :: sfc_updateActive
      logical, pointer :: lbcsActive
      logical, pointer :: gwd_stage_inActive
      logical, pointer :: gwd_stage_outActive
      logical, pointer :: vertical_stage_inActive
      logical, pointer :: vertical_stage_outActive
      logical, pointer :: met_stage_inActive
      logical, pointer :: met_stage_outActive
      logical, pointer :: first_guess_fieldActive


      integer :: numConstituents

      nullify(newSubPool)
      group_counter = -1
      group_started = .false.
      group_start = -1
      call mpas_pool_get_package(packagePool, 'initial_condsActive', initial_condsActive)
      call mpas_pool_get_package(packagePool, 'sfc_updateActive', sfc_updateActive)
      call mpas_pool_get_package(packagePool, 'lbcsActive', lbcsActive)
      call mpas_pool_get_package(packagePool, 'gwd_stage_inActive', gwd_stage_inActive)
      call mpas_pool_get_package(packagePool, 'gwd_stage_outActive', gwd_stage_outActive)
      call mpas_pool_get_package(packagePool, 'vertical_stage_inActive', vertical_stage_inActive)
      call mpas_pool_get_package(packagePool, 'vertical_stage_outActive', vertical_stage_outActive)
      call mpas_pool_get_package(packagePool, 'met_stage_inActive', met_stage_inActive)
      call mpas_pool_get_package(packagePool, 'met_stage_outActive', met_stage_outActive)
      call mpas_pool_get_package(packagePool, 'first_guess_fieldActive', first_guess_fieldActive)

      allocate(newSubPool)
      call mpas_pool_create_pool(newSubPool)
      call mpas_pool_add_subpool(structPool, 'diag_physics', newSubPool)
      call mpas_pool_add_subpool(block % allStructs, 'diag_physics', newSubPool)

! Define variable precipw
      allocate(r1Ptr)

! Setting up time level 1
      r1Ptr % fieldName = 'precipw'
      r1Ptr % isVarArray = .false.
      r1Ptr % isDecomposed = .true.
      r1Ptr % hasTimeDimension = .true.
      r1Ptr % isPersistent = .true.
      r1Ptr % isActive = .false.
! Setting up dimensions
      r1Ptr % dimNames(1) = 'nCells'
      r1Ptr % defaultValue = 0.0
      allocate(r1Ptr % attLists(1))
      allocate(r1Ptr % attLists(1) % attList)
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'units', 'kg m^{-2}')
      call mpas_add_att(r1Ptr % attLists(1) % attList, 'long_name', 'precipitable water')
      r1Ptr % missingValue = MPAS_REAL_FILLVAL
      r1Ptr % block => block

      if (met_stage_outActive) then
         r1Ptr % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'precipw', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'precipw', r1Ptr)



      if (associated(newSubPool)) then
         call mpas_pool_add_config(newSubPool, 'on_a_sphere', block % domain % on_a_sphere)
         call mpas_pool_add_config(newSubPool, 'sphere_radius', block % domain % sphere_radius)
         call mpas_pool_add_config(newSubPool, 'is_periodic', block % domain % is_periodic)
         call mpas_pool_add_config(newSubPool, 'x_period', block % domain % x_period)
         call mpas_pool_add_config(newSubPool, 'y_period', block % domain % y_period)
      end if

   end subroutine init_atm_generate_pool_diag_physics


   subroutine init_atm_generate_structs(block, structPool, dimensionPool, packagePool)
      use mpas_derived_types
      use mpas_io_units
      implicit none
      type (block_type), pointer, intent(inout) :: block
      type (mpas_pool_type), intent(inout) :: structPool
      type (mpas_pool_type), intent(inout) :: dimensionPool
      type (mpas_pool_type), intent(in) :: packagePool

      call init_atm_generate_pool_mesh(block, structPool, dimensionPool, packagePool)

      call init_atm_generate_pool_state(block, structPool, dimensionPool, packagePool)

      call init_atm_generate_pool_lbc_state(block, structPool, dimensionPool, packagePool)

      call init_atm_generate_pool_fg(block, structPool, dimensionPool, packagePool)

      call init_atm_generate_pool_diag(block, structPool, dimensionPool, packagePool)

      call init_atm_generate_pool_diag_physics(block, structPool, dimensionPool, packagePool)

   end subroutine init_atm_generate_structs
