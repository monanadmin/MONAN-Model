function init_atm_setup_immutable_streams(manager) result(iErr)

   use MPAS_derived_types, only : MPAS_streamManager_type, &
                                  MPAS_STREAM_INPUT_OUTPUT, MPAS_STREAM_INPUT, &
                                  MPAS_STREAM_OUTPUT, MPAS_STREAM_NONE, MPAS_STREAM_PROPERTY_IMMUTABLE
   use MPAS_stream_manager, only : MPAS_stream_mgr_create_stream, MPAS_stream_mgr_set_property, &
                                   MPAS_stream_mgr_add_field, MPAS_stream_mgr_add_pool
   use mpas_io_units

   implicit none

   type (MPAS_streamManager_type), pointer :: manager
   character (len=StrKIND) :: packages
   integer :: iErr

   iErr = 0

   call MPAS_stream_mgr_create_stream(manager, 'input', MPAS_STREAM_INPUT, 'x1.40962.grid.nc', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'latCell', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'lonCell', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'xCell', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'yCell', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'zCell', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'indexToCellID', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'latEdge', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'lonEdge', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'xEdge', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'yEdge', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'zEdge', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'indexToEdgeID', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'latVertex', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'lonVertex', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'xVertex', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'yVertex', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'zVertex', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'indexToVertexID', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'cellsOnEdge', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'nEdgesOnCell', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'nEdgesOnEdge', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'edgesOnCell', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'edgesOnEdge', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'weightsOnEdge', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'dvEdge', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'dcEdge', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'angleEdge', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'areaCell', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'areaTriangle', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'cellsOnCell', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'verticesOnCell', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'verticesOnEdge', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'edgesOnVertex', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'cellsOnVertex', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'kiteAreasOnVertex', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'meshDensity', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'nominalMinDc', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'bdyMaskCell', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'bdyMaskEdge', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'bdyMaskVertex', ierr=ierr)
   write(packages,'(a)') 'gwd_stage_in;vertical_stage_in;met_stage_in'
   call MPAS_stream_mgr_add_field(manager, 'input', 'edgeNormalVectors', packages=packages, ierr=ierr)
   write(packages,'(a)') 'gwd_stage_in;vertical_stage_in;met_stage_in'
   call MPAS_stream_mgr_add_field(manager, 'input', 'localVerticalUnitVectors', packages=packages, ierr=ierr)
   write(packages,'(a)') 'gwd_stage_in;vertical_stage_in;met_stage_in'
   call MPAS_stream_mgr_add_field(manager, 'input', 'cellTangentPlane', packages=packages, ierr=ierr)
   write(packages,'(a)') 'gwd_stage_in;vertical_stage_in;met_stage_in'
   call MPAS_stream_mgr_add_field(manager, 'input', 'fEdge', packages=packages, ierr=ierr)
   write(packages,'(a)') 'gwd_stage_in;vertical_stage_in;met_stage_in'
   call MPAS_stream_mgr_add_field(manager, 'input', 'fVertex', packages=packages, ierr=ierr)
   write(packages,'(a)') 'gwd_stage_in;vertical_stage_in;met_stage_in'
   call MPAS_stream_mgr_add_field(manager, 'input', 'ter', packages=packages, ierr=ierr)
   write(packages,'(a)') 'gwd_stage_in;vertical_stage_in;met_stage_in'
   call MPAS_stream_mgr_add_field(manager, 'input', 'landmask', packages=packages, ierr=ierr)
   write(packages,'(a)') 'gwd_stage_in;vertical_stage_in;met_stage_in'
   call MPAS_stream_mgr_add_field(manager, 'input', 'mminlu', packages=packages, ierr=ierr)
   write(packages,'(a)') 'gwd_stage_in;vertical_stage_in;met_stage_in'
   call MPAS_stream_mgr_add_field(manager, 'input', 'isice_lu', packages=packages, ierr=ierr)
   write(packages,'(a)') 'gwd_stage_in;vertical_stage_in;met_stage_in'
   call MPAS_stream_mgr_add_field(manager, 'input', 'iswater_lu', packages=packages, ierr=ierr)
   write(packages,'(a)') 'gwd_stage_in;vertical_stage_in;met_stage_in'
   call MPAS_stream_mgr_add_field(manager, 'input', 'ivgtyp', packages=packages, ierr=ierr)
   write(packages,'(a)') 'gwd_stage_in;vertical_stage_in;met_stage_in'
   call MPAS_stream_mgr_add_field(manager, 'input', 'isltyp', packages=packages, ierr=ierr)
   write(packages,'(a)') 'gwd_stage_in;vertical_stage_in;met_stage_in'
   call MPAS_stream_mgr_add_field(manager, 'input', 'snoalb', packages=packages, ierr=ierr)
   write(packages,'(a)') 'gwd_stage_in;vertical_stage_in;met_stage_in'
   call MPAS_stream_mgr_add_field(manager, 'input', 'soiltemp', packages=packages, ierr=ierr)
   write(packages,'(a)') 'gwd_stage_in;vertical_stage_in;met_stage_in'
   call MPAS_stream_mgr_add_field(manager, 'input', 'greenfrac', packages=packages, ierr=ierr)
   write(packages,'(a)') 'gwd_stage_in;vertical_stage_in;met_stage_in'
   call MPAS_stream_mgr_add_field(manager, 'input', 'shdmin', packages=packages, ierr=ierr)
   write(packages,'(a)') 'gwd_stage_in;vertical_stage_in;met_stage_in'
   call MPAS_stream_mgr_add_field(manager, 'input', 'shdmax', packages=packages, ierr=ierr)
   write(packages,'(a)') 'gwd_stage_in;vertical_stage_in;met_stage_in'
   call MPAS_stream_mgr_add_field(manager, 'input', 'albedo12m', packages=packages, ierr=ierr)
   write(packages,'(a)') 'vertical_stage_in;met_stage_in'
   call MPAS_stream_mgr_add_field(manager, 'input', 'var2d', packages=packages, ierr=ierr)
   write(packages,'(a)') 'vertical_stage_in;met_stage_in'
   call MPAS_stream_mgr_add_field(manager, 'input', 'con', packages=packages, ierr=ierr)
   write(packages,'(a)') 'vertical_stage_in;met_stage_in'
   call MPAS_stream_mgr_add_field(manager, 'input', 'oa1', packages=packages, ierr=ierr)
   write(packages,'(a)') 'vertical_stage_in;met_stage_in'
   call MPAS_stream_mgr_add_field(manager, 'input', 'oa2', packages=packages, ierr=ierr)
   write(packages,'(a)') 'vertical_stage_in;met_stage_in'
   call MPAS_stream_mgr_add_field(manager, 'input', 'oa3', packages=packages, ierr=ierr)
   write(packages,'(a)') 'vertical_stage_in;met_stage_in'
   call MPAS_stream_mgr_add_field(manager, 'input', 'oa4', packages=packages, ierr=ierr)
   write(packages,'(a)') 'vertical_stage_in;met_stage_in'
   call MPAS_stream_mgr_add_field(manager, 'input', 'ol1', packages=packages, ierr=ierr)
   write(packages,'(a)') 'vertical_stage_in;met_stage_in'
   call MPAS_stream_mgr_add_field(manager, 'input', 'ol2', packages=packages, ierr=ierr)
   write(packages,'(a)') 'vertical_stage_in;met_stage_in'
   call MPAS_stream_mgr_add_field(manager, 'input', 'ol3', packages=packages, ierr=ierr)
   write(packages,'(a)') 'vertical_stage_in;met_stage_in'
   call MPAS_stream_mgr_add_field(manager, 'input', 'ol4', packages=packages, ierr=ierr)
   write(packages,'(a)') 'gwd_stage_in;vertical_stage_in;met_stage_in'
   call MPAS_stream_mgr_add_field(manager, 'input', 'deriv_two', packages=packages, ierr=ierr)
   write(packages,'(a)') 'gwd_stage_in;vertical_stage_in;met_stage_in'
   call MPAS_stream_mgr_add_field(manager, 'input', 'defc_a', packages=packages, ierr=ierr)
   write(packages,'(a)') 'gwd_stage_in;vertical_stage_in;met_stage_in'
   call MPAS_stream_mgr_add_field(manager, 'input', 'defc_b', packages=packages, ierr=ierr)
   write(packages,'(a)') 'gwd_stage_in;vertical_stage_in;met_stage_in'
   call MPAS_stream_mgr_add_field(manager, 'input', 'cell_gradient_coef_x', packages=packages, ierr=ierr)
   write(packages,'(a)') 'gwd_stage_in;vertical_stage_in;met_stage_in'
   call MPAS_stream_mgr_add_field(manager, 'input', 'cell_gradient_coef_y', packages=packages, ierr=ierr)
   write(packages,'(a)') 'gwd_stage_in;vertical_stage_in;met_stage_in'
   call MPAS_stream_mgr_add_field(manager, 'input', 'coeffs_reconstruct', packages=packages, ierr=ierr)
   write(packages,'(a)') 'met_stage_in'
   call MPAS_stream_mgr_add_field(manager, 'input', 'cf1', packages=packages, ierr=ierr)
   write(packages,'(a)') 'met_stage_in'
   call MPAS_stream_mgr_add_field(manager, 'input', 'cf2', packages=packages, ierr=ierr)
   write(packages,'(a)') 'met_stage_in'
   call MPAS_stream_mgr_add_field(manager, 'input', 'cf3', packages=packages, ierr=ierr)
   write(packages,'(a)') 'met_stage_in'
   call MPAS_stream_mgr_add_field(manager, 'input', 'zgrid', packages=packages, ierr=ierr)
   write(packages,'(a)') 'met_stage_in'
   call MPAS_stream_mgr_add_field(manager, 'input', 'rdzw', packages=packages, ierr=ierr)
   write(packages,'(a)') 'met_stage_in'
   call MPAS_stream_mgr_add_field(manager, 'input', 'dzu', packages=packages, ierr=ierr)
   write(packages,'(a)') 'met_stage_in'
   call MPAS_stream_mgr_add_field(manager, 'input', 'rdzu', packages=packages, ierr=ierr)
   write(packages,'(a)') 'met_stage_in'
   call MPAS_stream_mgr_add_field(manager, 'input', 'fzm', packages=packages, ierr=ierr)
   write(packages,'(a)') 'met_stage_in'
   call MPAS_stream_mgr_add_field(manager, 'input', 'fzp', packages=packages, ierr=ierr)
   write(packages,'(a)') 'met_stage_in'
   call MPAS_stream_mgr_add_field(manager, 'input', 'zxu', packages=packages, ierr=ierr)
   write(packages,'(a)') 'met_stage_in'
   call MPAS_stream_mgr_add_field(manager, 'input', 'zz', packages=packages, ierr=ierr)
   write(packages,'(a)') 'met_stage_in'
   call MPAS_stream_mgr_add_field(manager, 'input', 'zb', packages=packages, ierr=ierr)
   write(packages,'(a)') 'met_stage_in'
   call MPAS_stream_mgr_add_field(manager, 'input', 'zb3', packages=packages, ierr=ierr)
   write(packages,'(a)') 'met_stage_in'
   call MPAS_stream_mgr_add_field(manager, 'input', 'dss', packages=packages, ierr=ierr)
   call MPAS_stream_mgr_set_property(manager, 'input', MPAS_STREAM_PROPERTY_IMMUTABLE, .true., ierr=ierr)

   call MPAS_stream_mgr_create_stream(manager, 'output', MPAS_STREAM_OUTPUT, 'x1.40962.init.nc', ierr=ierr)
   write(packages,'(a)') 'met_stage_out'
   call MPAS_stream_mgr_add_field(manager, 'output', 'scalars', packages=packages, ierr=ierr)
   write(packages,'(a)') 'met_stage_out'
   call MPAS_stream_mgr_add_field(manager, 'output', 'initial_time', packages=packages, ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'output', 'xtime', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'output', 'latCell', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'output', 'lonCell', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'output', 'xCell', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'output', 'yCell', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'output', 'zCell', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'output', 'indexToCellID', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'output', 'latEdge', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'output', 'lonEdge', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'output', 'xEdge', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'output', 'yEdge', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'output', 'zEdge', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'output', 'indexToEdgeID', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'output', 'latVertex', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'output', 'lonVertex', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'output', 'xVertex', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'output', 'yVertex', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'output', 'zVertex', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'output', 'indexToVertexID', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'output', 'cellsOnEdge', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'output', 'nEdgesOnCell', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'output', 'nEdgesOnEdge', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'output', 'edgesOnCell', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'output', 'edgesOnEdge', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'output', 'weightsOnEdge', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'output', 'dvEdge', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'output', 'dcEdge', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'output', 'angleEdge', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'output', 'areaCell', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'output', 'areaTriangle', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'output', 'cellsOnCell', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'output', 'verticesOnCell', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'output', 'verticesOnEdge', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'output', 'edgesOnVertex', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'output', 'cellsOnVertex', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'output', 'kiteAreasOnVertex', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'output', 'meshDensity', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'output', 'nominalMinDc', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'output', 'bdyMaskCell', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'output', 'bdyMaskEdge', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'output', 'bdyMaskVertex', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'output', 'edgeNormalVectors', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'output', 'localVerticalUnitVectors', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'output', 'cellTangentPlane', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'output', 'fEdge', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'output', 'fVertex', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'output', 'ter', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'output', 'landmask', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'output', 'mminlu', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'output', 'isice_lu', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'output', 'iswater_lu', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'output', 'ivgtyp', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'output', 'isltyp', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'output', 'snoalb', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'output', 'soiltemp', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'output', 'greenfrac', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'output', 'shdmin', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'output', 'shdmax', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'output', 'albedo12m', ierr=ierr)
   write(packages,'(a)') 'gwd_stage_out;vertical_stage_out;met_stage_out'
   call MPAS_stream_mgr_add_field(manager, 'output', 'var2d', packages=packages, ierr=ierr)
   write(packages,'(a)') 'gwd_stage_out;vertical_stage_out;met_stage_out'
   call MPAS_stream_mgr_add_field(manager, 'output', 'con', packages=packages, ierr=ierr)
   write(packages,'(a)') 'gwd_stage_out;vertical_stage_out;met_stage_out'
   call MPAS_stream_mgr_add_field(manager, 'output', 'oa1', packages=packages, ierr=ierr)
   write(packages,'(a)') 'gwd_stage_out;vertical_stage_out;met_stage_out'
   call MPAS_stream_mgr_add_field(manager, 'output', 'oa2', packages=packages, ierr=ierr)
   write(packages,'(a)') 'gwd_stage_out;vertical_stage_out;met_stage_out'
   call MPAS_stream_mgr_add_field(manager, 'output', 'oa3', packages=packages, ierr=ierr)
   write(packages,'(a)') 'gwd_stage_out;vertical_stage_out;met_stage_out'
   call MPAS_stream_mgr_add_field(manager, 'output', 'oa4', packages=packages, ierr=ierr)
   write(packages,'(a)') 'gwd_stage_out;vertical_stage_out;met_stage_out'
   call MPAS_stream_mgr_add_field(manager, 'output', 'ol1', packages=packages, ierr=ierr)
   write(packages,'(a)') 'gwd_stage_out;vertical_stage_out;met_stage_out'
   call MPAS_stream_mgr_add_field(manager, 'output', 'ol2', packages=packages, ierr=ierr)
   write(packages,'(a)') 'gwd_stage_out;vertical_stage_out;met_stage_out'
   call MPAS_stream_mgr_add_field(manager, 'output', 'ol3', packages=packages, ierr=ierr)
   write(packages,'(a)') 'gwd_stage_out;vertical_stage_out;met_stage_out'
   call MPAS_stream_mgr_add_field(manager, 'output', 'ol4', packages=packages, ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'output', 'deriv_two', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'output', 'defc_a', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'output', 'defc_b', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'output', 'cell_gradient_coef_x', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'output', 'cell_gradient_coef_y', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'output', 'coeffs_reconstruct', ierr=ierr)
   write(packages,'(a)') 'vertical_stage_out;met_stage_out'
   call MPAS_stream_mgr_add_field(manager, 'output', 'cf1', packages=packages, ierr=ierr)
   write(packages,'(a)') 'vertical_stage_out;met_stage_out'
   call MPAS_stream_mgr_add_field(manager, 'output', 'cf2', packages=packages, ierr=ierr)
   write(packages,'(a)') 'vertical_stage_out;met_stage_out'
   call MPAS_stream_mgr_add_field(manager, 'output', 'cf3', packages=packages, ierr=ierr)
   write(packages,'(a)') 'vertical_stage_out;met_stage_out'
   call MPAS_stream_mgr_add_field(manager, 'output', 'zgrid', packages=packages, ierr=ierr)
   write(packages,'(a)') 'vertical_stage_out;met_stage_out'
   call MPAS_stream_mgr_add_field(manager, 'output', 'rdzw', packages=packages, ierr=ierr)
   write(packages,'(a)') 'vertical_stage_out;met_stage_out'
   call MPAS_stream_mgr_add_field(manager, 'output', 'dzu', packages=packages, ierr=ierr)
   write(packages,'(a)') 'vertical_stage_out;met_stage_out'
   call MPAS_stream_mgr_add_field(manager, 'output', 'rdzu', packages=packages, ierr=ierr)
   write(packages,'(a)') 'vertical_stage_out;met_stage_out'
   call MPAS_stream_mgr_add_field(manager, 'output', 'fzm', packages=packages, ierr=ierr)
   write(packages,'(a)') 'vertical_stage_out;met_stage_out'
   call MPAS_stream_mgr_add_field(manager, 'output', 'fzp', packages=packages, ierr=ierr)
   write(packages,'(a)') 'vertical_stage_out;met_stage_out'
   call MPAS_stream_mgr_add_field(manager, 'output', 'zxu', packages=packages, ierr=ierr)
   write(packages,'(a)') 'vertical_stage_out;met_stage_out'
   call MPAS_stream_mgr_add_field(manager, 'output', 'zz', packages=packages, ierr=ierr)
   write(packages,'(a)') 'vertical_stage_out;met_stage_out'
   call MPAS_stream_mgr_add_field(manager, 'output', 'zb', packages=packages, ierr=ierr)
   write(packages,'(a)') 'vertical_stage_out;met_stage_out'
   call MPAS_stream_mgr_add_field(manager, 'output', 'zb3', packages=packages, ierr=ierr)
   write(packages,'(a)') 'vertical_stage_out;met_stage_out'
   call MPAS_stream_mgr_add_field(manager, 'output', 'dss', packages=packages, ierr=ierr)
   write(packages,'(a)') 'met_stage_out'
   call MPAS_stream_mgr_add_field(manager, 'output', 'u_init', packages=packages, ierr=ierr)
   write(packages,'(a)') 'met_stage_out'
   call MPAS_stream_mgr_add_field(manager, 'output', 'v_init', packages=packages, ierr=ierr)
   write(packages,'(a)') 'met_stage_out'
   call MPAS_stream_mgr_add_field(manager, 'output', 't_init', packages=packages, ierr=ierr)
   write(packages,'(a)') 'met_stage_out'
   call MPAS_stream_mgr_add_field(manager, 'output', 'qv_init', packages=packages, ierr=ierr)
   write(packages,'(a)') 'met_stage_out'
   call MPAS_stream_mgr_add_field(manager, 'output', 'u', packages=packages, ierr=ierr)
   write(packages,'(a)') 'met_stage_out'
   call MPAS_stream_mgr_add_field(manager, 'output', 'w', packages=packages, ierr=ierr)
   write(packages,'(a)') 'met_stage_out'
   call MPAS_stream_mgr_add_field(manager, 'output', 'dz', packages=packages, ierr=ierr)
   write(packages,'(a)') 'met_stage_out'
   call MPAS_stream_mgr_add_field(manager, 'output', 'dzs', packages=packages, ierr=ierr)
   write(packages,'(a)') 'met_stage_out'
   call MPAS_stream_mgr_add_field(manager, 'output', 'zs', packages=packages, ierr=ierr)
   write(packages,'(a)') 'met_stage_out'
   call MPAS_stream_mgr_add_field(manager, 'output', 'sh2o', packages=packages, ierr=ierr)
   write(packages,'(a)') 'met_stage_out'
   call MPAS_stream_mgr_add_field(manager, 'output', 'smois', packages=packages, ierr=ierr)
   write(packages,'(a)') 'met_stage_out'
   call MPAS_stream_mgr_add_field(manager, 'output', 'tslb', packages=packages, ierr=ierr)
   write(packages,'(a)') 'met_stage_out'
   call MPAS_stream_mgr_add_field(manager, 'output', 'tmn', packages=packages, ierr=ierr)
   write(packages,'(a)') 'met_stage_out'
   call MPAS_stream_mgr_add_field(manager, 'output', 'skintemp', packages=packages, ierr=ierr)
   write(packages,'(a)') 'met_stage_out'
   call MPAS_stream_mgr_add_field(manager, 'output', 'sst', packages=packages, ierr=ierr)
   write(packages,'(a)') 'met_stage_out'
   call MPAS_stream_mgr_add_field(manager, 'output', 'snow', packages=packages, ierr=ierr)
   write(packages,'(a)') 'met_stage_out'
   call MPAS_stream_mgr_add_field(manager, 'output', 'snowc', packages=packages, ierr=ierr)
   write(packages,'(a)') 'met_stage_out'
   call MPAS_stream_mgr_add_field(manager, 'output', 'snowh', packages=packages, ierr=ierr)
   write(packages,'(a)') 'met_stage_out'
   call MPAS_stream_mgr_add_field(manager, 'output', 'xice', packages=packages, ierr=ierr)
   write(packages,'(a)') 'met_stage_out'
   call MPAS_stream_mgr_add_field(manager, 'output', 'seaice', packages=packages, ierr=ierr)
   write(packages,'(a)') 'met_stage_out'
   call MPAS_stream_mgr_add_field(manager, 'output', 'vegfra', packages=packages, ierr=ierr)
   write(packages,'(a)') 'met_stage_out'
   call MPAS_stream_mgr_add_field(manager, 'output', 'sfc_albbck', packages=packages, ierr=ierr)
   write(packages,'(a)') 'met_stage_out'
   call MPAS_stream_mgr_add_field(manager, 'output', 'xland', packages=packages, ierr=ierr)
   write(packages,'(a)') 'met_stage_out'
   call MPAS_stream_mgr_add_field(manager, 'output', 'rho', packages=packages, ierr=ierr)
   write(packages,'(a)') 'met_stage_out'
   call MPAS_stream_mgr_add_field(manager, 'output', 'theta', packages=packages, ierr=ierr)
   write(packages,'(a)') 'met_stage_out'
   call MPAS_stream_mgr_add_field(manager, 'output', 'relhum', packages=packages, ierr=ierr)
   write(packages,'(a)') 'met_stage_out'
   call MPAS_stream_mgr_add_field(manager, 'output', 'rho_base', packages=packages, ierr=ierr)
   write(packages,'(a)') 'met_stage_out'
   call MPAS_stream_mgr_add_field(manager, 'output', 'theta_base', packages=packages, ierr=ierr)
   write(packages,'(a)') 'met_stage_out'
   call MPAS_stream_mgr_add_field(manager, 'output', 'surface_pressure', packages=packages, ierr=ierr)
   write(packages,'(a)') 'met_stage_out'
   call MPAS_stream_mgr_add_field(manager, 'output', 'precipw', packages=packages, ierr=ierr)
   write(packages,'(a)') 'met_stage_out'
   call MPAS_stream_mgr_add_field(manager, 'output', 'h_oml_initial', packages=packages, ierr=ierr)
   write(packages,'(a)') 'met_stage_out'
   call MPAS_stream_mgr_add_field(manager, 'output', 'u10', packages=packages, ierr=ierr)
   write(packages,'(a)') 'met_stage_out'
   call MPAS_stream_mgr_add_field(manager, 'output', 'v10', packages=packages, ierr=ierr)
   write(packages,'(a)') 'met_stage_out'
   call MPAS_stream_mgr_add_field(manager, 'output', 'q2', packages=packages, ierr=ierr)
   write(packages,'(a)') 'met_stage_out'
   call MPAS_stream_mgr_add_field(manager, 'output', 'rh2', packages=packages, ierr=ierr)
   write(packages,'(a)') 'met_stage_out'
   call MPAS_stream_mgr_add_field(manager, 'output', 't2m', packages=packages, ierr=ierr)
   call MPAS_stream_mgr_set_property(manager, 'output', MPAS_STREAM_PROPERTY_IMMUTABLE, .true., ierr=ierr)

   call MPAS_stream_mgr_create_stream(manager, 'surface', MPAS_STREAM_OUTPUT, 'x1.40962.sfc_update.nc', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'surface', 'xtime', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'surface', 'sst', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'surface', 'xice', ierr=ierr)
   call MPAS_stream_mgr_set_property(manager, 'surface', MPAS_STREAM_PROPERTY_IMMUTABLE, .true., ierr=ierr)

   call MPAS_stream_mgr_create_stream(manager, 'lbc', MPAS_STREAM_OUTPUT, 'lbc.$Y-$M-$D_$h.$m.$s.nc', ierr=ierr)
   call MPAS_stream_mgr_add_pool(manager, 'lbc', 'lbc_state', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'lbc', 'xtime', ierr=ierr)
   call MPAS_stream_mgr_set_property(manager, 'lbc', MPAS_STREAM_PROPERTY_IMMUTABLE, .true., ierr=ierr)

end function init_atm_setup_immutable_streams


