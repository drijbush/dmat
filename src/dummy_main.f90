Program dummy_main

  Use mpi

  Use k_point_matrix_module
  Use mapping_module
  
  Implicit None

  Integer, Dimension( : ), Allocatable :: i_hold
  Class( mapping ), Dimension( : ), Allocatable :: split_map
  Integer :: error

  Call mpi_init( error )
  Call mapping_init( MPI_COMM_WORLD )
  Call mapping_base_map%print()
  Call mapping_base_map%split( [ 1, 2, 1, 2 ], 'k split', split_map, i_hold )
  Call mapping_finalise
  Call mpi_finalize( error )
  
End Program dummy_main
