Program dummy_main

  Use mpi

  Use k_point_matrix_module
  Use proc_mapping_module
  Use matrix_mapping_module
  
  Implicit None

  Integer, Dimension( : ), Allocatable :: i_hold
  Type( proc_mapping ), Dimension( : ), Allocatable :: split_map
  Integer :: error

  Call mpi_init( error )
  Call proc_mapping_init( MPI_COMM_WORLD )
  Call proc_mapping_base%print()
  Call proc_mapping_base%split( [ 1, 2, 1, 2 ], 'k split', split_map, i_hold )
  Call proc_mapping_finalise
  Call mpi_finalize( error )
  
End Program dummy_main
