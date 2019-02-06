Program dummy_main

  Use mpi

!!$  Use k_point_matrix_module
  Use proc_mapping_module
  Use matrix_mapping_module
  Use distributed_matrix_module
  
  Implicit None

  Integer, Dimension( : ), Allocatable :: i_hold
  Type( proc_mapping ), Dimension( : ), Allocatable :: split_proc_map
  Type( matrix_mapping ), Dimension( : ), Allocatable :: split_matrix_map
  Integer :: i
  Integer :: error

  Call mpi_init( error )
  Call matrix_mapping_init( MPI_COMM_WORLD )
  Call matrix_mapping_base%print()
  Call proc_mapping_base%split( [ 1, 2, 1, 2 ], 'k split', split_proc_map, i_hold )
  Call matrix_mapping_base%split( [ 1, 2, 1, 2 ], 'k split', split_matrix_map, i_hold )
  Do i = 1, Size( split_matrix_map )
     Call split_matrix_map( i )%print()
  End Do
  Call matrix_mapping_finalise
  Call mpi_finalize( error )
  
End Program dummy_main
