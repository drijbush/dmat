Program dummy_main

  Use mpi

  Use distributed_k_module
  Use proc_mapping_module
  Use matrix_mapping_module
  Use distributed_matrix_module
  
  Implicit None

  Integer, Dimension( : ), Allocatable :: i_hold
  Type( proc_mapping ), Dimension( : ), Allocatable :: split_proc_map
  Type( matrix_mapping ), Dimension( : ), Allocatable :: split_matrix_map
  Type( real_distributed_matrix ) :: base_matrix

  Integer :: i
  Integer :: error

  Type( real_distributed_matrix ) :: A

  Call mpi_init( error )

  Call distributed_matrix_init( MPI_COMM_WORLD, base_matrix )
  Call distributed_matrix_set_default_blocking( 3 )

!!$  Call matrix_mapping_base%print()
!!$  Call proc_mapping_base%split( [ 1, 2, 1, 2 ], 'k split', split_proc_map, i_hold )
!!$  Call matrix_mapping_base%split( [ 1, 2, 1, 2 ], 'k split', split_matrix_map, i_hold )
!!$  Do i = 1, Size( split_matrix_map )
!!$     Call split_matrix_map( i )%print()
!!$  End Do

  Call A%create( 20, 15, base_matrix )
  
  Call distributed_matrix_finalise
  
  Call mpi_finalize( error )
  
  
End Program dummy_main
