Program dummy_main

  Use mpi

  Use distributed_k_module
  Use proc_mapping_module
  Use matrix_mapping_module
  Use distributed_matrix_module
  
  Implicit None

!!$  Integer, Dimension( : ), Allocatable :: i_hold
!!$  Type( proc_mapping ), Dimension( : ), Allocatable :: split_proc_map
!!$  Type( matrix_mapping ), Dimension( : ), Allocatable :: split_matrix_map
  Type( real_distributed_matrix ) :: base_matrix

  Real( wp ), Dimension( : ), Allocatable :: E

  Integer :: i
  Integer :: error

  Type( real_distributed_matrix ) :: A
  Class( real_distributed_matrix ), Allocatable :: Q
  Integer, Parameter :: n = 8
  Real( wp ), Dimension( 1:n, 1:n ) :: A_global
  Real( wp ), Dimension( : ), Allocatable :: work, ev
  Integer :: rank
  
  Call mpi_init( error )
  Call mpi_comm_rank( mpi_comm_world, rank, error )

  Call distributed_matrix_set_default_blocking( 3 )
  Call distributed_matrix_init( MPI_COMM_WORLD, base_matrix )

!!$  Call matrix_mapping_base%print()
!!$  Call proc_mapping_base%split( [ 1, 2, 1, 2 ], 'k split', split_proc_map, i_hold )
!!$  Call matrix_mapping_base%split( [ 1, 2, 1, 2 ], 'k split', split_matrix_map, i_hold )
!!$  Do i = 1, Size( split_matrix_map )
!!$     Call split_matrix_map( i )%print()
!!$  End Do

!!$  A_global = 0.0_wp
!!$  Do i = 1, Size( A_global, Dim = 1 )
!!$     A_global( i, i ) = Real( i, wp )
!!$  End Do
  Call Random_number( A_global )
  A_global = A_global + Transpose( A_global )
  
  Call A%create( n, n, base_matrix )
  Call A%set_by_global( 1, n, 1, n, A_global )
  Call A%diag( Q, E )

  Allocate( work( 1:64 * n ) )
  Allocate( ev( 1:n ) )
  Call dsyev( 'v', 'l', n, a_global, n, ev, work, Size( work ), error )

  If( rank == 0 ) Then
     Do i = 1, n
        Write( *, '( i4, 2( g30.16, 1x ), g24.16 )' ) i, E( i ), ev( i ), E( i ) - ev( i )
     End Do
     Write( *, '( a, 1x, g24.16 )' ) 'Max absolute difference: ', Maxval( Abs( E - ev ) )
  End If
  
  Call distributed_matrix_finalise
  
  Call mpi_finalize( error )
  
  
End Program dummy_main
