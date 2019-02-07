Program dummy_main

  Use mpi

  Use distributed_k_module
  Use proc_mapping_module
  Use matrix_mapping_module
  Use distributed_matrix_module
  
  Implicit None

  Type( real_distributed_matrix ) :: base_matrix

  Integer :: n, nb
  Integer :: i
  Integer :: error

  Integer :: rank
  
  Call mpi_init( error )
  Call mpi_comm_rank( mpi_comm_world, rank, error )

  If( rank == 0 ) Then
     Write( *, * ) 'n = ?'
     Read ( *, * ) n
     Write( *, * ) 'nb = ?'
     Read ( *, * ) nb
  End If
  Call mpi_bcast( n , 1, MPI_INTEGER, 0, MPI_COMM_WORLD, error )
  Call mpi_bcast( nb, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, error )
  
  Call distributed_matrix_set_default_blocking( 3 )
  Call distributed_matrix_init( MPI_COMM_WORLD, base_matrix )

  Call test_real()
  Call test_complex()
  
  Call distributed_matrix_finalise
  
  Call mpi_finalize( error )

Contains

  Subroutine test_real()

    Type ( real_distributed_matrix )              :: A
    Class( real_distributed_matrix ), Allocatable :: Q
    Real( wp ), Dimension( : )      , Allocatable :: E
    
    Real( wp ), Dimension( :, : ), Allocatable :: A_global
    Real( wp ), Dimension( :    ), Allocatable :: work, ev

    Integer :: unit = 10

    Allocate( A_global( 1:n, 1:n ) )
    
    Call Random_number( A_global )
    A_global = A_global + Transpose( A_global )
    
    Call A%create( n, n, base_matrix )
    Call A%set_by_global( 1, n, 1, n, A_global )
    Call A%diag( Q, E )
    
    Allocate( work( 1:64 * n ) )
    Allocate( ev( 1:n ) )
    Call dsyev( 'v', 'l', n, a_global, n, ev, work, Size( work ), error )

    If( rank == 0 ) Then
       Open( unit, file = 'real_eval_diff.dat' )
       Do i = 1, n
          Write( unit, '( i4, 2( g30.16, 1x ), g24.16 )' ) i, E( i ), ev( i ), E( i ) - ev( i )
       End Do
       Write( unit, '( a, 1x, g24.16 )' ) 'Max absolute difference: ', Maxval( Abs( E - ev ) )
       Close( unit )
       Write( *, '( a, 1x, g24.16 )' ) 'Real    Case: Max absolute difference: ', Maxval( Abs( E - ev ) )
    End If
  
  End Subroutine test_real
  
  Subroutine test_complex()

    Type ( complex_distributed_matrix )              :: A
    Class( complex_distributed_matrix ), Allocatable :: Q
    Real( wp ), Dimension( : )         , Allocatable :: E
    
    Complex( wp ), Dimension( :, : ), Allocatable :: A_global
    Complex( wp ), Dimension( :    ), Allocatable :: cwork

    Real( wp ), Dimension( :, : ), Allocatable :: tmp1, tmp2

    Real( wp ), Dimension( : ), Allocatable :: rwork
    Real( wp ), Dimension( : ), Allocatable :: ev

    Integer :: unit = 10

    Allocate( A_global( 1:n, 1:n ) )
    Allocate( tmp1( 1:n, 1:n ) )
    Allocate( tmp2( 1:n, 1:n ) )
    
    Call Random_number( tmp1 )
    Call Random_number( tmp2 )
    A_global = Cmplx( tmp1, tmp2, wp )
    A_global = A_global + Conjg( Transpose( A_global ) )
    
    Call A%create( n, n, base_matrix )
    Call A%set_by_global( 1, n, 1, n, A_global )
    Call A%diag( Q, E )
    
    Allocate( cwork( 1:64 * n ) )
    Allocate( rwork( 1:3 * n - 2 ) )
    Allocate( ev( 1:n ) )
    Call zheev( 'v', 'l', n, a_global, n, ev, cwork, Size( cwork ), rwork, error )

    If( rank == 0 ) Then
       Open( unit, file = 'complex_eval_diff.dat' )
       Do i = 1, n
          Write( unit, '( i4, 2( g30.16, 1x ), g24.16 )' ) i, E( i ), ev( i ), E( i ) - ev( i )
       End Do
       Write( unit, '( a, 1x, g24.16 )' ) 'Max absolute difference: ', Maxval( Abs( E - ev ) )
       Close( unit )
       Write( *, '( a, 1x, g24.16 )' ) 'Complex Case: Max absolute difference: ', Maxval( Abs( E - ev ) )
    End If
    
  End Subroutine test_complex
  
  
End Program dummy_main
