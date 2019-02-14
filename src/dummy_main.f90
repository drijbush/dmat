Program dummy_main

  Use mpi

  Use numbers_module, Only : wp
  Use distributed_k_module
  Use proc_mapping_module
  Use matrix_mapping_module
  Use distributed_matrix_module
  
  Implicit None
  
  Type( real_distributed_matrix ) :: base_matrix

  Integer :: n, nb
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

  Call test_matmul_real()
  Call test_matmul_real_ops() ! Internal compiler error in gcc 5.4
  Call test_matmul_complex()

  Call test_diag_real()
  Call test_diag_complex()
  
  Call distributed_matrix_finalise

  Call test_diag_k_real()
  Call test_diag_k_complex()
  Call test_diag_k_real_nm()
  Call test_diag_extract_complex()
  
  Call mpi_finalize( error )

Contains

  Subroutine test_diag_real()

    Type( real_distributed_matrix )              :: A
    Type( real_distributed_matrix ), Allocatable :: Q
    Real( wp ), Dimension( : )     , Allocatable :: E

    Type( real_distributed_matrix )              :: QT, B, C
    
    Real( wp ), Dimension( :, : ), Allocatable :: A_global
    Real( wp ), Dimension( :, : ), Allocatable :: tmp
    Real( wp ), Dimension( :    ), Allocatable :: work, ev

    Integer :: unit = 10
    Integer :: i

    Call create_global_real( n, A_global, A )
    A_global = A_global + Transpose( A_global )

    Call A%set_by_global( 1, n, 1, n, A_global )
    Call A%diag( Q, E )

    QT = Q
    QT = QT%dagger()
    B = QT%multiply( A )
    C = B%multiply( Q )
    Allocate( tmp( 1:n, 1:n ) )
    Call C%get_by_global( 1, n, 1, n, tmp )
    Do i = 1, n
       tmp( i, i ) = tmp( i, i ) - E( i )
    End Do
    
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
       Write( *, '( a, t64, g24.16 )' ) 'Diag  :Real    Case:             :Max absolute eval diff : ', Maxval( Abs( E - ev ) )
       Write( *, '( a, t64, g24.16 )' ) 'Diag  :Real    Case:             :Max sim tran trace test: ', Maxval( Abs( tmp ) )
    End If
  
  End Subroutine test_diag_real
  
  Subroutine test_diag_complex()

    Type( complex_distributed_matrix )              :: A
    Type( complex_distributed_matrix ), Allocatable :: Q
    Real( wp ), Dimension( : )        , Allocatable :: E
    
    Type( complex_distributed_matrix )              :: QT, B, C

    Complex( wp ), Dimension( :, : ), Allocatable :: A_global
    Complex( wp ), Dimension( :, : ), Allocatable :: tmp
    Complex( wp ), Dimension( :    ), Allocatable :: cwork

    Real( wp ), Dimension( : ), Allocatable :: rwork
    Real( wp ), Dimension( : ), Allocatable :: ev

    Integer :: unit = 10
    Integer :: i

    Call create_global_complex( n, A_global, A )
    
    A_global = A_global + Conjg( Transpose( A_global ) )
    
    Call A%set_by_global( 1, n, 1, n, A_global )
    Call A%diag( Q, E )
    
    QT = Q
    QT = QT%dagger()
    B = QT%multiply( A )
    C = B%multiply( Q )
    Allocate( tmp( 1:n, 1:n ) )
    Call C%get_by_global( 1, n, 1, n, tmp )
    Do i = 1, n
       tmp( i, i ) = tmp( i, i ) - E( i )
    End Do

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
       Write( *, '( a, t64, g24.16 )' ) 'Diag  :Complex Case:             :Max absolute eval diffe: ', Maxval( Abs( E - ev ) )
       Write( *, '( a, t64, g24.16 )' ) 'Diag  :Complex Case:             :Max sim tran trace test: ', Maxval( Abs( tmp ) )
    End If
    
  End Subroutine test_diag_complex

  Subroutine test_matmul_real()

    Type ( real_distributed_matrix ) :: A
    Type ( real_distributed_matrix ) :: B
    Type ( real_distributed_matrix ) :: C
    Type ( real_distributed_matrix ) :: D
    
    Real( wp ), Dimension( :, : ), Allocatable :: A_global
    Real( wp ), Dimension( :, : ), Allocatable :: B_global
    Real( wp ), Dimension( :, : ), Allocatable :: C_global
    Real( wp ), Dimension( :, : ), Allocatable :: D_global

    Call create_global_real( n, A_global, A )
    Call create_global_real( n, B_global, B )

    Allocate( D_global( 1:n, 1:n ) )

    ! Checks on Various Tranposes

    !NN
    C_global = Matmul( A_global, B_global )
    C = A%multiply( B )
    Call C%get_by_global( 1, n, 1, n, D_global )
    If( rank == 0 ) Then
       Write( *, '( a, t64, g24.16 )' ) 'Matmul:Real    Case:Transposes NN:Max absolute difference: ', &
            Maxval( Abs( C_global - D_global ) )
    End If

    !NT
    C_global = Matmul( A_global, Transpose( B_global ) )
    C = A%multiply( B%dagger() )
    Call C%get_by_global( 1, n, 1, n, D_global )
    If( rank == 0 ) Then
       Write( *, '( a, t64, g24.16 )' ) 'Matmul:Real    Case:Transposes NT:Max absolute difference: ', &
            Maxval( Abs( C_global - D_global ) )
    End If

    !TN
    C_global = Matmul( Transpose( A_global ), B_global )
    D = A%dagger()
    C = D%multiply( B )
    Call C%get_by_global( 1, n, 1, n, D_global )
    If( rank == 0 ) Then
       Write( *, '( a, t64, g24.16 )' ) 'Matmul:Real    Case:Transposes TN:Max absolute difference: ', &
            Maxval( Abs( C_global - D_global ) )
    End If

    !TT
    C_global = Matmul( Transpose( A_global ), Transpose( B_global ) )
    D = A%dagger()
    C = D%multiply( B%dagger() )
    Call C%get_by_global( 1, n, 1, n, D_global )
    If( rank == 0 ) Then
       Write( *, '( a, t64, g24.16 )' ) 'Matmul:Real    Case:Transposes TT:Max absolute difference: ', &
            Maxval( Abs( C_global - D_global ) )
    End If

  End Subroutine test_matmul_real
  
  ! Broken in gcc 5.4 - internal compiler error
  Subroutine test_matmul_real_ops()

    Type ( real_distributed_matrix ) :: A
    Type ( real_distributed_matrix ) :: B
    Type ( real_distributed_matrix ) :: C
    Type ( real_distributed_matrix ) :: D
    
    Real( wp ), Dimension( :, : ), Allocatable :: A_global
    Real( wp ), Dimension( :, : ), Allocatable :: B_global
    Real( wp ), Dimension( :, : ), Allocatable :: C_global
    Real( wp ), Dimension( :, : ), Allocatable :: D_global

    Call create_global_real( n, A_global, A )
    Call create_global_real( n, B_global, B )

    Allocate( D_global( 1:n, 1:n ) )

    ! Checks on Various Tranposes

    !NN
    C_global = Matmul( A_global, B_global )
    C = A * B
    Call C%get_by_global( 1, n, 1, n, D_global )
    If( rank == 0 ) Then
       Write( *, '( a, t64, g24.16 )' ) 'Matmul:Real    Case:Transposes NN:Max absolute difference: ', &
            Maxval( Abs( C_global - D_global ) )
    End If

    !NT
    C_global = Matmul( A_global, Transpose( B_global ) )
    C = A * .Dagger. B
    Call C%get_by_global( 1, n, 1, n, D_global )
    If( rank == 0 ) Then
       Write( *, '( a, t64, g24.16 )' ) 'Matmul:Real    Case:Transposes NT:Max absolute difference: ', &
            Maxval( Abs( C_global - D_global ) )
    End If

    !TN
    C_global = Matmul( Transpose( A_global ), B_global )
    D = .Dagger. A
    C = D * B
    Call C%get_by_global( 1, n, 1, n, D_global )
    If( rank == 0 ) Then
       Write( *, '( a, t64, g24.16 )' ) 'Matmul:Real    Case:Transposes TN:Max absolute difference: ', &
            Maxval( Abs( C_global - D_global ) )
    End If

    !TT
    C_global = Matmul( Transpose( A_global ), Transpose( B_global ) )
    D = .Dagger. A
    C = D * ( .Dagger. B )
    Call C%get_by_global( 1, n, 1, n, D_global )
    If( rank == 0 ) Then
       Write( *, '( a, t64, g24.16 )' ) 'Matmul:Real    Case:Transposes TT:Max absolute difference: ', &
            Maxval( Abs( C_global - D_global ) )
    End If

  End Subroutine test_matmul_real_ops
  
  Subroutine test_matmul_complex()

    Type ( complex_distributed_matrix ) :: A
    Type ( complex_distributed_matrix ) :: B
    Type ( complex_distributed_matrix ) :: C
    Type ( complex_distributed_matrix ) :: D
    
    Complex( wp ), Dimension( :, : ), Allocatable :: A_global
    Complex( wp ), Dimension( :, : ), Allocatable :: B_global
    Complex( wp ), Dimension( :, : ), Allocatable :: C_global
    Complex( wp ), Dimension( :, : ), Allocatable :: D_global

    Call create_global_complex( n, A_global, A )
    Call create_global_complex( n, B_global, B )

    Allocate( D_global( 1:n, 1:n ) )

    ! Checks on Various Daggers

    !NN
    C_global = Matmul( A_global, B_global )
    C = A%multiply( B )
    Call C%get_by_global( 1, n, 1, n, D_global )
    If( rank == 0 ) Then
       Write( *, '( a, t64, g24.16 )' ) 'Matmul:Complex Case:Daggers    NN:Max absolute difference: ', &
            Maxval( Abs( C_global - D_global ) )
    End If

    !NC
    C_global = Matmul( A_global, Conjg( Transpose( B_global ) ) )
    C = A%multiply( B%dagger() )
    Call C%get_by_global( 1, n, 1, n, D_global )
    If( rank == 0 ) Then
       Write( *, '( a, t64, g24.16 )' ) 'Matmul:Complex Case:Daggers    NC:Max absolute difference: ', &
            Maxval( Abs( C_global - D_global ) )
    End If

    !CN
    C_global = Matmul( Conjg( Transpose( A_global ) ), B_global )
    D = A%dagger()
    C = D%multiply( B )
    Call C%get_by_global( 1, n, 1, n, D_global )
    If( rank == 0 ) Then
       Write( *, '( a, t64, g24.16 )' ) 'Matmul:Complex Case:Daggers    CN:Max absolute difference: ', &
            Maxval( Abs( C_global - D_global ) )
    End If

    !CC
    C_global = Matmul( Conjg( Transpose( A_global ) ), Conjg( Transpose( B_global ) ) )
    D = A%dagger()
    C = D%multiply( B%dagger() )
    Call C%get_by_global( 1, n, 1, n, D_global )
    If( rank == 0 ) Then
       Write( *, '( a, t64, g24.16 )' ) 'Matmul:Complex Case:Daggers    CC:Max absolute difference: ', &
            Maxval( Abs( C_global - D_global ) )
    End If
  End Subroutine test_matmul_complex
  
  Subroutine create_global_real( n, A_g, A )

    Integer                                     , Intent( In    ) :: n
    Real( wp ), Dimension( :, :   ), Allocatable, Intent(   Out ) :: A_g
    Type( real_distributed_matrix )             , Intent(   Out ) :: A

    Allocate( A_g( 1:n, 1:n ) )
    Call Random_number( A_g )
    
    Call A%create( n, n, base_matrix )
    Call A%set_by_global( 1, n, 1, n, A_g )
    
  End Subroutine create_global_real

  Subroutine create_global_complex( n, A_g, A )

    Integer                                        , Intent( In    ) :: n
    Complex( wp ), Dimension( :, :   ), Allocatable, Intent(   Out ) :: A_g
    Type( complex_distributed_matrix )             , Intent(   Out ) :: A

    Real( wp ), Dimension( :, : ), Allocatable :: tmp

    Allocate( tmp( 1:n, 1:n ) )
    
    Allocate( A_g( 1:n, 1:n ) )
    Call Random_number( tmp )
    A_g = tmp
    Call Random_number( tmp )
    A_g = A_g + Cmplx( 0.0_wp, tmp, wp )
    
    Call A%create( n, n, base_matrix )
    Call A%set_by_global( 1, n, 1, n, A_g )
    
  End Subroutine create_global_complex

  Subroutine test_diag_k_real()

    Type( distributed_k_matrix ) :: A
    Type( distributed_k_matrix ) :: Q
    Real( wp ), Dimension( : )  , Allocatable :: E

    Type( distributed_k_matrix ) :: QT, B, C
    Type( distributed_k_matrix ) :: base_k
    
    Real( wp ), Dimension( :, : ), Allocatable :: A_global
    Real( wp ), Dimension( :, : ), Allocatable :: tmp
    Real( wp ), Dimension( :    ), Allocatable :: work, ev

    Integer :: unit = 10
    Integer :: i

    Allocate( A_global( 1:n, 1:n ) )
    Allocate( tmp( 1:n, 1:n ) )

    Call random_number( A_global )
    A_global = A_global + Transpose( A_global )
    
    Call distributed_k_matrix_init( MPI_COMM_WORLD, base_k )
    
    Call A%create( .False., 1, [ 0, 0, 0 ], n, n, base_k )
    Call A%set_by_global( 1, n, 1, n, A_global )
    Call A%diag( Q, E )

    QT = .Dagger. Q
    B = QT * A
    C = B  * Q
    Call C%get_by_global( 1, n, 1, n, tmp )
    Do i = 1, n
       tmp( i, i ) = tmp( i, i ) - E( i )
    End Do
    Call distributed_k_matrix_finalise
    
    Allocate( work( 1:64 * n ) )
    Allocate( ev( 1:n ) )
    Call dsyev( 'v', 'l', n, a_global, n, ev, work, Size( work ), error )

    If( rank == 0 ) Then
       Open( unit, file = 'real_eval_k_diff.dat' )
       Do i = 1, n
          Write( unit, '( i4, 2( g30.16, 1x ), g24.16 )' ) i, E( i ), ev( i ), E( i ) - ev( i )
       End Do
       Write( unit, '( a, 1x, g24.16 )' ) 'Max absolute difference: ', Maxval( Abs( E - ev ) )
       Close( unit )
       Write( *, '( a, t64, g24.16 )' ) 'Diag  :Real    Case:             :Max absolute eval diff : ', Maxval( Abs( E - ev ) )
       Write( *, '( a, t64, g24.16 )' ) 'Diag  :Real    Case:             :Max sim tran trace test: ', Maxval( Abs( tmp ) )
    End If
  
  End Subroutine test_diag_k_real

  Subroutine test_diag_k_complex()

    Type( distributed_k_matrix )              :: A
    Type( distributed_k_matrix )              :: Q
    Real( wp ), Dimension( : )  , Allocatable :: E
    
    Type( distributed_k_matrix ) :: QT, B, C
    Type( distributed_k_matrix ) :: base_k

    Complex( wp ), Dimension( :, : ), Allocatable :: A_global
    Complex( wp ), Dimension( :, : ), Allocatable :: tmp
    Complex( wp ), Dimension( :    ), Allocatable :: cwork

    Real( wp ), Dimension( :, : ), Allocatable :: rtmp

    Real( wp ), Dimension( : ), Allocatable :: rwork
    Real( wp ), Dimension( : ), Allocatable :: ev

    Integer :: unit = 10
    Integer :: i

    Allocate( A_global( 1:n, 1:n ) )
    Allocate( tmp( 1:n, 1:n ) )

    Allocate( rtmp( 1:n, 1:n ) )
    Call random_number( rtmp )
    A_global = rtmp
    Call random_number( rtmp )
    A_global = A_global + Cmplx( 0.0_wp, rtmp, wp )
    A_global = A_global + Conjg( Transpose( A_global ) )
    
    Call distributed_k_matrix_init( MPI_COMM_WORLD, base_k )
    
    Call A%create( .True., 1, [ 0, 0, 0 ], n, n, base_k )
    Call A%set_by_global( 1, n, 1, n, A_global )
    Call A%diag( Q, E )
    
    QT = Q
    QT = .Dagger. Q
    B = QT * A
    C = B  * Q
    Call C%get_by_global( 1, n, 1, n, tmp )
    Do i = 1, n
       tmp( i, i ) = tmp( i, i ) - E( i )
    End Do
    Call distributed_k_matrix_finalise

    Allocate( cwork( 1:64 * n ) )
    Allocate( rwork( 1:3 * n - 2 ) )
    Allocate( ev( 1:n ) )
    Call zheev( 'v', 'l', n, a_global, n, ev, cwork, Size( cwork ), rwork, error )

    If( rank == 0 ) Then
       Open( unit, file = 'complex_eval_k_diff.dat' )
       Do i = 1, n
          Write( unit, '( i4, 2( g30.16, 1x ), g24.16 )' ) i, E( i ), ev( i ), E( i ) - ev( i )
       End Do
       Write( unit, '( a, 1x, g24.16 )' ) 'Max absolute difference: ', Maxval( Abs( E - ev ) )
       Close( unit )
       Write( *, '( a, t64, g24.16 )' ) 'Diag  :Complex Case:             :Max absolute eval diffe: ', Maxval( Abs( E - ev ) )
       Write( *, '( a, t64, g24.16 )' ) 'Diag  :Complex Case:             :Max sim tran trace test: ', Maxval( Abs( tmp ) )
    End If
    
  End Subroutine test_diag_k_complex

  Subroutine test_diag_k_real_nm()

    Type( distributed_k_matrix ) :: A_nm
    Type( distributed_k_matrix ) :: A_nn
    Type( distributed_k_matrix ) :: Q
    Real( wp ), Dimension( : )  , Allocatable :: E

    Type( distributed_k_matrix ) :: QT, B, C
    Type( distributed_k_matrix ) :: base_k
    
    Real( wp ), Dimension( :, : ), Allocatable :: A_global_nm
    Real( wp ), Dimension( :, : ), Allocatable :: A_global
    Real( wp ), Dimension( :, : ), Allocatable :: tmp
    Real( wp ), Dimension( :    ), Allocatable :: work, ev

    Integer :: unit = 10
    Integer :: i
    Integer :: m

    m = Nint( n * 0.75_wp )

    Allocate( A_global_nm( 1:n, 1:m ) )
    Allocate( A_global( 1:n, 1:n ) )
    Allocate( tmp( 1:n, 1:n ) )

    A_global_nm = 0.01_wp
    Do i = 1, m
       A_global_nm( i, i ) = 1.0_wp
    End Do
    A_global = Matmul( A_global_nm, Transpose( A_global_nm ) )
    
    Call distributed_k_matrix_init( MPI_COMM_WORLD, base_k )
    
    Call A_nm%create( .False., 1, [ 0, 0, 0 ], n, m, base_k )
    Call A_nm%set_by_global( 1, n, 1, m, A_global_nm )
    A_nn = A_nm * .Dagger. A_nm
    Call A_nn%diag( Q, E )

    QT = .Dagger. Q
    B = QT * A_nn
    C = B  * Q
    Call C%get_by_global( 1, n, 1, n, tmp )
    Do i = 1, n
       tmp( i, i ) = tmp( i, i ) - E( i )
    End Do
    Call distributed_k_matrix_finalise
    
    Allocate( work( 1:64 * n ) )
    Allocate( ev( 1:n ) )
    Call dsyev( 'v', 'l', n, a_global, n, ev, work, Size( work ), error )

    If( rank == 0 ) Then
       Open( unit, file = 'real_eval_k_nm_diff.dat' )
       Do i = 1, n
          Write( unit, '( i4, 2( g30.16, 1x ), g24.16 )' ) i, E( i ), ev( i ), E( i ) - ev( i )
       End Do
       Write( unit, '( a, 1x, g24.16 )' ) 'Max absolute difference: ', Maxval( Abs( E - ev ) )
       Close( unit )
       Write( *, '( a, t64, g24.16 )' ) 'Diag  :Real    Case:             :Max absolute eval diff : ', Maxval( Abs( E - ev ) )
       Write( *, '( a, t64, g24.16 )' ) 'Diag  :Real    Case:             :Max sim tran trace test: ', Maxval( Abs( tmp ) )
    End If
  
  End Subroutine test_diag_k_real_nm

  Subroutine test_diag_extract_complex()

    Type( distributed_k_matrix )              :: A
    Type( distributed_k_matrix )              :: Q
    Real( wp ), Dimension( : )  , Allocatable :: E
    
    Type( distributed_k_matrix ) :: Qe, QeT, B, C
    Type( distributed_k_matrix ) :: base_k

    Complex( wp ), Dimension( :, : ), Allocatable :: A_global
    Complex( wp ), Dimension( :, : ), Allocatable :: tmp
    Complex( wp ), Dimension( :    ), Allocatable :: cwork

    Real( wp ), Dimension( :, : ), Allocatable :: rtmp

    Real( wp ), Dimension( : ), Allocatable :: rwork
    Real( wp ), Dimension( : ), Allocatable :: ev

    Integer :: m
    Integer :: unit = 10
    Integer :: i

    m = Nint( n * 0.75_wp )
    
    Allocate( A_global( 1:n, 1:n ) )
    Allocate( tmp( 1:m, 1:m ) )

    Allocate( rtmp( 1:n, 1:n ) )
    Call random_number( rtmp )
    A_global = rtmp
    Call random_number( rtmp )
    A_global = A_global + Cmplx( 0.0_wp, rtmp, wp )
    A_global = A_global + Conjg( Transpose( A_global ) )
    
    Call distributed_k_matrix_init( MPI_COMM_WORLD, base_k )
    
    Call A%create( .True., 1, [ 0, 0, 0 ], n, n, base_k )
    Call A%set_by_global( 1, n, 1, n, A_global )
    Call A%diag( Q, E )

    Call Q%extract_cols( 1, m, Qe )
    QeT = .Dagger. Qe
    B = QeT * A
    C = B  * Qe
    Call C%get_by_global( 1, m, 1, m, tmp )
    Do i = 1, m
       tmp( i, i ) = tmp( i, i ) - E( i )
    End Do
    Call distributed_k_matrix_finalise

    Allocate( cwork( 1:64 * n ) )
    Allocate( rwork( 1:3 * n - 2 ) )
    Allocate( ev( 1:n ) )
    Call zheev( 'v', 'l', n, a_global, n, ev, cwork, Size( cwork ), rwork, error )

    If( rank == 0 ) Then
       Open( unit, file = 'complex_eval_extract_diff.dat' )
       Do i = 1, m
          Write( unit, '( i4, 2( g30.16, 1x ), g24.16 )' ) i, E( i ), ev( i ), E( i ) - ev( i )
       End Do
       Write( unit, '( a, 1x, g24.16 )' ) 'Max absolute difference: ', Maxval( Abs( E - ev ) )
       Close( unit )
       Write( *, '( a, t64, g24.16 )' ) 'Diag  :Complex Case:             :Max absolute eval diffe: ', Maxval( Abs( E - ev ) )
       Write( *, '( a, t64, g24.16 )' ) 'Diag  :Complex Case:             :Max sim tran trace test: ', Maxval( Abs( tmp ) )
    End If
    
  End Subroutine test_diag_extract_complex

End Program dummy_main
