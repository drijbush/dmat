Module init_matrix_module

  Use numbers_module      , Only : wp

  Implicit None

  Interface init_sym_matrix
     Procedure init_sym_matrix_complex
  End Interface init_sym_matrix

  Interface init_pos_def_matrix
     Procedure init_sym_matrix_complex
  End Interface init_pos_def_matrix

Contains

  Subroutine init_sym_matrix_complex( A )
    Complex( wp ), Dimension( :, : ), Intent( Out ) :: A
    Real( wp ), Dimension( :, : ), Allocatable :: tmp
    Integer :: n
    Integer :: i
    n = Size( A, Dim = 1 )
    Allocate( tmp( 1:n, 1:n  ) )
    Call Random_number( tmp )
    A = tmp
    Call Random_number( tmp )
    A = A + Cmplx( 0.0_wp, tmp, Kind = wp )
    A = A - 0.5_wp
    A = 0.5_wp * ( A + Conjg( Transpose( A ) ) )
    A = A / 15.0_wp
    Do i = 1, Size( A, Dim = 1 )
       A( i, i ) = 1.0_wp
    End Do
  End Subroutine init_sym_matrix_complex

  Subroutine init_pos_def_matrix_complex( A )
    Complex( wp ), Dimension( :, : ), Intent( Out ) :: A
    Call init_sym_matrix( A )
    A = Matmul( Transpose( Conjg( A ) ), A )
  End Subroutine init_pos_def_matrix_complex

End Module init_matrix_module


Program main

  Use mpi

  Use numbers_module      , Only : wp
  Use init_matrix_module  , Only : init_sym_matrix, init_pos_def_matrix

  Implicit None

  Complex( wp ), Dimension( :, : ), Allocatable :: A, S, P, P_dist, Q, S_copy

  Complex( wp ), Dimension( : ), Allocatable :: cwork

  Real( wp ), Dimension( : ), Allocatable :: Evals_lapack
  Real( wp ), Dimension( : ), Allocatable :: Evals_dist
  Real( wp ), Dimension( : ), Allocatable :: rwork
  
  Real( wp ) :: Energy_lapack, Energy_dist

  Integer, Dimension( : ), Allocatable :: iwork
  
  Integer :: n, nb, ne
  Integer :: rank
  Integer :: error

  Call mpi_init( error )
  Call mpi_comm_rank( mpi_comm_world, rank, error )

  If( rank == 0 ) Then
     Write( *, * ) 'n = ?'
     Read ( *, * ) n
     Write( *, * ) 'ne = ?'
     Read ( *, * ) ne
     Write( *, * ) 'nb = ?'
     Read ( *, * ) nb
  End If
  Call mpi_bcast( n , 1, MPI_INTEGER, 0, MPI_COMM_WORLD, error )
  Call mpi_bcast( ne, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, error )
  Call mpi_bcast( nb, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, error )

  Allocate( A( 1:n, 1:n ) )
  Allocate( S( 1:n, 1:n ) )
  Allocate( Evals_lapack( 1:n ) )

  Call init_sym_matrix( A )
  Call init_pos_def_matrix( S )

  Q = A
  S_copy = S
  Allocate( cwork( 1:2 * n + n * n ) )
  Allocate( rwork( 1:1 + 6 * n + 2 * n * n ) )
  Allocate( iwork(  1:3 + 5 * n ) )
  Call zhegvd( 1, 'V', 'U', n, Q, n, S_copy, n, evals_lapack, &
       cwork, Size( cwork ), rwork, Size( rwork ), iwork, Size( iwork ), error )
  P = Matmul( Q( :, 1:ne ), Transpose( Conjg( Q( :, 1:ne ) ) ) )
  If( error /= 0 ) Then
     If( rank == 0 ) Then
        Write( *, * ) 'zhegvd failed with error = ', error
     End If
     Call mpi_finalize( error )
     Stop
  End If

  Energy_lapack = Sum( evals_lapack( 1:ne ) )
  If( rank == 0 ) Then
     Write( *, '( a, g26.16 )' ) 'Lapack    Energy: ', Energy_lapack
  End If

  Call distributed( ne, A, S, P_dist, evals_dist )
  Energy_dist = Sum( evals_dist( 1:ne ) )
  If( rank == 0 ) Then
     Write( *, '( a, g26.16 )' ) 'Scalapack Energy: ', Energy_dist
     Write( *, '( a, g26.16 )' ) 'Max diff in evals ', Maxval( Abs( evals_dist - evals_lapack ) )
     Write( *, '( a, g26.16 )' ) 'Max diff in P     ', Maxval( Abs( P_dist     - P ) )
  End If
  
  Call distributed_by_diag( ne, A, S, P_dist, evals_dist )
  Energy_dist = Sum( evals_dist( 1:ne ) )
  If( rank == 0 ) Then
     Write( *, '( a, g26.16 )' ) 'Scalapack Energy: ', Energy_dist
     Write( *, '( a, g26.16 )' ) 'Max diff in evals ', Maxval( evals_dist - evals_lapack( 1:Size( evals_dist ) ) )
     If( Size( P_dist ) == Size( P ) ) Then
        Write( *, '( a, g26.16 )' ) 'Max diff in P     ', Maxval( Abs( P_dist     - P ) )
     End If
  End If

  Call mpi_finalize( error )

Contains

  Subroutine distributed( ne, A_rep, S_rep, P_rep, evals )

    Use distributed_k_module, Only : distributed_k_matrix, &
         distributed_k_matrix_init, distributed_k_matrix_finalise

    Integer                         ,              Intent( In    ) :: ne
    Complex( wp ), Dimension( :, : ),              Intent( In    ) :: A_rep
    Complex( wp ), Dimension( :, : ),              Intent( In    ) :: S_rep
    Complex( wp ), Dimension( :, : ), Allocatable, Intent(   Out ) :: P_rep
    Real( wp )   , Dimension( :    ), Allocatable, Intent(   Out ) :: evals
    
    Type( distributed_k_matrix ) :: A
    Type( distributed_k_matrix ) :: S
    Type( distributed_k_matrix ) :: L, L_inv, L_inv_T
    Type( distributed_k_matrix ) :: Q, Qe
    Type( distributed_k_matrix ) :: P
    Type( distributed_k_matrix ) :: Identity
    Type( distributed_k_matrix ) :: base_k

    Integer :: n

    n = Size( A_rep, Dim = 1 )

    Call distributed_k_matrix_init( MPI_COMM_WORLD, base_k )

    ! Put A, S into distributed matrices
    Call A%create( .True., 1, [ 0, 0, 0 ], n, n, base_k )
    Call A%set_by_global( 1, n, 1, n, A_rep )
    Call S%create( .True., 1, [ 0, 0, 0 ], n, n, base_k )
    Call S%set_by_global( 1, n, 1, n, S_rep )

    ! Generate the orthogonalising matrix
    Call Identity%create( .True., 1, [ 0, 0, 0 ], n, n, base_k )
    Call Identity%set_to_identity()
    L     = S%Choleski() 
    L_inv = L%solve( Identity )
    L_inv_T = .Dagger. L_inv

    ! Transform AO -> MO
    A = L_inv * A * L_inv_T

    ! Diagonalise in the MO basis
    Call A%diag( Q, evals )

    ! Transformt the evecs MO -> AO
    Q = L_inv_T * Q

    ! Get the ne occupied eigenstates
    Qe = Q%extract( 1, n, 1, ne )

    ! Form the Density matrix
    P = Qe * .Dagger. Qe

    ! Re-replicate into standard array
    Allocate( P_rep( 1:n, 1:n ) )
    Call P%get_by_global( 1, n, 1, n, P_rep )

    Call distributed_k_matrix_finalise()
    
  End Subroutine distributed

  Subroutine distributed_by_diag( ne, A_rep, S_rep, P_rep, evals )

    Use distributed_k_module, Only : distributed_k_matrix, &
         distributed_k_matrix_init, distributed_k_matrix_finalise
    
    Integer                         ,              Intent( In    ) :: ne
    Complex( wp ), Dimension( :, : ),              Intent( In    ) :: A_rep
    Complex( wp ), Dimension( :, : ),              Intent( In    ) :: S_rep
    Complex( wp ), Dimension( :, : ), Allocatable, Intent(   Out ) :: P_rep
    Real   ( wp ), Dimension( :    ), Allocatable, Intent(   Out ) :: evals
    
    Type( distributed_k_matrix ) :: A
    Type( distributed_k_matrix ) :: S
    Type( distributed_k_matrix ) :: Q, Qe, QeT, Qproj
    Type( distributed_k_matrix ) :: P
    Type( distributed_k_matrix ) :: base_k

    Integer :: n, nb
    Integer :: start

    Real( wp ), Parameter :: tol = 0.2_wp

    n = Size( A_rep, Dim = 1 )

    Call distributed_k_matrix_init( MPI_COMM_WORLD, base_k )

    ! Put A, S into distributed matrices
    Call A%create( .True., 1, [ 0, 0, 0 ], n, n, base_k )
    Call A%set_by_global( 1, n, 1, n, A_rep )
    Call S%create( .True., 1, [ 0, 0, 0 ], n, n, base_k )
    Call S%set_by_global( 1, n, 1, n, S_rep )

    ! Generate orthogonalising transform
    Call S%diag( Q, evals )
    Do start = 1, n
       If( evals( start ) > tol ) Then
          Exit
       End If
    End Do
    If( rank == 0 ) Then
       Write( *, * ) 'Ignoring ', start - 1, ' evecs of the S matrix'
       Write( *, '( a, g24.16 )' ) 'Smallest eval: ', evals( 1 )
       If( start - 1 > 0 ) Then
          Write( *, '( a, g24.16 )' ) 'Biggest ignored: ', evals( start - 1 )
       End If
    End If
    Qe = Q%extract( 1, n, start, n )
    Qe = Qe * ( 1.0_wp / Sqrt( evals( start:n ) ) )
    QeT = .Dagger. Qe

    ! Transform AO -> projected MO
    A = QeT * A * Qe

    ! Diagonalise in the projected MO basis
    Call A%diag( Qproj, evals )

    ! Transformt the evecs MO -> AO
    Q = Qe * Qproj

    ! Get the ne occupied eigenstates
    nb = n - start + 1
    Qe = Q%extract( 1, nb, 1, ne )

    ! Form the Density matrix
    P = Qe * .Dagger. Qe

    ! Re-replicate into standard array
    Allocate( P_rep( 1:nb, 1:nb ) )
    Call P%get_by_global( 1, nb, 1, nb, P_rep )

    Call distributed_k_matrix_finalise()
    
  End Subroutine distributed_by_diag
  
End Program main
