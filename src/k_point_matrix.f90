Module k_point_matrix_module

  Use numbers_module      , Only : wp
  Use distributed_k_module, Only : distributed_k_matrix_init, &
       distributed_k_matrix_finalise, distributed_k_matrix

  Implicit None

  Integer, Public, Parameter :: K_POINT_REAL = 0
  Integer, Public, Parameter :: K_POINT_COMPLEX = 1

  Integer, Private, Parameter :: INVALID = -1
  Integer, Private, Parameter :: NOT_ME  = -2

  Type k_point_info
     Integer                              :: spin
     Integer, Dimension( : ), Allocatable :: k_indices
  End type k_point_info

  Type, Private :: k_point_matrices
     ! External label - e.g. irrep number
     Integer                      :: label = INVALID
     Type( distributed_k_matrix ) :: matrix
  End type k_point_matrices

  Type, Private :: k_point
     Type( k_point_info     )                              :: info
     ! This splitting allows irreps
     Type( k_point_matrices ), Dimension( : ), Allocatable :: data
     ! Want to hide eventually
     Integer                                               :: communicator = INVALID
  End type k_point

  Type, Public :: ks_array
     Type( k_point_info ), Dimension( : ), Allocatable :: all_k_point_info
     ! this splitting allows multiple k points on this process
     Type( k_point      ), Dimension( : ), Allocatable :: my_k_points
     ! Want to hide eventually
     Integer                                           :: parent_communicator = INVALID
   Contains
     Procedure :: create => matrices_create
     Procedure :: diag   => matrices_diag
     Procedure, Private :: get_all_ks_index
     Procedure, Private :: get_my_ks_index
  End type ks_array
  
  Type, Public :: eval_storage
     Real( wp ), Dimension( : ), Allocatable :: evals
  End type eval_storage

  Private

  Public :: k_matrices_init
  Public :: k_matrices_finalise
  
Contains

  Subroutine k_matrices_init( comm, base_matrix )

    ! Want to think about what kind of thing is returned here ...
    
    Integer                        , Intent( In    ) :: comm
    Type   ( distributed_k_matrix ), Intent(   Out ) :: base_matrix

    Call  distributed_k_matrix_init( comm, base_matrix ) 

  End Subroutine k_matrices_init

  Subroutine k_matrices_finalise

    Call distributed_k_matrix_finalise
    
  End Subroutine k_matrices_finalise

  Subroutine matrices_create( matrices, n_spin, k_point_type, k_points, m, n, source )

    ! Create a matrix in all k point mode with no irreps //ism
    
    ! Also want to think about what kind of object should be used as a source
    
    Class( ks_array )              , Intent(   Out ) :: matrices
    Integer                        , Intent( In    ) :: n_spin
    Integer, Dimension( :    )     , Intent( In    ) :: k_point_type
    Integer, Dimension( :, : )     , Intent( In    ) :: k_points
    Integer                        , Intent( In    ) :: m
    Integer                        , Intent( In    ) :: n
    Type ( distributed_k_matrix   ), Intent( In    ) :: source

    Integer :: n_k_points
    Integer :: s, k, ks
    
    n_k_points = Size( k_point_type)

    ! Set up the all k point data structure
    Allocate( matrices%all_k_point_info( 1:n_spin * n_k_points ) )
    ks = 0
    Do s = 1, n_spin
       Do k = 1, n_k_points
          ks = ks + 1
          matrices%all_k_point_info( ks )%spin      = s
          matrices%all_k_point_info( ks )%k_indices = k_points( :, k )
       End Do
    End Do

    ! Now my k points
    Allocate( matrices%my_k_points( 1:n_spin * n_k_points ) )
    ks = 0
    Do s = 1, n_spin
       Do k = 1, n_k_points
          ks = ks + 1
          matrices%my_k_points( ks )%info%spin    = s
          matrices%my_k_points( ks )%info%k_indices = k_points( :, k )
          Allocate( matrices%my_k_points( ks )%data( 1:1 ) )
          matrices%my_k_points( ks )%data( 1 )%label = 1
          Call matrices%my_k_points( ks )%data( 1 )%matrix%create( k_point_type( k ) == K_POINT_COMPLEX, &
               s, k_points( :, k ), m, n, source )
          matrices%my_k_points( ks )%communicator = source%get_comm()
       End Do
    End Do

    matrices%parent_communicator = source%get_comm()
    
  End Subroutine matrices_create

  Subroutine matrices_diag( A, Q, E )

    Use mpi

    Class( ks_array     ),                 Intent( In    ) :: A
    Type ( ks_array     ),                 Intent(   Out ) :: Q
    Type ( eval_storage ), Dimension( : ), Intent(   Out ) :: E

    Real( wp ) :: rdum

    Integer, Dimension( 1:2 ) :: buff_send, buff_recv
    
    Integer :: me, me_parent
    Integer :: ks_root, nb
    Integer :: error
    Integer :: request
    Integer :: rsize, handle
    Integer :: my_ks, ks

    Logical :: sending_data
    
    ! Make Q have the same set up as A
    Q = A

    Do my_ks = 1, Size( A%my_k_points )
       ks = A%get_all_ks_index( my_ks )
       Call A%my_k_points( my_ks )%data( 1 )%matrix%diag( Q%my_k_points( my_ks )%data( 1 )%matrix, E( ks )%evals )
    End Do

    ! Replicate evals
    Call mpi_comm_rank( A%parent_communicator, me_parent, error )
    Do ks = 1, Size( A%all_k_point_info )
       my_ks = A%get_my_ks_index( ks )
       ! Work out who holds this set of evals and send how many there are
       ! the root node of the communicator holding them back to the root node of the parent communicator
       sending_data = .False.
       If( my_ks /= NOT_ME ) Then
          Call mpi_comm_rank( A%my_k_points( my_ks )%communicator, me, error )
          sending_data = me == 0
          If( sending_data ) then
             buff_send( 1 ) = me
             buff_send( 2 ) = Size( E( ks )%evals )
             Call mpi_isend( buff_send, 2, MPI_INTEGER, 0, 10, A%parent_communicator, request, error )
          End If
       End If
       If( me_parent == 0 ) Then
          Call mpi_recv( buff_recv, 2, MPI_INTEGER, MPI_ANY_SOURCE, 10, A%parent_communicator, MPI_STATUS_IGNORE, error )
       End If
       If( sending_data ) Then
          Call mpi_wait( request, MPI_STATUS_IGNORE, error )
       End If
       ! Now on root of parent bcast back to all
       Call mpi_bcast( buff_recv, 2, MPI_INTEGER, 0, A%parent_communicator, error )
       ks_root = buff_recv( 1 )
       nb      = buff_recv( 2 )
       ! Now know how many evals we will recv - allocate memory if haven't done so already
       If( .Not. Allocated( E( ks )%evals ) ) Then
          Allocate( E( ks )%evals( 1:nb ) )
       End If
       ! And finally bcast out the values from the root node for this set of evals
       Call mpi_sizeof( rdum, rsize, error )
       Call mpi_type_match_size( MPI_TYPECLASS_REAL, rsize, handle, error )
       Call mpi_bcast( E( ks )%evals, nb, handle, ks_root, A%parent_communicator, error )          
    End Do
    
  End Subroutine matrices_diag

  Pure Function get_all_ks_index( A, my_ks ) Result( ks )

    Integer :: ks

    Class( ks_array ), Intent( In ) :: A
    Integer          , Intent( In ) :: my_ks

    Do ks = 1, Size( A%all_k_point_info )
       If( A%all_k_point_info( ks )%spin == A%my_k_points( my_ks )%info%spin ) Then
          If( All( A%all_k_point_info( ks )%k_indices == A%my_k_points( my_ks )%info%k_indices ) ) Then
             Exit
          End If
       End If
    End Do
    
  End Function get_all_ks_index

  Pure Function get_my_ks_index( A, ks ) Result( my_ks )

    Integer :: my_ks

    Class( ks_array ), Intent( In ) :: A
    Integer          , Intent( In ) :: ks

    Integer :: k
    
    my_ks = NOT_ME

    Do k = 1, Size( A%my_k_points )
       If( A%all_k_point_info( ks )%spin == A%my_k_points( k )%info%spin ) Then
          If( All( A%all_k_point_info( ks )%k_indices == A%my_k_points( k )%info%k_indices ) ) Then
             my_ks = k
             Exit
          End If
       End If
    End Do

  End Function get_my_ks_index

End Module k_point_matrix_module
