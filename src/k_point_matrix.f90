Module k_point_matrix_module

  Use numbers_module      , Only : wp
  Use distributed_k_module, Only : distributed_k_matrix_init, &
       distributed_k_matrix_finalise, distributed_k_matrix

  Implicit None

  Integer, Public, Parameter :: K_POINT_REAL = 0
  Integer, Public, Parameter :: K_POINT_COMPLEX = 1

  Integer, Private, Parameter :: INVALID = -1

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
!!$     Integer                 , Dimension( : ), Allocatable :: spin
!!$     Integer                 , Dimension( : ), Allocatable :: k
     Type( k_point_info     )                              :: info
     ! This splitting allows irreps
     Type( k_point_matrices ), Dimension( : ), Allocatable :: data
     ! Want to hide eventually
     Integer                                               :: communicator = INVALID
  End type k_point

  Type, Public :: distributed_k_matrices
     Type( k_point_info ), Dimension( : ), Allocatable :: all_k_point_info
     ! this splitting allows multiple k points on this process
     Type( k_point      ), Dimension( : ), Allocatable :: my_k_points
     ! Want to hide eventually
     Integer                                           :: parent_communicator = INVALID
   Contains
     Procedure :: create => matrices_create
     Procedure :: diag   => matrices_diag
     Procedure, Private :: get_sk_index
  End type distributed_k_matrices
  
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
    
    Class( distributed_k_matrices ), Intent(   Out ) :: matrices
    Integer                        , Intent( In    ) :: n_spin
    Integer, Dimension( :    )     , Intent( In    ) :: k_point_type
    Integer, Dimension( :, : )     , Intent( In    ) :: k_points
    Integer                        , Intent( In    ) :: m
    Integer                        , Intent( In    ) :: n
    Type ( distributed_k_matrix   ), Intent( In    ) :: source

    Integer :: n_k_points
    Integer :: s, k, sk
    
    n_k_points = Size( k_point_type)

    ! Set up the all k point data structure
    Allocate( matrices%all_k_point_info( 1:n_spin * n_k_points ) )
    sk = 0
    Do s = 1, n_spin
       Do k = 1, n_k_points
          sk = sk + 1
          matrices%all_k_point_info( sk )%spin      = s
          matrices%all_k_point_info( sk )%k_indices = k_points( :, k )
       End Do
    End Do

    ! Now my k points
    Allocate( matrices%my_k_points( 1:n_spin * n_k_points ) )
    sk = 0
    Do s = 1, n_spin
       Do k = 1, n_k_points
          sk = sk + 1
          matrices%my_k_points( sk )%info%spin    = s
          matrices%my_k_points( sk )%info%k_indices = k_points( :, k )
          Allocate( matrices%my_k_points( sk )%data( 1:1 ) )
          matrices%my_k_points( sk )%data( 1 )%label = 1
          Call matrices%my_k_points( sk )%data( 1 )%matrix%create( k_point_type( k ) == K_POINT_COMPLEX, &
               s, k_points( :, k ), m, n, source )
          matrices%my_k_points( sk )%communicator = source%get_comm()
       End Do
    End Do

    matrices%parent_communicator = source%get_comm()
    
  End Subroutine matrices_create

  Subroutine matrices_diag( A, Q, E )

    Class( distributed_k_matrices ),                 Intent( In    ) :: A
    Type ( distributed_k_matrices ),                 Intent(   Out ) :: Q
    Type ( eval_storage           ), Dimension( : ), Intent(   Out ) :: E

    Integer :: my_ks, ks

    ! Make Q have the same set up as A
    Q = A

    Do my_ks = 1, Size( A%my_k_points )
       ! Probably want to combin ks and make evals a 1D array
       ks = A%get_sk_index( my_ks )
       If( ks /= INVALID ) Then
          Call A%my_k_points( my_ks )%data( 1 )%matrix%diag( Q%my_k_points( my_ks )%data( 1 )%matrix, E( ks )%evals )
       End If
    End Do

    ! Replicate evals
    
  End Subroutine matrices_diag

  Pure Function get_sk_index( A, my_sk ) Result( sk )

    Integer :: sk

    Class( distributed_k_matrices ), Intent( In ) :: A
    Integer                        , Intent( In ) :: my_sk

    Do sk = 1, Size( A%all_k_point_info )
       If( A%all_k_point_info( sk )%spin == A%my_k_points( my_sk )%info%spin ) Then
          If( All( A%all_k_point_info( sk )%k_indices == A%my_k_points( my_sk )%info%k_indices ) ) Then
             Exit
          End If
       End If
    End Do
    
  End Function get_sk_index

End Module k_point_matrix_module
