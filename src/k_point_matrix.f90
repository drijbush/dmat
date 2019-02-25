Module k_point_matrix_module

  Use numbers_module      , Only : wp
  Use distributed_k_module, Only : distributed_k_matrix_init, &
       distributed_k_matrix_finalise, distributed_k_matrix

  Implicit None

  Integer, Public, Parameter :: K_POINT_REAL = 0
  Integer, Public, Parameter :: K_POINT_COMPLEX = 1

  Integer, Private, Parameter :: INVALID = -1

  Type k_point_data
     Integer                              :: spin
     Integer, Dimension( : ), Allocatable :: indices
  End type k_point_data

  Type, Private :: k_point_matrices
     ! External label - e.g. irrep number
     Integer                      :: label = INVALID
     Type( distributed_k_matrix ) :: matrix
  End type k_point_matrices

  Type, Private :: k_point
     Integer                 , Dimension( : ), Allocatable :: spin
     Integer                 , Dimension( : ), Allocatable :: k
     ! This splitting allows irreps
     Type( k_point_matrices ), Dimension( : ), Allocatable :: data
     ! Want to hide eventually
     Integer                                               :: communicator = INVALID
  End type k_point

  Type, Public :: distributed_k_matrices
     Type( k_point_data ), Dimension( : ), Allocatable :: all_k_point_data
     ! this splitting allows multiple k points on this process
     Type( k_point      ), Dimension( : ), Allocatable :: my_k_points
     ! Want to hide eventually
     Integer                                           :: parent_communicator = INVALID
   Contains
     Procedure :: create => matrices_create
  End type distributed_k_matrices
  
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
    Allocate( matrices%all_k_point_data( 1:n_spin * n_k_points ) )
    sk = 0
    Do s = 1, n_spin
       Do k = 1, n_k_points
          sk = sk + 1
          matrices%all_k_point_data( sk )%spin    = s
          matrices%all_k_point_data( sk )%indices = k_points( :, k )
       End Do
    End Do

    ! Now my k points
    Allocate( matrices%my_k_points( 1:n_spin * n_k_points ) )
    sk = 0
    Do s = 1, n_spin
       Do k = 1, n_k_points
          sk = sk + 1
          matrices%my_k_points( sk )%spin = s
          matrices%my_k_points( sk )%k    = k
          Allocate( matrices%my_k_points( sk )%data( 1:1 ) )
          matrices%my_k_points( sk )%data( 1 )%label = 1
          Call matrices%my_k_points( sk )%data( 1 )%matrix%create( k_point_type( k ) == K_POINT_COMPLEX, &
               s, k_points( :, k ), m, n, source )
          matrices%my_k_points( sk )%communicator = source%get_comm()
       End Do
    End Do

    matrices%parent_communicator = source%get_comm()
    
  End Subroutine matrices_create

End Module k_point_matrix_module
