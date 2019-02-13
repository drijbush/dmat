Module distributed_k_module

  ! Wrapper for dist matrix types so can make opaque array - only 1 k point at this level

  Use numbers_module       , Only : wp
  Use distributed_matrix_module, Only : distributed_matrix, real_distributed_matrix, complex_distributed_matrix, &
       distributed_matrix_init, distributed_matrix_finalise

  Implicit None

  Type, Private :: k_point_matrix
     Integer                   :: this_spin
     Integer, Dimension( 1:3 ) :: this_k_point
     Class( distributed_matrix ), Allocatable, Private :: matrix
  End Type k_point_matrix

  Type, Extends( k_point_matrix ), Private :: k_wave_function
     Real( wp ), Dimension( : ), Allocatable, Private :: evals
  End Type k_wave_function

  Type, Public :: distributed_k_matrix
     Class( k_point_matrix ), Allocatable :: k_point
   Contains
     Procedure :: create        => distributed_k_matrix_create
!!$     Procedure :: diag          => distributed_k_matrix_diag
  End Type distributed_k_matrix
  
  ! a) Set up base matrix which is 1 matrix across all
  ! b) Provide splits

  Integer, Parameter :: INVALID = -1

Contains

  Subroutine distributed_k_matrix_init( comm, base_matrix )

    Integer                        , Intent( In    ) :: comm
    Type   ( distributed_k_matrix ), Intent(   Out ) :: base_matrix

    base_matrix%k_point%this_spin    = INVALID
    base_matrix%k_point%this_k_point = INVALID

    Allocate( real_distributed_matrix:: base_matrix%k_point%matrix )

    Call distributed_matrix_init( comm, base_matrix%k_point%matrix )
    
  End Subroutine distributed_k_matrix_init

  Subroutine distributed_k_matrix_create( matrix, is_complex, this_spin, this_k_point, m, n, source_matrix )

    Class  ( distributed_k_matrix ), Intent(   Out ) :: matrix
    Logical                        , Intent( In    ) :: is_complex
    Integer                        , Intent( In    ) :: this_spin
    Integer, Dimension( 1:3 )      , Intent( In    ) :: this_k_point
    Integer                        , Intent( In    ) :: m
    Integer                        , Intent( In    ) :: n
    Type   ( distributed_k_matrix ), Intent(   Out ) :: source_matrix

    If( is_complex ) Then
      Allocate( complex_distributed_matrix :: matrix%k_point%matrix ) 
    Else
      Allocate( real_distributed_matrix    :: matrix%k_point%matrix ) 
    End If

    matrix%k_point%this_spin    = this_spin
    matrix%k_point%this_k_point = this_k_point
    
    Call matrix%k_point%matrix%create( m, n, source_matrix%k_point%matrix )
    
  End Subroutine distributed_k_matrix_create

  Subroutine distributed_k_matrix_diag( A, Q )
    
    Class( distributed_k_matrix ), Intent( In    ) :: A
    Type ( distributed_k_matrix ), Intent(   Out ) :: Q

    Type( real_distributed_matrix ), Allocatable :: Br
    Type( complex_distributed_matrix ), Allocatable :: Bc

    Real( wp ), Dimension( : ), Allocatable  :: evals
    
    Allocate( k_wave_function  :: Q%k_point )

    Associate( kA => A%k_point, kQ => Q%k_point, &
         kAm => A%k_point%matrix, kQm => Q%k_point%matrix )

      Select Type( kQ ) !! hack!!!!

      Class Default
         Stop "Illegal type  in distributed_k_matrix_diag"

      Class is ( k_wave_function  )

         Select Type( kAm )

         Class Default
            Stop "Illegal type  in distributed_k_matrix_diag"

         Class is ( real_distributed_matrix )

            Select Type( kQm )
            Class Default
               Stop "Illegal type  in distributed_k_matrix_diag"
            Class is ( real_distributed_matrix )
               Call kAm%diag( Br, evals )
               Allocate( Q%k_point%matrix, source = Br )
               kQ%evals  = evals
            End Select

         Class is ( complex_distributed_matrix )
            Select Type( kQm )
            Class Default
               Stop "Illegal type  in distributed_k_matrix_diag"
            Class is ( real_distributed_matrix )
               Call kAm%diag( Bc, evals )
               Allocate( Q%k_point%matrix, source = Bc )
               kQ%evals  = evals
            End Select
         End Select
      End Select
      
    End Associate

  End Subroutine distributed_k_matrix_diag
  
  Subroutine distributed_k_matrix_finalize

    Call distributed_matrix_finalise
    
  End Subroutine distributed_k_matrix_finalize
  
End Module distributed_k_module
