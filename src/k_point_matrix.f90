Module k_point_matrix_module

  Use numbers_module           , Only : wp
  Use distributed_matrix_module, Only : distributed_matrix, real_distributed_matrix, complex_distributed_matrix

  Implicit None

  Type, Public :: k_point_matrix
     Integer                   :: spin
     Integer, Dimension( 1:3 ) :: k_point
     Class( distributed_matrix ), Allocatable, Private :: matrix
  End Type k_point_matrix

  Type, Extends( k_point_matrix ), Public :: k_wave_function
     Real( wp ), Dimension( : ), Allocatable, Private :: evals
  End Type k_wave_function
  
  Private

Contains

  Subroutine testit
    Type( k_wave_function ) :: a
    Class( distributed_matrix ), Allocatable :: b
    Allocate( complex_distributed_matrix :: a%matrix )
    Call a%matrix%diag( b, a%evals )
  End Subroutine testit

End Module k_point_matrix_module


