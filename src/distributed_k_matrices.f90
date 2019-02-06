Module distributed_k_module

  Use numbers_module       , Only : wp
  Use distributed_matrix_module, Only : distributed_matrix, real_distributed_matrix, complex_distributed_matrix

  Implicit None

  Type, Private :: k_point_matrix
     Integer                   :: spin
     Integer, Dimension( 1:3 ) :: k_point
     Class( distributed_matrix ), Allocatable, Private :: matrix
  End Type k_point_matrix

  Type, Extends( k_point_matrix ), Private :: k_wave_function
     Real( wp ), Dimension( : ), Allocatable, Private :: evals
  End Type k_wave_function

  Type, Public :: distributed_k_matrix
     Class( k_point_matrix ), Dimension( : ), Allocatable :: k_points
  End Type distributed_k_matrix
  
  ! a) Set up base matrix which is 1 matrix across all
  ! b) Provide splits
  
End Module distributed_k_module
