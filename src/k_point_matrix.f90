Module k_point_matrix_module

  Use numbers_module           , Only : wp
  Use distributed_matrix_module, Only : distributed_matrix_ops, real_distributed_matrix, complex_distributed_matrix

  Implicit None

  Type, Public :: k_point_matrix
     Class( distributed_matrix_ops ), Allocatable, Private :: data
  End Type k_point_matrix

  Type, Extends( k_point_matrix ), Public :: k_wave_function
     Real( wp ), Dimension( : ), Allocatable, Private :: evals
  End Type k_wave_function
  
  Private

End Module k_point_matrix_module


