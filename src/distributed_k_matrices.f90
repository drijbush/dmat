Module distributed_k_module

  Use numbers_module       , Only : wp
  Use k_point_matrix_module, Only : k_point_matrix, k_wave_function

  Implicit None

  Type, Public :: distributed_k_matrix
     Type( k_point_matrix ), Dimension( : ), Allocatable :: k_points
  End Type distributed_k_matrix
  
  Type, Public :: distributed_k_vectors
     Type( k_point_matrix ), Dimension( : ), Allocatable :: k_points
  End Type distributed_k_vectors


  ! a) Set up base matrix which is 1 matrix across all
  ! b) Provide splits
  
End Module distributed_k_module
