Module distributed_matrix_module

  Use numbers_module, Only : wp
  
  Implicit None

  Type, Abstract, Public :: distributed_matrix_ops
     Integer, Dimension( : ), Allocatable :: descriptor
     Integer                              :: communicator
   Contains
     Procedure( diag_interface ), Deferred :: diag
  End Type distributed_matrix_ops

  Type, Extends( distributed_matrix_ops ), Public :: real_distributed_matrix
     Real( wp ), Dimension( : ), Allocatable :: data
   Contains
     Procedure :: diag => real_distributed_matrix_diag
  End type real_distributed_matrix

  Private  

  Abstract Interface
     Subroutine diag_interface( A, Q, E )
       Import :: distributed_matrix_ops
       Import :: wp
       Implicit None
       Class( distributed_matrix_ops ), Intent( In    ) :: A
       Class( distributed_matrix_ops ), Intent(   Out ) :: Q
       Real( wp ), Dimension( : )     , Intent(   Out ) :: E
     End Subroutine diag_interface
  End Interface

Contains

  Subroutine real_distributed_matrix_diag( A, Q, E )
    Implicit None
    Class( real_distributed_matrix ), Intent( In    ) :: A
    Class( distributed_matrix_ops  ), Intent(   Out ) :: Q
    Real( wp ), Dimension( : )      , Intent(   Out ) :: E
  End Subroutine real_distributed_matrix_diag
  
End Module distributed_matrix_module
