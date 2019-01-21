Module distributed_matrix_module

  Use numbers_module, Only : wp
  
  Implicit None

  Type, Abstract :: distributed_matrix_ops
   Contains
     Procedure( diag_interface ), Deferred :: diag
  End Type distributed_matrix_ops

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
  
End Module distributed_matrix_module
