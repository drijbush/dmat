Module distributed_matrix_module

  Use numbers_module, Only : wp
  
  Implicit None

  Type, Abstract, Public :: distributed_matrix_ops
     Integer, Dimension( : ), Allocatable, Private :: descriptor
     Integer                             , Private :: communicator
   Contains
     Procedure( diag_interface ), Deferred :: diag
  End Type distributed_matrix_ops

  Type, Extends( distributed_matrix_ops ), Public :: real_distributed_matrix
     Real( wp ), Dimension( : ), Allocatable, Private :: data
   Contains
     Procedure  :: diag => real_distributed_matrix_diag
  End Type real_distributed_matrix

  Type, Extends( distributed_matrix_ops ), Public :: complex_distributed_matrix
     Complex( wp ), Dimension( : ), Allocatable, Private :: data
   Contains
     Procedure :: diag => complex_distributed_matrix_diag
  End Type complex_distributed_matrix

  Private  

  Abstract Interface
     Subroutine diag_interface( A, Q, E )
       Import :: distributed_matrix_ops
       Import :: wp
       Implicit None
       Class( distributed_matrix_ops ),              Intent( In    ) :: A
       Class( distributed_matrix_ops ), Allocatable, Intent(   Out ) :: Q
       Real( wp ), Dimension( : )     , Allocatable, Intent(   Out ) :: E
     End Subroutine diag_interface
  End Interface

  Integer, Parameter :: desc_n = 1

Contains

  Subroutine real_distributed_matrix_diag( A, Q, E )

    Use numbers_module, Only : wp

    Implicit None

    Class( real_distributed_matrix ),              Intent( In    ) :: A
    Class( distributed_matrix_ops  ), Allocatable, Intent(   Out ) :: Q
    Real( wp ), Dimension( : )      , Allocatable, Intent(   Out ) :: E

    Integer :: n

    Allocate( Q, Source = A )

    ! Want to change to function geting things from descriptor eventually
    n = A%descriptor( desc_n )
    
    Allocate( E( 1:n ) )

  End Subroutine real_distributed_matrix_diag
  
  Subroutine complex_distributed_matrix_diag( A, Q, E )

    Use numbers_module, Only : wp
    
    Implicit None

    Class( complex_distributed_matrix ),              Intent( In    ) :: A
    Class( distributed_matrix_ops     ), Allocatable, Intent(   Out ) :: Q
    Real( wp ), Dimension( : )         , Allocatable, Intent(   Out ) :: E

    Integer :: n

    Allocate( Q, Source = A )

    ! Want to change to function geting things from descriptor eventually
    n = A%descriptor( desc_n )
    
    Allocate( E( 1:n ) )

  End Subroutine complex_distributed_matrix_diag
  
End Module distributed_matrix_module
