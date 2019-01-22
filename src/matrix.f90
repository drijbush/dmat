Module distributed_matrix_module

  Use numbers_module, Only : wp 
  Use mapping_module, Only : mapping, mapping_get_data, mapping_get_global_n_row
  
  Implicit None

  Type, Abstract, Public :: distributed_matrix
     Type( mapping ), Private :: matrix_map
   Contains
     Procedure( diag_interface ), Deferred :: diag
  End Type distributed_matrix

  Type, Extends( distributed_matrix ), Public :: real_distributed_matrix
     Real( wp ), Dimension( : ), Allocatable, Private :: data
   Contains
     Procedure  :: diag => real_distributed_matrix_diag
  End Type real_distributed_matrix

  Type, Extends( distributed_matrix ), Public :: complex_distributed_matrix
     Complex( wp ), Dimension( : ), Allocatable, Private :: data
   Contains
     Procedure :: diag => complex_distributed_matrix_diag
  End Type complex_distributed_matrix

  Private  

  Abstract Interface
     Subroutine diag_interface( A, Q, E )
       Import :: distributed_matrix
       Import :: wp
       Implicit None
       Class( distributed_matrix ),              Intent( In    ) :: A
       Class( distributed_matrix ), Allocatable, Intent(   Out ) :: Q
       Real( wp ), Dimension( : ) , Allocatable, Intent(   Out ) :: E
     End Subroutine diag_interface
  End Interface

  Integer, Parameter :: desc_n = 1

Contains

  Subroutine real_distributed_matrix_diag( A, Q, E )

    Use numbers_module, Only : wp

    Implicit None

    Class( real_distributed_matrix ),              Intent( In    ) :: A
    Class( distributed_matrix      ), Allocatable, Intent(   Out ) :: Q
    Real( wp ), Dimension( : )      , Allocatable, Intent(   Out ) :: E

    Integer :: n

    Allocate( Q, Source = A )

    ! Want to change to function geting things from descriptor eventually
    ! When write mappings properly belongs in there
    Call mapping_get_data( A%matrix_map, mapping_get_global_n_row, n )
    
    Allocate( E( 1:n ) )

  End Subroutine real_distributed_matrix_diag
  
  Subroutine complex_distributed_matrix_diag( A, Q, E )

    Use numbers_module, Only : wp
    
    Implicit None

    Class( complex_distributed_matrix ),              Intent( In    ) :: A
    Class( distributed_matrix         ), Allocatable, Intent(   Out ) :: Q
    Real( wp ), Dimension( : )         , Allocatable, Intent(   Out ) :: E

    Integer :: n

    Allocate( Q, Source = A )

    ! Want to change to function geting things from descriptor eventually
    ! When write mappings properly belongs in there
    Call mapping_get_data( A%matrix_map, mapping_get_global_n_row, n )
    
    Allocate( E( 1:n ) )

  End Subroutine complex_distributed_matrix_diag

End Module distributed_matrix_module
