Module distributed_matrix_module

  Use numbers_module       , Only : wp 
  Use matrix_mapping_module, Only : matrix_mapping
  
  Implicit None

  Type, Abstract, Private :: matrix_data
!!$   Contains
     ! Will put in some printing stuff etc. here eventually
  End type matrix_data

  Type, Extends( matrix_data ), Private :: real_matrix_data
     Real( wp ), Dimension( :, : ), Allocatable :: data
  End type real_matrix_data

  Type, Extends( matrix_data ), Private :: complex_matrix_data
     Complex( wp ), Dimension( :, : ), Allocatable :: data
  End type complex_matrix_data

  Type, Public :: distributed_matrix
     Type ( matrix_mapping ),              Private :: matrix_map
     Class( matrix_data    ), Allocatable, Private :: data
   Contains
     Procedure :: diag => matrix_diag
  End type distributed_matrix

  Interface diag
     Procedure complex_diag
     Procedure real_diag
  End Interface diag

Contains

  Subroutine matrix_diag( A, Q, E )

    Use numbers_module, Only : wp

    Implicit None

    Class( distributed_matrix ),              Intent( In    ) :: A
    Class( distributed_matrix ), Allocatable, Intent(   Out ) :: Q
    Real( wp ), Dimension( : ) , Allocatable, Intent(   Out ) :: E

    Integer :: m, n

    Allocate( Q, Source = A )

    Call A%matrix_map%get_data( m, n )

    ! To stop whinging before implemented
    Allocate( E( 1:m ) )

    Call diag( a%data%data, q%data%data, e )

  End Subroutine matrix_diag

  Subroutine complex_diag( a, q, e )
    Complex( wp ), Dimension( :, : ), Intent( In    ) :: a
    Complex( wp ), Dimension( :, : ), Intent(   Out ) :: q
    Real   ( wp ), Dimension( :    ), Intent(   Out ) :: e
  End Subroutine complex_diag
     
  Subroutine real_diag( a, q, e )
    Real( wp ), Dimension( :, : ), Intent( In    ) :: a
    Real( wp ), Dimension( :, : ), Intent(   Out ) :: q
    Real( wp ), Dimension( :    ), Intent(   Out ) :: e
  End Subroutine real_diag
     
End Module distributed_matrix_module
