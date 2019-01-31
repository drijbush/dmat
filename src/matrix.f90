Module distributed_matrix_module

  Use numbers_module       , Only : wp 
  Use matrix_mapping_module, Only : matrix_mapping
  
  Implicit None

  Type, Abstract, Public :: distributed_matrix
     Type( matrix_mapping ) :: matrix_map
     ! Will put in some printing stuff etc. here eventually
  End type distributed_matrix

  Type, Extends( distributed_matrix ), Public :: real_distributed_matrix
     Real( wp ), Dimension( :, : ), Allocatable :: data
   Contains
     Procedure :: diag => matrix_diag_real
  End type real_distributed_matrix

  Type, Extends( distributed_matrix ), Public :: complex_distributed_matrix
     Complex( wp ), Dimension( :, : ), Allocatable :: data
   Contains
     Procedure :: diag => matrix_diag_complex
  End type complex_distributed_matrix

  Private

  Integer, Parameter :: diag_work_size_fiddle_factor = 4 ! From experience Scalapack sometimes returns too small a work size
  
Contains

  Subroutine matrix_diag_real( A, Q, E )

    Use numbers_module, Only : wp

    Implicit None

    Class( real_distributed_matrix ),              Intent( In    ) :: A
    Class( real_distributed_matrix ), Allocatable, Intent(   Out ) :: Q
    Real( wp ), Dimension( : )      , Allocatable, Intent(   Out ) :: E

    Real( wp ), Dimension( : ), Allocatable :: work
    
    Real( wp ) :: work_size

    Integer, Dimension( : ), Allocatable :: iwork
    
    Integer :: nwork
    Integer :: npcol
    Integer :: m, n
    Integer :: info

    Allocate( Q, Source = A )

    Call A%matrix_map%get_data( m = m, n = n, npcol = npcol )

    ! Do we need to allocate Q%data or does the source take care of that???????

    Allocate( E( 1:m ) )

    ! Workspace size enquiry
    Call pdsyevd( 'V', 'U', m, A%data, 1, 1, A%matrix_map%get_descriptor(), E, Q%data, 1, 1, Q%matrix_map%get_descriptor(), &
         work_size, -1, iwork, 0, info )
    nwork = Nint( work( 1 ) )
    nwork = nwork * diag_work_size_fiddle_factor ! From experience ...
    Allocate(  work( 1:nwork ) )
    Allocate( iwork( 1:7 * m + 8 * npcol + 2 ) )
    ! Do the diag
    Call pdsyevd( 'V', 'U', m, A%data, 1, 1, A%matrix_map%get_descriptor(), E, Q%data, 1, 1, Q%matrix_map%get_descriptor(), &
         work, -1, iwork, Size( iwork ), info )

  End Subroutine matrix_diag_real

  Subroutine matrix_diag_complex( A, Q, E )

    Use numbers_module, Only : wp

    Implicit None

    Class( complex_distributed_matrix ),              Intent( In    ) :: A
    Class( complex_distributed_matrix ), Allocatable, Intent(   Out ) :: Q
    Real( wp ), Dimension( : )         , Allocatable, Intent(   Out ) :: E

    Integer :: m, n

    Allocate( Q, Source = A )

    Call A%matrix_map%get_data( m, n )

    Allocate( E( 1:m ) )

  End Subroutine matrix_diag_complex

     
End Module distributed_matrix_module
