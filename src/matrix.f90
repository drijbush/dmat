Module distributed_matrix_module

  Use numbers_module       , Only : wp
  Use Scalapack_interfaces 
  Use matrix_mapping_module, Only : matrix_mapping, matrix_mapping_base_start, matrix_mapping_base, &
       matrix_mapping_init, matrix_mapping_finalise
  
  Implicit None

  Type, Abstract, Public :: distributed_matrix
     Type( matrix_mapping ) :: matrix_map
     ! Will put in some printing stuff etc. here eventually
  End type distributed_matrix

  Type, Extends( distributed_matrix ), Public :: base_distributed_matrix
   Contains
     Procedure :: create => matrix_create
  End type base_distributed_matrix

  Type( base_distributed_matrix ), Parameter :: base_matrix_start = &
       base_distributed_matrix( matrix_map = matrix_mapping_base_start )
  Type( base_distributed_matrix ), Public :: base_matrix = base_matrix_start

  Type, Extends( base_distributed_matrix ), Public :: real_distributed_matrix
     Real( wp ), Dimension( :, : ), Allocatable :: data
   Contains
     Procedure :: diag   => matrix_diag_real
  End type real_distributed_matrix

  Type, Extends( base_distributed_matrix ), Public :: complex_distributed_matrix
     Complex( wp ), Dimension( :, : ), Allocatable :: data
   Contains
     Procedure :: diag => matrix_diag_complex
  End type complex_distributed_matrix

  Public :: distributed_matrix_init
  Public :: distributed_matrix_finalise
  
  Private

  Integer, Parameter :: diag_work_size_fiddle_factor = 4 ! From experience Scalapack sometimes returns too small a work size
  
  Integer, Parameter, Private :: default_block_fac = 96
  Integer,            Private :: block_fac = default_block_fac

Contains

  Subroutine distributed_matrix_init( comm )

    Integer, Intent( In ) :: comm

    Call matrix_mapping_init( comm )

    base_matrix = base_matrix_start

  End Subroutine distributed_matrix_init

  Subroutine distributed_matrix_finalise

    Call matrix_mapping_finalise
    
    base_matrix = base_matrix_start

  End Subroutine distributed_matrix_finalise

  Subroutine matrix_create( matrix, m, n, source_matrix )

    Class( base_distributed_matrix ), Intent(   Out ) :: matrix
    Integer                         , Intent( In    ) :: m
    Integer                         , Intent( In    ) :: n
    Class( base_distributed_matrix ), Intent( In    ) :: source_matrix

    Integer :: nprow, myprow, mb, lda
    Integer :: npcol, mypcol, nb, sda

    matrix%matrix_map = source_matrix%matrix_map

    ! Need to fix if n, m smaller than blocking fac
    mb = block_fac
    nb = block_fac
    mb = Min( mb, nb )
    nb = mb

    Call matrix%matrix_map%get_data( nprow = nprow, myprow = myprow )
    lda = numroc( m, mb, myprow, 0, nprow )

    Call matrix%matrix_map%get_data( npcol = npcol, mypcol = mypcol )
    sda = numroc( n, nb, mypcol, 0, npcol )

    Call matrix%matrix_map%set( matrix%matrix_map%proc_mapping, m, n, mb, nb, 0, 0, lda )

    Select Type( matrix )
    Class Default
       Stop "Illegal type in matrix_create"
    Class is ( real_distributed_matrix )
       Allocate( matrix%data( 1:lda, 1:sda  ) )
    Class is ( complex_distributed_matrix )
       Allocate( matrix%data( 1:lda, 1:sda  ) )
    End Select

  End Subroutine matrix_create
  
  Subroutine matrix_diag_real( A, Q, E )

    Use numbers_module, Only : wp

    Implicit None

    Class( real_distributed_matrix ),              Intent( In    ) :: A
    Class( real_distributed_matrix ), Allocatable, Intent(   Out ) :: Q
    Real( wp ), Dimension( : )      , Allocatable, Intent(   Out ) :: E

    Real( wp ), Dimension( :, : ), Allocatable :: tmp_a

    Real( wp ), Dimension( : ), Allocatable :: work

    Integer, Dimension( : ), Allocatable :: iwork
    
    Integer :: nwork
    Integer :: npcol
    Integer :: m, n
    Integer :: info

    Allocate( Q, Source = A )

    Call A%matrix_map%get_data( m = m, n = n, npcol = npcol )

    Allocate( E( 1:m ) )

    ! The diag overwrites the matrix. Horrible so use a temporary
    tmp_A = A%data

    ! Workspace size enquiry
    Call pdsyevd( 'V', 'U', m, tmp_A, 1, 1, A%matrix_map%get_descriptor(), E, Q%data, 1, 1, Q%matrix_map%get_descriptor(), &
         work, -1, iwork, 0, info )
    nwork = Nint( work( 1 ) )
    nwork = nwork * diag_work_size_fiddle_factor ! From experience ...
    Allocate(  work( 1:nwork ) )
    Allocate( iwork( 1:7 * m + 8 * npcol + 2 ) )
    ! Do the diag
    Call pdsyevd( 'V', 'U', m, tmp_A, 1, 1, A%matrix_map%get_descriptor(), E, Q%data, 1, 1, Q%matrix_map%get_descriptor(), &
         work, Size( work ), iwork, Size( iwork ), info )

  End Subroutine matrix_diag_real

  Subroutine matrix_diag_complex( A, Q, E )

    Use numbers_module, Only : wp

    Implicit None

    Class( complex_distributed_matrix ),              Intent( In    ) :: A
    Class( complex_distributed_matrix ), Allocatable, Intent(   Out ) :: Q
    Real( wp ), Dimension( : )         , Allocatable, Intent(   Out ) :: E

    Complex( wp ), Dimension( :, : ), Allocatable :: tmp_a

    Complex( wp ), Dimension( : ), Allocatable :: cwork

    Real( wp ), Dimension( : ), Allocatable :: rwork

    Integer, Dimension( : ), Allocatable :: iwork
    
    Integer :: ncwork, nrwork
    Integer :: npcol
    Integer :: m, n
    Integer :: info

    Allocate( Q, Source = A )

    Call A%matrix_map%get_data( m = m, n = n, npcol = npcol )

    Allocate( E( 1:m ) )

    ! The diag overwrites the matrix. Horrible so use a temporary
    tmp_A = A%data

    ! Workspace size enquiry
    Call pzheevd( 'V', 'U', m, tmp_A, 1, 1, A%matrix_map%get_descriptor(), E, Q%data, 1, 1, Q%matrix_map%get_descriptor(), &
         cwork, -1, rwork, -1, iwork, 0, info )
    ncwork = Nint( Real( cwork( 1 ), wp ) )
    ncwork = ncwork * diag_work_size_fiddle_factor ! From experience ...
    Allocate( cwork( 1:ncwork ) )
    nrwork = Nint( rwork( 1 ) )
    nrwork = nrwork * diag_work_size_fiddle_factor ! From experience ...
    Allocate( rwork( 1:nrwork ) )
    Allocate( iwork( 1:7 * m + 8 * npcol + 2 ) )
    ! Do the diag
    Call pzheevd( 'V', 'U', m, tmp_A, 1, 1, A%matrix_map%get_descriptor(), E, Q%data, 1, 1, Q%matrix_map%get_descriptor(), &
         cwork, Size( cwork ), rwork, Size( rwork ), iwork, Size( iwork ), info )

  End Subroutine matrix_diag_complex

     
End Module distributed_matrix_module
 
