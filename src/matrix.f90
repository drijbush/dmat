Module distributed_matrix_module

  Use numbers_module       , Only : wp
  Use Scalapack_interfaces 
  Use matrix_mapping_module, Only : matrix_mapping, matrix_mapping_init, matrix_mapping_finalise
  
  Implicit None

  Integer, Parameter :: distributed_matrix_INVALID = -1
  Integer, Parameter :: distributed_matrix_NOT_ME  = -2
  
  Type, Abstract, Public :: distributed_matrix
     Type( matrix_mapping ) :: matrix_map
     Integer, Dimension( : ), Allocatable :: global_to_local_rows
     Integer, Dimension( : ), Allocatable :: global_to_local_cols
     Integer, Dimension( : ), Allocatable :: local_to_global_rows
     Integer, Dimension( : ), Allocatable :: local_to_global_cols
   Contains
     Procedure :: create => matrix_create
  End type distributed_matrix

  Type, Extends( distributed_matrix ), Public :: real_distributed_matrix
     Real( wp ), Dimension( :, : ), Allocatable :: data
   Contains
     Procedure :: diag          => matrix_diag_real
     Procedure :: set_by_global => matrix_set_global_real
     Procedure :: set_by_local  => matrix_set_local_real
  End type real_distributed_matrix

  Type, Extends( distributed_matrix ), Public :: complex_distributed_matrix
     Complex( wp ), Dimension( :, : ), Allocatable :: data
   Contains
     Procedure :: diag          => matrix_diag_complex
     Procedure :: set_by_global => matrix_set_global_complex
     Procedure :: set_by_local  => matrix_set_local_complex
  End type complex_distributed_matrix

  Public :: distributed_matrix_init
  Public :: distributed_matrix_finalise
  Public :: distributed_matrix_set_default_blocking
  
  Private

  Integer, Parameter, Private :: diag_work_size_fiddle_factor = 4 ! From experience Scalapack sometimes returns too small a work size
  
  Integer, Parameter, Private :: default_block_fac = 4
  Integer,            Private :: block_fac = default_block_fac

Contains

  Subroutine distributed_matrix_init( comm, base_matrix )

    Class  ( distributed_matrix ), Intent(   Out ) :: base_matrix 
    Integer                      , Intent( In    ) :: comm

    Type( matrix_mapping ) :: base_matrix_mapping
    
    Call matrix_mapping_init( comm, base_matrix_mapping )

    base_matrix%matrix_map = base_matrix_mapping

    base_matrix%global_to_local_rows = [ distributed_matrix_INVALID ]
    base_matrix%global_to_local_cols = [ distributed_matrix_INVALID ]
    base_matrix%local_to_global_rows = [ distributed_matrix_INVALID ]
    base_matrix%local_to_global_cols = [ distributed_matrix_INVALID ]
    
  End Subroutine distributed_matrix_init

  Subroutine distributed_matrix_finalise

    Call matrix_mapping_finalise
    
  End Subroutine distributed_matrix_finalise

  Subroutine distributed_matrix_set_default_blocking( bfac )

    Integer, Intent( In ) :: bfac

    block_fac = bfac
    
  End Subroutine distributed_matrix_set_default_blocking

  Subroutine matrix_set_global_real( matrix, m, n, p, q, data )

    ! Sets the data ( m:n, p:q ) in the global matrix

    Class( real_distributed_matrix ), Intent( InOut ) :: matrix
    Integer                         , Intent( In    ) :: m
    Integer                         , Intent( In    ) :: n
    Integer                         , Intent( In    ) :: p
    Integer                         , Intent( In    ) :: q
    Real( wp ), Dimension( m:, p: ) , Intent( In    ) :: data

    Integer :: i_glob, j_glob
    Integer :: i_loc , j_loc
    
    ! THIS NEEDS OPTIMISATION!!

    Do j_glob = p, q
       j_loc = matrix%global_to_local_cols( j_glob )
       If( j_loc == distributed_matrix_NOT_ME ) Cycle
       Do i_glob = m, n
          i_loc = matrix%global_to_local_rows( i_glob )
          If( i_loc == distributed_matrix_NOT_ME ) Cycle
          matrix%data( i_loc, j_loc ) = data( j_glob, i_glob )
       End Do
    End Do
       
  End Subroutine matrix_set_global_real

  Subroutine matrix_set_local_real( matrix, m, n, p, q, data )

    ! Sets the data ( m:n, p:q ) in the local matrix

    Class( real_distributed_matrix ), Intent( InOut ) :: matrix
    Integer                         , Intent( In    ) :: m
    Integer                         , Intent( In    ) :: n
    Integer                         , Intent( In    ) :: p
    Integer                         , Intent( In    ) :: q
    Real( wp ), Dimension( m:, p: ) , Intent( In    ) :: data

    matrix%data( m:n, p:q ) = data( m:n, p:q )
    
  End Subroutine matrix_set_local_real

  Subroutine matrix_set_global_complex( matrix, m, n, p, q, data )

    ! Sets the data ( m:n, p:q ) in the global matrix

    Class( complex_distributed_matrix ), Intent( InOut ) :: matrix
    Integer                            , Intent( In    ) :: m
    Integer                            , Intent( In    ) :: n
    Integer                            , Intent( In    ) :: p
    Integer                            , Intent( In    ) :: q
    Complex( wp ), Dimension( m:, p: ) , Intent( In    ) :: data

    Integer :: i_glob, j_glob
    Integer :: i_loc , j_loc
    
    ! THIS NEEDS OPTIMISATION!!

    Do j_glob = p, q
       j_loc = matrix%global_to_local_cols( j_glob )
       If( j_loc == distributed_matrix_NOT_ME ) Cycle
       Do i_glob = m, n
          i_loc = matrix%global_to_local_rows( i_glob )
          If( i_loc == distributed_matrix_NOT_ME ) Cycle
          matrix%data( i_loc, j_loc ) = data( j_glob, i_glob )
       End Do
    End Do
       
  End Subroutine matrix_set_global_complex
  
  Subroutine matrix_set_local_complex( matrix, m, n, p, q, data )

    ! Sets the data ( m:n, p:q ) in the local matrix

    Class( complex_distributed_matrix ), Intent( InOut ) :: matrix
    Integer                            , Intent( In    ) :: m
    Integer                            , Intent( In    ) :: n
    Integer                            , Intent( In    ) :: p
    Integer                            , Intent( In    ) :: q
    Complex( wp ), Dimension( m:, p: ) , Intent( In    ) :: data

    matrix%data( m:n, p:q ) = data( m:n, p:q )
    
  End Subroutine matrix_set_local_complex

  Subroutine matrix_create( matrix, m, n, source_matrix )

    Class( distributed_matrix ), Intent(   Out ) :: matrix
    Integer                    , Intent( In    ) :: m
    Integer                    , Intent( In    ) :: n
    Class( distributed_matrix ), Intent( In    ) :: source_matrix

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

    Call set_local_to_global( matrix%local_to_global_rows, m, mb, myprow, nprow, lda )
    Call set_local_to_global( matrix%local_to_global_cols, n, nb, mypcol, npcol, sda )
    Write( *, * ) 'l->g', m, mb, myprow, nprow, lda, matrix%local_to_global_rows
    Write( *, * ) 'l->g', n, nb, mypcol, npcol, sda, matrix%local_to_global_cols

    Call set_global_to_local( matrix%global_to_local_rows, m, mb, myprow, nprow )
    Call set_global_to_local( matrix%global_to_local_cols, n, nb, mypcol, npcol )
    Write( *, * ) 'g->l', m, mb, myprow, nprow, lda, matrix%global_to_local_rows
    Write( *, * ) 'g->l', n, nb, mypcol, npcol, sda, matrix%global_to_local_cols

  End Subroutine matrix_create

  Subroutine set_local_to_global( loc_to_glob, n, nb, myp, np, da )

    Integer, Dimension( : ), Allocatable, Intent(   Out ) :: loc_to_glob
    Integer                             , Intent( In    ) :: n
    Integer                             , Intent( In    ) :: nb
    Integer                             , Intent( In    ) :: myp
    Integer                             , Intent( In    ) :: np
    Integer                             , Intent( In    ) :: da

    Integer :: i_glob, i_loc, skip, start

    Allocate( loc_to_glob( 1:da ) )

    skip =  np * nb

    i_loc = 1
    start = myp * nb + 1
    Do While( start <= n )
       Do i_glob = start, Min( start + nb - 1, n )
          loc_to_glob( i_loc ) = i_glob
          i_loc = i_loc + 1
       End Do
       start = start + skip
    End Do

  End Subroutine set_local_to_global
  
  Subroutine set_global_to_local( glob_to_loc, n, nb, myp, np )

    Integer, Dimension( : ), Allocatable, Intent(   Out ) :: glob_to_loc
    Integer                             , Intent( In    ) :: n
    Integer                             , Intent( In    ) :: nb
    Integer                             , Intent( In    ) :: myp
    Integer                             , Intent( In    ) :: np

    Integer :: i_glob, i_loc, skip, start

    Allocate( glob_to_loc( 1:n ) )

    glob_to_loc = distributed_matrix_NOT_ME
    
    skip =  np * nb

    i_loc = 1
    start = myp * nb + 1
    Do While( start <= n )
       Do i_glob = start, Min( start + nb - 1, n )
          glob_to_loc( i_glob ) = i_loc
          i_loc = i_loc + 1
       End Do
       start = start + skip
    End Do

  End Subroutine set_global_to_local
  
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
 
