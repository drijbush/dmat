Module distributed_matrix_module

  ! General thoughts: Container for data and overload for unified data type???

  Use mpi
  
  Use numbers_module       , Only : wp
  Use Scalapack_interfaces 
  Use matrix_mapping_module, Only : matrix_mapping, matrix_mapping_init, matrix_mapping_finalise

  
  Implicit None

  Integer, Parameter :: distributed_matrix_INVALID = -1
  Integer, Parameter :: distributed_matrix_NOT_ME  = -2
  
  Type, Abstract, Public :: distributed_matrix
     Type( matrix_mapping )               :: matrix_map
     Integer, Dimension( : ), Allocatable :: global_to_local_rows
     Integer, Dimension( : ), Allocatable :: global_to_local_cols
     Integer, Dimension( : ), Allocatable :: local_to_global_rows
     Integer, Dimension( : ), Allocatable :: local_to_global_cols
     Logical                              :: daggered = .False. 
   Contains
     Procedure          :: create        => matrix_create
     Procedure          :: get_maps      => matrix_get_maps
!!$     Procedure, Private :: dummy
!!$     Procedure, Private :: dummy_f
!!$     Generic            :: diag          => dummy
!!$     Generic            :: dagger        => dummy_f
!!$     Generic            :: multiply      => dummy_f
!!$     Generic            :: set_by_global => dummy
!!$     Generic            :: set_by_local  => dummy
!!$     Generic            :: get_by_global => dummy
!!$     Generic            :: get_by_local  => dummy
!!$     Generic            :: extract_cols  => dummy
  End type distributed_matrix

  Type, Extends( distributed_matrix ), Public :: real_distributed_matrix
     Real( wp ), Dimension( :, : ), Allocatable :: data
   Contains
     Procedure, Private   :: diag_r               => matrix_diag_real
     Generic              :: diag                 => diag_r
     Procedure, Private   :: dagger_r             => matrix_dagger_real
     Generic              :: dagger               => dagger_r
     Generic              :: Operator( .Dagger. ) => dagger_r
     Procedure, Private   :: multiply_r           => matrix_multiply_real
     Generic              :: multiply             => multiply_r
     Generic              :: Operator( * )        => multiply_r
     Procedure, Pass( A ) :: pre_scale            => matrix_pre_scale_real
     Procedure            :: post_scale           => matrix_post_scale_real
     Generic              :: Operator( * )        => pre_scale, post_scale
     Procedure            :: add                  => matrix_add_real
     Generic              :: Operator( + )        => add
     Procedure, Private   :: set_by_global_r      => matrix_set_global_real
     Procedure, Private   :: set_by_local_r       => matrix_set_local_real
     Procedure, Private   :: get_by_global_r      => matrix_get_global_real
     Procedure, Private   :: get_by_local_r       => matrix_get_local_real
     Generic              :: set_by_global        => set_by_global_r
     Generic              :: set_by_local         => set_by_local_r
     Generic              :: get_by_global        => get_by_global_r
     Generic              :: get_by_local         => get_by_local_r
     Procedure, Private   :: extract_cols_r       => matrix_extract_cols_real
     Generic              :: extract_cols         => extract_cols_r
  End type real_distributed_matrix

  Type, Extends( distributed_matrix ), Public :: complex_distributed_matrix
     Complex( wp ), Dimension( :, : ), Allocatable :: data
   Contains
     Procedure, Private   :: diag_c               => matrix_diag_complex
     Generic              :: diag                 => diag_c
     Procedure, Private   :: dagger_c             => matrix_dagger_complex
     Generic              :: dagger               => dagger_c 
     Generic              :: Operator( .Dagger. ) => dagger_c
     Procedure, Private   :: multiply_c           => matrix_multiply_complex
     Generic              :: multiply             => multiply_c
     Generic              :: Operator( * )        => multiply_c
     Procedure, Pass( A ) :: pre_scale            => matrix_pre_scale_complex
     Procedure            :: post_scale           => matrix_post_scale_complex
     Generic              :: Operator( * )        => pre_scale, post_scale
     Procedure            :: add                  => matrix_add_complex
     Generic              :: Operator( + )        => add
     Procedure, Private   :: set_by_global_c      => matrix_set_global_complex
     Procedure, Private   :: set_by_local_c       => matrix_set_local_complex
     Procedure, Private   :: get_by_global_c      => matrix_get_global_complex
     Procedure, Private   :: get_by_local_c       => matrix_get_local_complex
     Generic              :: set_by_global        => set_by_global_c
     Generic              :: set_by_local         => set_by_local_c
     Generic              :: get_by_global        => get_by_global_c
     Generic              :: get_by_local         => get_by_local_c
     Procedure, Private   :: extract_cols_c       => matrix_extract_cols_complex
     Generic              :: extract_cols         => extract_cols_c
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

    Integer                      , Intent( In    ) :: comm
    Class  ( distributed_matrix ), Intent(   Out ) :: base_matrix 

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

  Subroutine matrix_create( matrix, m, n, source_matrix )

    Class( distributed_matrix ), Intent(   Out ) :: matrix
    Integer                    , Intent( In    ) :: m
    Integer                    , Intent( In    ) :: n
    Class( distributed_matrix ), Intent( In    ) :: source_matrix

    Integer :: nprow, myprow, mb, lda
    Integer :: npcol, mypcol, nb, sda
    Integer :: ctxt

    ! Need to fix if n, m smaller than blocking fac
    mb = block_fac
    nb = block_fac
    mb = Min( mb, nb )
    nb = mb

    Call source_matrix%matrix_map%get_data( nprow = nprow, myprow = myprow )
    lda = numroc( m, mb, myprow, 0, nprow )

    Call source_matrix%matrix_map%get_data( npcol = npcol, mypcol = mypcol )
    sda = numroc( n, nb, mypcol, 0, npcol )

    Call source_matrix%matrix_map%get_data( ctxt = ctxt )

    Call matrix%matrix_map%set( matrix%matrix_map%proc_mapping, ctxt, m, n, mb, nb, 0, 0, lda )

    Call set_local_to_global( matrix%local_to_global_rows, m, mb, myprow, nprow, lda )
    Call set_local_to_global( matrix%local_to_global_cols, n, nb, mypcol, npcol, sda )

    Call set_global_to_local( matrix%global_to_local_rows, m, mb, myprow, nprow )
    Call set_global_to_local( matrix%global_to_local_cols, n, nb, mypcol, npcol )

    Select Type( matrix )
    Class Default
       Stop "Illegal type in matrix_create"
    Class is ( real_distributed_matrix )
       Allocate( matrix%data( 1:lda, 1:sda  ) )
    Class is ( complex_distributed_matrix )
       Allocate( matrix%data( 1:lda, 1:sda  ) )
    End Select

  Contains

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
      
      glob_to_loc = DISTRIBUTED_MATRIX_NOT_ME
      
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
    
  End Subroutine matrix_create

  Pure Function matrix_dagger_real( matrix ) Result( tm )

    Class( real_distributed_matrix ), Allocatable :: tm

    Class( real_distributed_matrix ), Intent( In ) :: matrix

    Allocate( tm, Source = matrix )
    tm%daggered = .Not. tm%daggered
    
  End Function matrix_dagger_real

  Pure Function matrix_dagger_complex( matrix ) Result( tm )

    Class( complex_distributed_matrix ), Allocatable :: tm

    Class( complex_distributed_matrix ), Intent( In ) :: matrix

    Allocate( tm, Source = matrix )
    tm%daggered = .Not. tm%daggered
    
  End Function matrix_dagger_complex

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
          matrix%data( i_loc, j_loc ) = data( i_glob, j_glob )
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
          matrix%data( i_loc, j_loc ) = data( i_glob, j_glob )
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

  Subroutine matrix_get_global_real( matrix, m, n, p, q, data )

    ! Gets the data ( m:n, p:q ) in the global matrix

    Class( real_distributed_matrix ), Intent( In    ) :: matrix
    Integer                         , Intent( In    ) :: m
    Integer                         , Intent( In    ) :: n
    Integer                         , Intent( In    ) :: p
    Integer                         , Intent( In    ) :: q
    Real( wp ), Dimension( m:, p: ) , Intent(   Out ) :: data

    Integer :: i_glob, j_glob
    Integer :: i_loc , j_loc
    Integer :: handle
    Integer :: error
    
    ! THIS NEEDS OPTIMISATION!!
    data = 0.0_wp
    Do j_glob = p, q
       j_loc = matrix%global_to_local_cols( j_glob )
       If( j_loc == distributed_matrix_NOT_ME ) Cycle
       Do i_glob = m, n
          i_loc = matrix%global_to_local_rows( i_glob )
          If( i_loc == distributed_matrix_NOT_ME ) Cycle
          data( i_glob, j_glob ) = matrix%data( i_loc, j_loc )
       End Do
    End Do
    ! Generate a portable MPI data type handle from the variable to be communicated
    Call MPI_Type_create_f90_real( Precision( data ), Range( data ), handle, error )
    ! Replicate the data
    Call MPI_Allreduce( MPI_IN_PLACE, data, Size( data ), handle, MPI_SUM, matrix%matrix_map%get_comm(), error )
       
  End Subroutine matrix_get_global_real

  Subroutine matrix_get_local_real( matrix, m, n, p, q, data )

    ! Gets the data ( m:n, p:q ) in the local matrix

    Class( real_distributed_matrix ), Intent( In    ) :: matrix
    Integer                         , Intent( In    ) :: m
    Integer                         , Intent( In    ) :: n
    Integer                         , Intent( In    ) :: p
    Integer                         , Intent( In    ) :: q
    Real( wp ), Dimension( m:, p: ) , Intent(   Out ) :: data

    data( m:n, p:q ) = matrix%data( m:n, p:q )
       
  End Subroutine matrix_get_local_real

  Subroutine matrix_get_global_complex( matrix, m, n, p, q, data )

    ! Gets the data ( m:n, p:q ) in the global matrix

    Class( complex_distributed_matrix ), Intent( In    ) :: matrix
    Integer                            , Intent( In    ) :: m
    Integer                            , Intent( In    ) :: n
    Integer                            , Intent( In    ) :: p
    Integer                            , Intent( In    ) :: q
    Complex( wp ), Dimension( m:, p: ) , Intent(   Out ) :: data

    Integer :: i_glob, j_glob
    Integer :: i_loc , j_loc
    Integer :: handle
    Integer :: error
    
    ! THIS NEEDS OPTIMISATION!!
    data = 0.0_wp
    Do j_glob = p, q
       j_loc = matrix%global_to_local_cols( j_glob )
       If( j_loc == distributed_matrix_NOT_ME ) Cycle
       Do i_glob = m, n
          i_loc = matrix%global_to_local_rows( i_glob )
          If( i_loc == distributed_matrix_NOT_ME ) Cycle
          data( i_glob, j_glob ) = matrix%data( i_loc, j_loc )
       End Do
    End Do
    ! Generate a portable MPI data type handle from the variable to be communicated
    Call MPI_Type_create_f90_complex( Precision( data ), Range( data ), handle, error )
    ! Replicate the data
    Call MPI_Allreduce( MPI_IN_PLACE, data, Size( data ), handle, MPI_SUM, matrix%matrix_map%get_comm(), error )
       
  End Subroutine matrix_get_global_complex

  Subroutine matrix_get_local_complex( matrix, m, n, p, q, data )

    ! Gets the data ( m:n, p:q ) in the local matrix

    Class( complex_distributed_matrix ), Intent( In    ) :: matrix
    Integer                            , Intent( In    ) :: m
    Integer                            , Intent( In    ) :: n
    Integer                            , Intent( In    ) :: p
    Integer                            , Intent( In    ) :: q
    Complex( wp ), Dimension( m:, p: ) , Intent(   Out ) :: data

    data( m:n, p:q ) = matrix%data( m:n, p:q )
       
  End Subroutine matrix_get_local_complex
  
  Subroutine matrix_get_maps( matrix, gl_rows, gl_cols, lg_rows, lg_cols )

    Class( distributed_matrix )         , Intent( In    ) :: matrix
    Integer, Dimension( : ), Allocatable, Intent(   Out ) :: gl_rows
    Integer, Dimension( : ), Allocatable, Intent(   Out ) :: gl_cols
    Integer, Dimension( : ), Allocatable, Intent(   Out ) :: lg_rows
    Integer, Dimension( : ), Allocatable, Intent(   Out ) :: lg_cols

    ! Note using allocate on set
    gl_rows = matrix%global_to_local_rows
    gl_cols = matrix%global_to_local_cols
    lg_rows = matrix%local_to_global_rows
    lg_cols = matrix%local_to_global_cols
    
  End Subroutine matrix_get_maps

  Subroutine matrix_diag_real( A, Q, E )

    Use numbers_module, Only : wp

    Implicit None

    Class( real_distributed_matrix ),              Intent( In    ) :: A
    Type ( real_distributed_matrix ), Allocatable, Intent(   Out ) :: Q
    Real( wp ), Dimension( : )      , Allocatable, Intent(   Out ) :: E

    Real( wp ), Dimension( :, : ), Allocatable :: tmp_a

    Real( wp ), Dimension( : ), Allocatable :: work

    Integer, Dimension( : ), Allocatable :: iwork
    
    Integer :: nwork
    Integer :: npcol
    Integer :: m, n
    Integer :: info

    ! Give Q the same mapping as A
    Q = A
    
    Call A%matrix_map%get_data( m = m, n = n, npcol = npcol )

    Allocate( E( 1:m ) )
    
    ! The diag overwrites the matrix. Horrible so use a temporary
    tmp_A = A%data
    
    ! Workspace size enquiry
    Allocate( work( 1:1 ), iwork( 1:1 ) )
    Call pdsyevd( 'V', 'U', m, tmp_A, 1, 1, A%matrix_map%get_descriptor(), E, Q%data, 1, 1, Q%matrix_map%get_descriptor(), &
         work, -1, iwork, 0, info )
    nwork = Nint( work( 1 ) )
    nwork = nwork * diag_work_size_fiddle_factor ! From experience ...
    Deallocate( work, iwork )
    Allocate(  work( 1:nwork ) )
    ! Scalapack recipe is behind the strange numbers
    Allocate( iwork( 1:7 * m + 8 * npcol + 2 ) )
    ! Do the diag
    Call pdsyevd( 'V', 'U', m, tmp_A, 1, 1, A%matrix_map%get_descriptor(), E, Q%data, 1, 1, Q%matrix_map%get_descriptor(), &
         work, Size( work ), iwork, Size( iwork ), info )

  End Subroutine matrix_diag_real

  Subroutine matrix_diag_complex( A, Q, E )

    Use numbers_module, Only : wp

    Implicit None

    Class( complex_distributed_matrix ),              Intent( In    ) :: A
    Type ( complex_distributed_matrix ), Allocatable, Intent(   Out ) :: Q
    Real( wp ), Dimension( : )         , Allocatable, Intent(   Out ) :: E

    Complex( wp ), Dimension( :, : ), Allocatable :: tmp_a

    Complex( wp ), Dimension( : ), Allocatable :: cwork

    Real( wp ), Dimension( : ), Allocatable :: rwork

    Integer, Dimension( : ), Allocatable :: iwork
    
    Integer :: ncwork, nrwork
    Integer :: npcol
    Integer :: m, n
    Integer :: info

    ! Give Q the same mapping as A
    Q = A
    
    Call A%matrix_map%get_data( m = m, n = n, npcol = npcol )

    Allocate( E( 1:m ) )

    ! The diag overwrites the matrix. Horrible so use a temporary
    tmp_A = A%data
       
    ! Workspace size enquiry
    Allocate( cwork( 1:1 ), rwork( 1:1 ), iwork( 1:1 ) )
    Call pzheevd( 'V', 'U', m, tmp_A, 1, 1, A%matrix_map%get_descriptor(), E, Q%data, 1, 1, Q%matrix_map%get_descriptor(), &
         cwork, -1, rwork, -1, iwork, 0, info )
    ncwork = Nint( Real( cwork( 1 ), wp ) )
    ncwork = ncwork * diag_work_size_fiddle_factor ! From experience ...
    nrwork = Nint( rwork( 1 ) )
    nrwork = nrwork * diag_work_size_fiddle_factor ! From experience ...
    Deallocate( cwork, rwork, iwork )
    Allocate( cwork( 1:ncwork ) )
    Allocate( rwork( 1:nrwork ) )
    ! Scalapack recipe is behind the strange numbers
    Allocate( iwork( 1:7 * m + 8 * npcol + 2 ) )
    ! Do the diag
    Call pzheevd( 'V', 'U', m, tmp_A, 1, 1, A%matrix_map%get_descriptor(), E, Q%data, 1, 1, Q%matrix_map%get_descriptor(), &
            cwork, Size( cwork ), rwork, Size( rwork ), iwork, Size( iwork ), info )

  End Subroutine matrix_diag_complex

  Function matrix_multiply_real( A, B ) Result( C )

    Class( real_distributed_matrix ), Allocatable :: C

    Class( real_distributed_matrix ), Intent( In ) :: A
    Class( real_distributed_matrix ), Intent( In ) :: B

    Integer :: ma, na
    Integer :: mb, nb
    Integer :: m, n, k

    Character :: t1, t2

    ! Give C the same mapping as A
    Allocate( C, Source = A )

    ! There must be a neater way ...
    Deallocate( C%data )
    Deallocate( C%local_to_global_rows )
    Deallocate( C%local_to_global_cols )
    Deallocate( C%global_to_local_rows )
    Deallocate( C%global_to_local_cols )
    C%daggered = .False.
    
    t1 = Merge( 'T', 'N', A%daggered )
    t2 = Merge( 'T', 'N', B%daggered )
    
    Call A%matrix_map%get_data( m = ma, n = na )
    Call B%matrix_map%get_data( m = mb, n = nb )
    
    If( t1 == 'N' .And. t2 == 'N' ) Then
       m = ma
       n = nb
       k = na
    Else If( t1 == 'T' .And. t2 == 'N' ) Then
       m = na
       n = nb
       k = ma
    Else If( t1 == 'N' .And. t2 == 'T' ) Then
       m = ma
       n = mb
       k = na
    Else If( t1 == 'T' .And. t2 == 'T' ) Then
       m = na
       n = mb
       k = ma
    Else
       Stop 'How did we get here in matrix_multiply_real???'
    End If
    
    Call matrix_create( C, m, n, A )
    
    Call pdgemm( t1, t2, m, n, k, 1.0_wp, A%data, 1, 1, A%matrix_map%get_descriptor(), &
                                          B%data, 1, 1, B%matrix_map%get_descriptor(), &
                                  0.0_wp, C%data, 1, 1, C%matrix_map%get_descriptor() )
          
  End Function matrix_multiply_real
     
  Function matrix_multiply_complex( A, B ) Result( C )

    Class( complex_distributed_matrix ), Allocatable :: C

    Class( complex_distributed_matrix ), Intent( In ) :: A
    Class( complex_distributed_matrix ), Intent( In ) :: B

    Integer :: ma, na
    Integer :: mb, nb
    Integer :: m, n, k

    Character :: t1, t2

    ! Give C the same mapping as A
    Allocate( C, Source = A )

    ! There must be a neater way ...
    Deallocate( C%data )
    Deallocate( C%local_to_global_rows )
    Deallocate( C%local_to_global_cols )
    Deallocate( C%global_to_local_rows )
    Deallocate( C%global_to_local_cols )
    C%daggered = .False.
    
    t1 = Merge( 'C', 'N', A%daggered )
    t2 = Merge( 'C', 'N', B%daggered )
       
    Call A%matrix_map%get_data( m = ma, n = na )
    Call B%matrix_map%get_data( m = mb, n = nb )

    If( t1 == 'N' .And. t2 == 'N' ) Then
       m = ma
       n = nb
       k = na
    Else If( t1 == 'C' .And. t2 == 'N' ) Then
       m = na
       n = nb
       k = ma
    Else If( t1 == 'N' .And. t2 == 'C' ) Then
       m = ma
       n = mb
       k = na
    Else If( t1 == 'C' .And. t2 == 'C' ) Then
       m = na
       n = mb
       k = ma
    Else
       Stop 'How did we get here in matrix_multiply_complex???'
    End If
    
    Call matrix_create( C, m, n, A )

    Call pzgemm( t1, t2, m, n, k, ( 1.0_wp, 0.0_wp ), A%data, 1, 1, A%matrix_map%get_descriptor(), &
                                                      B%data, 1, 1, B%matrix_map%get_descriptor(), &
                                  ( 0.0_wp, 0.0_wp ), C%data, 1, 1, C%matrix_map%get_descriptor() )

  End Function matrix_multiply_complex

  Subroutine matrix_extract_cols_real( A, c1, c2, B )

    ! Needs extending!!! Only does first c2 cols currently
    ! Also need to think about transposes - probably easiest once equivalent row functionality
    ! implemented
    
    Class( real_distributed_matrix ), Intent( In    ) :: A
    Integer                         , Intent( In    ) :: c1 
    Integer                         , Intent( In    ) :: c2
    Type ( real_distributed_matrix ), Intent(   Out ) :: B

    Integer :: m
    Integer :: c
    Integer :: c1_loc, c2_loc

    Call A%matrix_map%get_data( m = m )
    Call matrix_create( B, m, c2 - c1 + 1, A )

    If( c1 /= 1 ) Then
       Stop "General extraction of cols, not implemented"
    End If
    c1_loc = 1

    ! Map global c2 to a local index
    c = c2
    c2_loc = 0
    Do While( c > 0 )
       If( A%global_to_local_cols( c ) /= DISTRIBUTED_MATRIX_NOT_ME ) Then
          c2_loc = A%global_to_local_cols( c )
          Exit
       Else
          c = c - 1
       End If
    End Do

    B%data = A%data( :, c1_loc:c2_loc )

  End Subroutine matrix_extract_cols_real

  Subroutine matrix_extract_cols_complex( A, c1, c2, B )

    ! Needs extending!!! Only does first c2 cols currently
    ! Also need to think about transposes - probably easiest once equivalent row functionality
    ! implemented
    
    Class( complex_distributed_matrix ), Intent( In    ) :: A
    Integer                            , Intent( In    ) :: c1 
    Integer                            , Intent( In    ) :: c2
    Type ( complex_distributed_matrix ), Intent(   Out ) :: B

    Integer :: m
    Integer :: c
    Integer :: c1_loc, c2_loc

    Call A%matrix_map%get_data( m = m )
    Call matrix_create( B, m, c2 - c1 + 1, A )
    
    If( c1 /= 1 ) Then
       Stop "General extraction of cols, not implemented"
    End If
    c1_loc = 1

    ! Map global c2 to a local index
    c = c2
    c2_loc = 0
    Do While( c > 0 )
       If( A%global_to_local_cols( c ) /= DISTRIBUTED_MATRIX_NOT_ME ) Then
          c2_loc = A%global_to_local_cols( c )
          Exit
       Else
          c = c - 1
       End If
    End Do

    B%data = A%data( :, c1_loc:c2_loc )

  End Subroutine matrix_extract_cols_complex

  Function  matrix_post_scale_real( A, s ) Result( B )

    Class( real_distributed_matrix ), Allocatable :: B

    Class( real_distributed_matrix ), Intent( In ) :: A
    Real( wp )                      , Intent( In ) :: s

    Allocate( B, Source = A )
    B%data = s * A%data
    
  End Function matrix_post_scale_real

  Function  matrix_pre_scale_real( s, A ) Result( B )

    Class( real_distributed_matrix ), Allocatable :: B

    Real( wp )                      , Intent( In ) :: s
    Class( real_distributed_matrix ), Intent( In ) :: A

    Allocate( B, Source = A )
    B%data = s * A%data
    
  End Function matrix_pre_scale_real

  Function  matrix_post_scale_complex( A, s ) Result( B )

    Class( complex_distributed_matrix ), Allocatable :: B

    Class( complex_distributed_matrix ), Intent( In ) :: A
    Complex( wp )                      , Intent( In ) :: s

    Allocate( B, Source = A )
    B%data = s * A%data
    
  End Function matrix_post_scale_complex

  Function  matrix_pre_scale_complex( s, A ) Result( B )

    Class( complex_distributed_matrix ), Allocatable :: B

    Complex( wp )                      , Intent( In ) :: s
    Class( complex_distributed_matrix ), Intent( In ) :: A

    Allocate( B, Source = A )
    B%data = s * A%data
    
  End Function matrix_pre_scale_complex

  Function matrix_add_real( A, B ) Result( C )

    ! Note in an effort to avoid communication through transposes the addition occurs in
    ! the form with A NOT transposed, and then the result is indicate as requiring transposition
    ! or not as required byt this.

    Class( real_distributed_matrix ), Allocatable :: C

    Class( real_distributed_matrix ), Intent( In ) :: A
    Class( real_distributed_matrix ), Intent( In ) :: B

    Integer :: m, n
    
    Character :: tA, tB

    tA = Merge( 'T', 'N', A%daggered )
    tB = Merge( 'T', 'N', B%daggered )
    Call A%matrix_map%get_data( m = m, n = n )
    Allocate( real_distributed_matrix :: C )
    Call matrix_create( C, m, n, A )
    C%data = B%data
    Call pdgeadd( tB, m, n, 1.0_wp, A%data, 1, 1, A%matrix_map%get_descriptor(), &
                            1.0_wp, C%data, 1, 1, C%matrix_map%get_descriptor() )
  
    Select Case( tA )
    Case Default
       Stop "How did we get here in matrix_add_real"
    Case( "N" )
       C%daggered = .False.
    Case( "T" )
       C%daggered = .True.
    End Select
              
  End Function matrix_add_real
     
  Function matrix_add_complex( A, B ) Result( C )

    ! Note in an effort to avoid communication through transposes the addition occurs in
    ! the form with A NOT transposed, and then the result is indicate as requiring transposition
    ! or not as required byt this.

    Class( complex_distributed_matrix ), Allocatable :: C

    Class( complex_distributed_matrix ), Intent( In ) :: A
    Class( complex_distributed_matrix ), Intent( In ) :: B

    Integer :: m, n
    
    Character :: tA, tB

    tA = Merge( 'T', 'N', A%daggered )
    tB = Merge( 'T', 'N', B%daggered )
    Call A%matrix_map%get_data( m = m, n = n )
    Allocate( complex_distributed_matrix :: C )
    Call matrix_create( C, m, n, A )
    C%data = B%data
    Call pzgeadd( tB, m, n, ( 1.0_wp, 0.0_wp ), A%data, 1, 1, A%matrix_map%get_descriptor(), &
                            ( 1.0_wp, 0.0_wp ), C%data, 1, 1, C%matrix_map%get_descriptor() )
  
    Select Case( tA )
    Case Default
       Stop "How did we get here in matrix_add_complex"
    Case( "N" )
       C%daggered = .False.
    Case( "T" )
       C%daggered = .True.
    End Select
              
  End Function matrix_add_complex
     
!!$  Subroutine dummy( A )
!!$    Class( distributed_matrix ), Intent( In ) :: A
!!$    Stop "Should never get here"
!!$    Write( *, * ) A%daggered
!!$  End Subroutine dummy
!!$  
!!$  Logical Function dummy_f( A, rubbish )
!!$    Class( distributed_matrix ), Intent( In ) :: A
!!$    Integer                    , Intent( In ) :: rubbish
!!$    Stop "Should never get here"
!!$    dummy_f = A%daggered
!!$    Write( *, * ) rubbish
!!$  End Function dummy_f
  
End Module distributed_matrix_module
 
