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
     Procedure          :: create          => matrix_create
     Procedure          :: get_maps        => matrix_get_maps
     Procedure          :: global_to_local => matrix_global_to_local
     Procedure          :: local_to_global => matrix_local_to_global
     Procedure          :: size            => matrix_size
     Procedure          :: local_size      => matrix_local_size
     Procedure          :: get_comm        => matrix_communicator
  End type distributed_matrix

  Type, Extends( distributed_matrix ), Public :: real_distributed_matrix
     Real( wp ), Dimension( :, : ), Allocatable :: data
   Contains
     Procedure, Private   :: diag_r               => matrix_diag_real
     Generic              :: diag                 => diag_r
     Procedure, Private   :: dagger_r             => matrix_dagger_real
     Generic              :: dagger               => dagger_r
     Generic              :: Operator( .Dagger. ) => dagger_r
     Procedure            :: multiply             => matrix_multiply_real
     Procedure, Pass( A ) :: pre_scale            => matrix_pre_scale_real
     Procedure            :: post_scale           => matrix_post_scale_real
     Procedure, Pass( A ) :: pre_mult_diag        => matrix_pre_mult_diag_real
     Procedure            :: post_mult_diag       => matrix_post_mult_diag_real
     Generic              :: Operator( * )        => multiply, pre_scale, post_scale
     Generic              :: Operator( * )        => pre_mult_diag, post_mult_diag
     Procedure            :: add                  => matrix_add_real
     Procedure            :: post_add_diag        => matrix_post_add_diag_real
     Procedure, Pass( A ) :: pre_add_diag         => matrix_pre_add_diag_real
     Generic              :: Operator( + )        => add, post_add_diag, pre_add_diag
     Procedure            :: subtract             => matrix_subtract_real
     Procedure            :: subtract_diag        => matrix_post_subtract_diag_real
     Generic              :: Operator( - )        => subtract, subtract_diag
     Procedure            :: Choleski             => matrix_choleski_real
     Procedure            :: Solve                => matrix_solve_real
     Procedure            :: set_to_identity      => matrix_set_to_identity_real
     Procedure, Private   :: set_by_global_r      => matrix_set_global_real
     Procedure, Private   :: set_by_local_r       => matrix_set_local_real
     Procedure, Private   :: get_by_global_r      => matrix_get_global_real
     Procedure, Private   :: get_by_local_r       => matrix_get_local_real
     Generic              :: set_by_global        => set_by_global_r
     Generic              :: set_by_local         => set_by_local_r
     Generic              :: get_by_global        => get_by_global_r
     Generic              :: get_by_local         => get_by_local_r
     Procedure, Private   :: extract_r            => matrix_extract_real
     Generic              :: extract              => extract_r
  End type real_distributed_matrix

  Type, Extends( distributed_matrix ), Public :: complex_distributed_matrix
     Complex( wp ), Dimension( :, : ), Allocatable :: data
   Contains
     Procedure, Private   :: diag_c               => matrix_diag_complex
     Generic              :: diag                 => diag_c
     Procedure, Private   :: dagger_c             => matrix_dagger_complex
     Generic              :: dagger               => dagger_c 
     Generic              :: Operator( .Dagger. ) => dagger_c
     Procedure            :: multiply             => matrix_multiply_complex
     Procedure, Pass( A ) :: pre_scale            => matrix_pre_scale_complex
     Procedure            :: post_scale           => matrix_post_scale_complex
     Procedure, Pass( A ) :: pre_mult_diag        => matrix_pre_mult_diag_complex
     Procedure            :: post_mult_diag       => matrix_post_mult_diag_complex
     Generic              :: Operator( * )        => multiply, pre_scale, post_scale
     Generic              :: Operator( * )        => pre_mult_diag, post_mult_diag
     Procedure            :: add                  => matrix_add_complex
     Procedure            :: post_add_diag        => matrix_post_add_diag_complex
     Procedure, Pass( A ) :: pre_add_diag         => matrix_pre_add_diag_complex
     Generic              :: Operator( + )        => add, post_add_diag, pre_add_diag
     Procedure            :: subtract             => matrix_subtract_complex
     Procedure            :: subtract_diag        => matrix_post_subtract_diag_complex
     Generic              :: Operator( - )        => subtract, subtract_diag
     Procedure            :: Choleski             => matrix_choleski_complex
     Procedure            :: Solve                => matrix_solve_complex
     Procedure            :: set_to_identity      => matrix_set_to_identity_complex
     Procedure, Private   :: set_by_global_c      => matrix_set_global_complex
     Procedure, Private   :: set_by_local_c       => matrix_set_local_complex
     Procedure, Private   :: get_by_global_c      => matrix_get_global_complex
     Procedure, Private   :: get_by_local_c       => matrix_get_local_complex
     Generic              :: set_by_global        => set_by_global_c
     Generic              :: set_by_local         => set_by_local_c
     Generic              :: get_by_global        => get_by_global_c
     Generic              :: get_by_local         => get_by_local_c
     Procedure, Private   :: extract_c            => matrix_extract_complex
     Generic              :: extract              => extract_c
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

    Use, Intrinsic :: ieee_arithmetic

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

    Call matrix%matrix_map%set( source_matrix%matrix_map%proc_mapping, ctxt, m, n, mb, nb, 0, 0, lda )

    Call set_local_to_global( matrix%local_to_global_rows, m, mb, myprow, nprow, lda )
    Call set_local_to_global( matrix%local_to_global_cols, n, nb, mypcol, npcol, sda )

    Call set_global_to_local( matrix%global_to_local_rows, m, mb, myprow, nprow )
    Call set_global_to_local( matrix%global_to_local_cols, n, nb, mypcol, npcol )

    Select Type( matrix )
    Class Default
       Stop "Illegal type in matrix_create"
    Class is ( real_distributed_matrix )
       Allocate( matrix%data( 1:lda, 1:sda  ) )
       matrix%data = ieee_value( matrix%data, ieee_signaling_nan )
    Class is ( complex_distributed_matrix )
       Allocate( matrix%data( 1:lda, 1:sda  ) )
       ! No signalling nan for complex
       matrix%data = Cmplx( ieee_value( 0.0_wp, ieee_signaling_nan ), &
            ieee_value( 0.0_wp, ieee_signaling_nan ), wp )
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

  Function matrix_choleski_real( A ) Result( C )

    Class( real_distributed_matrix ), Allocatable :: C

    Class( real_distributed_matrix ), Intent( In ) :: A

    Integer :: m
    Integer :: i_glob, j_glob
    Integer :: i, j
    Integer :: error
    
    Allocate( C, Source = A )

    ! Transpose don't matter as A must be symmetric. So set it as untransposed
    ! so I don't get confused
    C%daggered = .False.

    ! Zero Upper half of C
    Do j = 1, Size( C%data, Dim = 2 )
       j_glob = C%local_to_global_cols( j )
       Do i = 1, Size( C%data, Dim = 1 )
          i_glob = C%local_to_global_rows( i )
          If( j_glob > i_glob ) Then
             C%data( i, j ) = 0.0_wp
          End If
       End Do
    End Do

    Call C%matrix_map%get_data( m = m )
    Call pdpotrf( 'L', m, C%data, 1, 1, C%matrix_map%get_descriptor(), error )
    If( error /= 0 ) Then
       Deallocate( C )
    End If
    
  End Function matrix_choleski_real

  Function matrix_choleski_complex( A ) Result( C )

    Class( complex_distributed_matrix ), Allocatable :: C

    Class( complex_distributed_matrix ), Intent( In ) :: A

    Integer :: m
    Integer :: i_glob, j_glob
    Integer :: i, j
    Integer :: error
    
    Allocate( C, Source = A )

    ! Transpose don't matter as A must be symmetric. So set it as untransposed
    ! so I don't get confused
    C%daggered = .False.

    ! Zero Upper half of C
    Do j = 1, Size( C%data, Dim = 2 )
       j_glob = C%local_to_global_cols( j )
       Do i = 1, Size( C%data, Dim = 1 )
          i_glob = C%local_to_global_rows( i )
          If( j_glob > i_glob ) Then
             C%data( i, j ) = 0.0_wp
          End If
       End Do
    End Do

    Call C%matrix_map%get_data( m = m )
    Call pzpotrf( 'L', m, C%data, 1, 1, C%matrix_map%get_descriptor(), error )
    
    If( error /= 0 ) Then
       Deallocate( C )
    End If
    
  End Function matrix_choleski_complex

  Function matrix_solve_real( A, B ) Result( C )

    ! Need to think tranposes!!!

    Class( real_distributed_matrix ), Allocatable :: C

    Class( real_distributed_matrix ), Intent( In ) :: A
    Class( real_distributed_matrix ), Intent( In ) :: B

    Class( real_distributed_matrix ), Allocatable :: R

    Integer, Dimension( : ), Allocatable :: pivots

    Integer :: m, mb, nrhs
    Integer :: error

    ! Otherwise A, B overwritten by pdgesv, Grrrrrr
    Allocate( R, Source = A )
    Allocate( C, Source = B )

    Call R%matrix_map%get_data( m = m, mb = mb )
    Call C%matrix_map%get_data( n = nrhs )
    Allocate( pivots( 1:m + mb ) )
    Call pdgesv( m, nrhs, R%data, 1, 1, R%matrix_map%get_descriptor(), pivots, &
                          C%data, 1, 1, C%matrix_map%get_descriptor(), error )
    If( error /= 0 ) Then
       Deallocate( C )
    End If
    
  End Function matrix_solve_real

  Function matrix_solve_complex( A, B ) Result( C )

    ! Need to think tranposes!!!

    Class( complex_distributed_matrix ), Allocatable :: C

    Class( complex_distributed_matrix ), Intent( In ) :: A
    Class( complex_distributed_matrix ), Intent( In ) :: B

    Class( complex_distributed_matrix ), Allocatable :: R

    Integer, Dimension( : ), Allocatable :: pivots

    Integer :: m, mb, nrhs
    Integer :: error

    ! Otherwise A, B overwritten by pdgesv, Grrrrrr
    Allocate( R, Source = A )
    Allocate( C, Source = B )

    Call R%matrix_map%get_data( m = m, mb = mb )
    Call C%matrix_map%get_data( n = nrhs )
    Allocate( pivots( 1:m + mb ) )
    Call pzgesv( m, nrhs, R%data, 1, 1, R%matrix_map%get_descriptor(), pivots, &
                          C%data, 1, 1, C%matrix_map%get_descriptor(), error )
    If( error /= 0 ) Then
       Deallocate( C )
    End If
    
  End Function matrix_solve_complex

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

    Use mpi
    
    ! Gets the data ( m:n, p:q ) in the global matrix

    Class( real_distributed_matrix ), Intent( In    ) :: matrix
    Integer                         , Intent( In    ) :: m
    Integer                         , Intent( In    ) :: n
    Integer                         , Intent( In    ) :: p
    Integer                         , Intent( In    ) :: q
    Real( wp ), Dimension( m:, p: ) , Intent(   Out ) :: data
    
    Real( wp ) :: rdum

    Integer :: i_glob, j_glob
    Integer :: i_loc , j_loc
    Integer :: handle
    Integer :: rsize, error
    
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
!!!!HACK TO WORK AROUND BUG IN MVAPICH2
!!$    Call MPI_Allreduce( MPI_IN_PLACE, data, Size( data ), handle, MPI_SUM, matrix%matrix_map%get_comm(), error )
    Call mpi_sizeof( rdum, rsize, error )
    Call mpi_type_match_size( MPI_TYPECLASS_REAL, rsize, handle, error )
    Call MPI_Allreduce( MPI_IN_PLACE, data, Size( data ), MPI_DOUBLE_PRECISION, MPI_SUM, matrix%matrix_map%get_comm(), error )
       
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

    Real( wp ) :: cdum

    Integer :: i_glob, j_glob
    Integer :: i_loc , j_loc
    Integer :: csize, handle
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
!!!!HACK TO WORK AROUND BUG IN MVAPICH2
!!$    Call MPI_Allreduce( MPI_IN_PLACE, data, Size( data ), handle, MPI_SUM, matrix%matrix_map%get_comm(), error )
    Call mpi_sizeof( cdum, csize, error )
    Call mpi_type_match_size( MPI_TYPECLASS_REAL, csize, handle, error )
    Call MPI_Allreduce( MPI_IN_PLACE, data, Size( data ), MPI_double_complex, MPI_SUM, matrix%matrix_map%get_comm(), error )
       
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
!!$    Allocate( Q, Source = A )
    
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

    If( info /= 0 ) Then
       Deallocate( Q )
       Deallocate( E )
    End If

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
!!$    Allocate( Q, Source = A )
    
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

    If( info /= 0 ) Then
       Deallocate( Q )
       Deallocate( E )
    End If

  End Subroutine matrix_diag_complex

  Function matrix_multiply_real( A, B ) Result( C )

    Use mpi

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

  Function matrix_post_mult_diag_real( A, d ) Result( B )

    Class( real_distributed_matrix ), Allocatable :: B

    Class( real_distributed_matrix ),                 Intent( In ) :: A
    Real( wp )                      , Dimension( : ), Intent( In ) :: d

    Integer :: ma, na
    Integer :: j_glob
    Integer :: j_loc

    ! TRANSPOSES!!!
    
    Call A%matrix_map%get_data( m = ma, n = na )

    Select Case( A%daggered )
    Case( .False. )
       ! Error check
       If( na == Size( d ) ) Then
          Allocate( B, Source = A )
          ! Check MATHS!!!
          Do j_loc = 1, Size( B%data, Dim = 2 )
             j_glob = B%local_to_global_cols( j_loc )
             B%data( :, j_loc ) = B%data( :, j_loc ) * d( j_glob )
          End Do
       End If
    Case( .True. )
       Stop "mult diag real not implemented tranposes"
    End Select
    
  End Function matrix_post_mult_diag_real

  Function matrix_post_mult_diag_complex( A, d ) Result( B )

    Class( complex_distributed_matrix ), Allocatable :: B

    Class( complex_distributed_matrix ),                 Intent( In ) :: A
    Complex( wp )                      , Dimension( : ), Intent( In ) :: d

    Integer :: ma, na
    Integer :: j_glob
    Integer :: j_loc

    ! TRANSPOSES!!!
    
    Call A%matrix_map%get_data( m = ma, n = na )

    Select Case( A%daggered )
    Case( .False. )
       ! Error check
       If( na == Size( d ) ) Then
          Allocate( B, Source = A )
          ! Check MATHS!!!
          Do j_loc = 1, Size( B%data, Dim = 2 )
             j_glob = B%local_to_global_cols( j_loc )
             B%data( :, j_loc ) = B%data( :, j_loc ) * d( j_glob )
          End Do
       End If
    Case( .True. )
       Stop "mult diag complex not implemented tranposes"
    End Select
    
  End Function matrix_post_mult_diag_complex

  Function matrix_pre_mult_diag_real( d, A ) Result( B )

    Class( real_distributed_matrix ), Allocatable :: B

    Real( wp )                      , Dimension( : ), Intent( In ) :: d
    Class( real_distributed_matrix ),                 Intent( In ) :: A

    Integer :: ma, na
    Integer :: i_glob
    Integer :: i_loc

    ! TRANSPOSES!!!
    
    Call A%matrix_map%get_data( m = ma, n = na )

    Select Case( A%daggered )
    Case( .False. )
       ! Error check
       If( ma == Size( d ) ) Then
          Allocate( B, Source = A )
          ! Check MATHS!!!
          Do i_loc = 1, Size( B%data, Dim = 1 )
             i_glob = B%local_to_global_rows( i_loc )
             B%data( i_loc, : ) = B%data( i_loc, : ) * d( i_glob )
          End Do
       End If
    Case( .True. )
       Stop "mult diag real not implemented tranposes"
    End Select
    
  End Function matrix_pre_mult_diag_real

  Function matrix_pre_mult_diag_complex( d, A ) Result( B )

    Class( complex_distributed_matrix ), Allocatable :: B

    Complex( wp )                      , Dimension( : ), Intent( In ) :: d
    Class( complex_distributed_matrix ),                 Intent( In ) :: A

    Integer :: ma, na
    Integer :: i_glob
    Integer :: i_loc

    ! TRANSPOSES!!!
    
    Call A%matrix_map%get_data( m = ma, n = na )

    Select Case( A%daggered )
    Case( .False. )
       ! Error check
       If( ma == Size( d ) ) Then
          Allocate( B, Source = A )
          ! Check MATHS!!!
          Do i_loc = 1, Size( B%data, Dim = 1 )
             i_glob = B%local_to_global_rows( i_loc )
             B%data( i_loc, : ) = B%data( i_loc, : ) * d( i_glob )
          End Do
       End If
    Case( .True. )
       Stop "mult diag complex not implemented tranposes"
    End Select
    
  End Function matrix_pre_mult_diag_complex

  Function matrix_extract_real( A, r1, r2, c1, c2 ) Result( B )

    Class( real_distributed_matrix ), Allocatable :: B

    ! ALSO NEED TO THINK ABOUT TRANSPOSES
    
    Class( real_distributed_matrix ), Intent( In    ) :: A
    Integer                         , Intent( In    ) :: r1 
    Integer                         , Intent( In    ) :: r2
    Integer                         , Intent( In    ) :: c1 
    Integer                         , Intent( In    ) :: c2

    Integer :: mb
    Integer :: nb
    Integer :: a_ctxt

    Allocate( real_distributed_matrix :: B )

    mb = r2 - r1 + 1
    nb = c2 - c1 + 1
    Call matrix_create( B, mb, nb, A )
    !!!TRANSPOSES!!!! 
    B%daggered = A%daggered

    Call A%matrix_map%get_data( ctxt = a_ctxt )
    Call pdgemr2d( mb, nb, A%data, r1, c1, A%matrix_map%get_descriptor(), &
                           B%data,  1,  1, B%matrix_map%get_descriptor(), a_ctxt )

  End Function matrix_extract_real

  Function matrix_extract_complex( A, r1, r2, c1, c2 ) Result( B )

    Class( complex_distributed_matrix ), Allocatable :: B

    ! ALSO NEED TO THINK ABOUT TRANSPOSES
    
    Class( complex_distributed_matrix ), Intent( In    ) :: A
    Integer                         , Intent( In    ) :: r1 
    Integer                         , Intent( In    ) :: r2
    Integer                         , Intent( In    ) :: c1 
    Integer                         , Intent( In    ) :: c2

    Integer :: mb
    Integer :: nb
    Integer :: a_ctxt

    Allocate( complex_distributed_matrix :: B )

    mb = r2 - r1 + 1
    nb = c2 - c1 + 1
    Call matrix_create( B, mb, nb, A )
    !!!TRANSPOSES!!!! 
    B%daggered = A%daggered
    
    Call A%matrix_map%get_data( ctxt = a_ctxt )
    Call pzgemr2d( mb, nb, A%data, r1, c1, A%matrix_map%get_descriptor(), &
                           B%data,  1,  1, B%matrix_map%get_descriptor(), a_ctxt )

  End Function matrix_extract_complex

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
     
  Function matrix_post_add_diag_real( A, d ) Result( B )

    Class( real_distributed_matrix ), Allocatable :: B

    Class( real_distributed_matrix ),                 Intent( In ) :: A
    Real( wp )                      , Dimension( : ), Intent( In ) :: d

    Integer :: m, n
    Integer :: i_glob
    Integer :: i_loc, j_loc
    
    Call A%matrix_map%get_data( m = m, n = n )

    If( m == n .And. Size( d ) == n ) Then
       Allocate( B, Source = A )
       Do i_glob = 1, n
          i_loc = A%global_to_local_rows( i_glob )
          j_loc = A%global_to_local_cols( i_glob )
          If(  i_loc /= distributed_matrix_NOT_ME .And. &
               j_loc /= distributed_matrix_NOT_ME ) Then
             B%data( i_loc, j_loc ) = A%data( i_loc, j_loc ) + d( i_glob )
          End If
       End Do
    End If
    
  End Function matrix_post_add_diag_real

  Function  matrix_post_add_diag_complex( A, d ) Result( B )

    Class( complex_distributed_matrix ), Allocatable :: B

    Class( complex_distributed_matrix ),                 Intent( In ) :: A
    Complex( wp )                      , Dimension( : ), Intent( In ) :: d

    Integer :: m, n
    Integer :: i_glob
    Integer :: i_loc, j_loc
    
    Call A%matrix_map%get_data( m = m, n = n )

    If( m == n .And. Size( d ) == n ) Then
       Allocate( B, Source = A )
       Do i_glob = 1, n
          i_loc = A%global_to_local_rows( i_glob )
          j_loc = A%global_to_local_cols( i_glob )
          If(  i_loc /= distributed_matrix_NOT_ME .And. &
               j_loc /= distributed_matrix_NOT_ME ) Then
             B%data( i_loc, j_loc ) = A%data( i_loc, j_loc ) + d( i_glob )
          End If
       End Do
    End If
    
  End Function matrix_post_add_diag_complex

  Function matrix_pre_add_diag_real( d, A ) Result( B )

    Class( real_distributed_matrix ), Allocatable :: B

    Real( wp )                      , Dimension( : ), Intent( In ) :: d
    Class( real_distributed_matrix ),                 Intent( In ) :: A


    Integer :: m, n
    Integer :: i_glob
    Integer :: i_loc, j_loc
    
    Call A%matrix_map%get_data( m = m, n = n )

    If( m == n .And. Size( d ) == n ) Then
       Allocate( B, Source = A )
       Do i_glob = 1, n
          i_loc = A%global_to_local_rows( i_glob )
          j_loc = A%global_to_local_cols( i_glob )
          If(  i_loc /= distributed_matrix_NOT_ME .And. &
               j_loc /= distributed_matrix_NOT_ME ) Then
             B%data( i_loc, j_loc ) = A%data( i_loc, j_loc ) + d( i_glob )
          End If
       End Do
    End If
    
  End Function matrix_pre_add_diag_real

  Function matrix_pre_add_diag_complex( d, A ) Result( B )

    Class( complex_distributed_matrix ), Allocatable :: B

    Complex( wp )                      , Dimension( : ), Intent( In ) :: d
    Class( complex_distributed_matrix ),                 Intent( In ) :: A

    Integer :: m, n
    Integer :: i_glob
    Integer :: i_loc, j_loc
    
    Call A%matrix_map%get_data( m = m, n = n )

    If( m == n .And. Size( d ) == n ) Then
       Allocate( B, Source = A )
       Do i_glob = 1, n
          i_loc = A%global_to_local_rows( i_glob )
          j_loc = A%global_to_local_cols( i_glob )
          If(  i_loc /= distributed_matrix_NOT_ME .And. &
               j_loc /= distributed_matrix_NOT_ME ) Then
             B%data( i_loc, j_loc ) = A%data( i_loc, j_loc ) + d( i_glob )
          End If
       End Do
    End If
    
  End Function matrix_pre_add_diag_complex

  Function matrix_subtract_real( A, B ) Result( C )

    ! Note in an effort to avoid communication through transposes the subtraction occurs in
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
    Call pdgeadd( tB, m, n,   1.0_wp, A%data, 1, 1, A%matrix_map%get_descriptor(), &
                            - 1.0_wp, C%data, 1, 1, C%matrix_map%get_descriptor() )
  
    Select Case( tA )
    Case Default
       Stop "How did we get here in matrix_subtract_real"
    Case( "N" )
       C%daggered = .False.
    Case( "T" )
       C%daggered = .True.
    End Select
              
  End Function matrix_subtract_real
     
  Function matrix_subtract_complex( A, B ) Result( C )

    ! Note in an effort to avoid communication through transposes the subtraction occurs in
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
    Call pzgeadd( tB, m, n, (   1.0_wp, 0.0_wp ), A%data, 1, 1, A%matrix_map%get_descriptor(), &
                            ( - 1.0_wp, 0.0_wp ), C%data, 1, 1, C%matrix_map%get_descriptor() )
  
    Select Case( tA )
    Case Default
       Stop "How did we get here in matrix_subtract_complex"
    Case( "N" )
       C%daggered = .False.
    Case( "T" )
       C%daggered = .True.
    End Select
              
  End Function matrix_subtract_complex
     
  Function matrix_post_subtract_diag_real( A, d ) Result( B )

    Class( real_distributed_matrix ), Allocatable :: B

    Class( real_distributed_matrix ),                 Intent( In ) :: A
    Real( wp )                      , Dimension( : ), Intent( In ) :: d

    Integer :: m, n
    Integer :: i_glob
    Integer :: i_loc, j_loc
    
    Call A%matrix_map%get_data( m = m, n = n )

    If( m == n .And. Size( d ) == n ) Then
       Allocate( B, Source = A )
       Do i_glob = 1, n
          i_loc = A%global_to_local_rows( i_glob )
          j_loc = A%global_to_local_cols( i_glob )
          If(  i_loc /= distributed_matrix_NOT_ME .And. &
               j_loc /= distributed_matrix_NOT_ME ) Then
             B%data( i_loc, j_loc ) = A%data( i_loc, j_loc ) + d( i_glob )
          End If
       End Do
    End If
    
  End Function matrix_post_subtract_diag_real

  Function  matrix_post_subtract_diag_complex( A, d ) Result( B )

    Class( complex_distributed_matrix ), Allocatable :: B

    Class( complex_distributed_matrix ),                 Intent( In ) :: A
    Complex( wp )                      , Dimension( : ), Intent( In ) :: d

    Integer :: m, n
    Integer :: i_glob
    Integer :: i_loc, j_loc
    
    Call A%matrix_map%get_data( m = m, n = n )

    If( m == n .And. Size( d ) == n ) Then
       Allocate( B, Source = A )
       Do i_glob = 1, n
          i_loc = A%global_to_local_rows( i_glob )
          j_loc = A%global_to_local_cols( i_glob )
          If(  i_loc /= distributed_matrix_NOT_ME .And. &
               j_loc /= distributed_matrix_NOT_ME ) Then
             B%data( i_loc, j_loc ) = A%data( i_loc, j_loc ) - d( i_glob )
          End If
       End Do
    End If
    
  End Function matrix_post_subtract_diag_complex

  Subroutine matrix_set_to_identity_real( A ) 

    Class( real_distributed_matrix ), Intent( InOut ) :: A

    Integer :: m
    Integer :: i_glob
    Integer :: i_loc, j_loc
    
    A%data = 0.0_wp

    Call A%matrix_map%get_data( m = m )
    Do i_glob = 1, m
       i_loc = A%global_to_local_rows( i_glob )
       j_loc = A%global_to_local_cols( i_glob )
       If(  i_loc /= distributed_matrix_NOT_ME .And. &
            j_loc /= distributed_matrix_NOT_ME ) Then
          A%data( i_loc, j_loc ) = 1.0_wp
       End If
    End Do

  End Subroutine matrix_set_to_identity_real

  Subroutine matrix_set_to_identity_complex( A ) 

    Class( complex_distributed_matrix ), Intent( InOut ) :: A

    Integer :: m
    Integer :: i_glob
    Integer :: i_loc, j_loc
    
    A%data = 0.0_wp

    Call A%matrix_map%get_data( m = m )
    Do i_glob = 1, m
       i_loc = A%global_to_local_rows( i_glob )
       j_loc = A%global_to_local_cols( i_glob )
       If(  i_loc /= distributed_matrix_NOT_ME .And. &
            j_loc /= distributed_matrix_NOT_ME ) Then
          A%data( i_loc, j_loc ) = 1.0_wp
       End If
    End Do

  End Subroutine matrix_set_to_identity_complex

  Function matrix_global_to_local( A, what ) Result( gl_indexing )

    Integer, Dimension( : ), Allocatable :: gl_indexing
    
    Class( distributed_matrix ), Intent( In ) :: A
    Character( Len = * )       , Intent( In ) :: what

    Select Case( what )
    Case Default
       Stop "Illegal WHAT in global_to_local"
    Case( 'R', 'r' )
       gl_indexing = A%global_to_local_rows
    Case( 'C', 'c' )
       gl_indexing = A%global_to_local_cols
    End Select

  End Function matrix_global_to_local
  
  Function matrix_local_to_global( A, what ) Result( lg_indexing )

    Integer, Dimension( : ), Allocatable :: lg_indexing
    
    Class( distributed_matrix ), Intent( In ) :: A
    Character( Len = * )       , Intent( In ) :: what

    Select Case( what )
    Case Default
       Stop "Illegal WHAT in local_to_global"
    Case( 'R', 'r' )
       lg_indexing = A%local_to_global_rows
    Case( 'C', 'c' )
       lg_indexing = A%local_to_global_cols
    End Select

  End Function matrix_local_to_global

  Function matrix_local_size( A, dim ) Result( n )

    Integer :: n

    Class( distributed_matrix ), Intent( In ) :: A
    Integer                    , Intent( In ) :: dim

    If( dim <= 2 ) Then

       Select Type( A )
       Class Default
          Stop "Illegal type in matrix_local_size"
       Class is ( real_distributed_matrix )
          n = Size( A%data, Dim = dim )
       Class is ( complex_distributed_matrix )
          n = Size( A%data, Dim = dim )
       End Select

    Else

       Stop "Ilegal dimension in matrix_local_size"

    End If
       
  End Function matrix_local_size
  
  Function matrix_size( A, dim ) Result( n )

    Integer :: n

    Class( distributed_matrix ), Intent( In ) :: A
    Integer                    , Intent( In ) :: dim

    Select Case( dim )
    Case Default
       Stop "Ilegal dimension in matrix_size"
    Case( 1 )
       Call A%matrix_map%get_data( m = n )
    Case( 2 )
       Call A%matrix_map%get_data( n = n )
    End Select
       
  End Function matrix_size
  
  Function matrix_communicator( A ) Result( c )

    Integer :: c

    Class( distributed_matrix ), Intent( In ) :: A

    c = A%matrix_map%get_comm()
    
  End Function matrix_communicator
  
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
 
