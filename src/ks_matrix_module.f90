Module ks_matrix_module

  ! Wrapper for dist matrix types so can make opaque array - only 1 k point at this level

  Use numbers_module           , Only : wp
  Use distributed_matrix_module, Only : distributed_matrix, real_distributed_matrix, complex_distributed_matrix, &
       distributed_matrix_init, distributed_matrix_comm_to_base, distributed_matrix_finalise, &
       distributed_matrix_remap_data
  Implicit None

  Type, Private :: k_point_matrix
     ! Should drop this spin and this_k_point in final implementation and
     ! only store in level above. But does no major harm so leave as going to
     ! rewrite anyway
     Integer                                 , Private :: this_spin
     Integer, Dimension( 1:3 )               , Private :: this_k_point
     Class( distributed_matrix ), Allocatable, Private :: matrix
  End Type k_point_matrix

  Type, Public :: ks_matrix
     Type( k_point_matrix ), Allocatable, Private :: k_point
   Contains
     ! Public methods
     Procedure                     :: create               => ks_matrix_create
     Generic                       :: Operator( .Dagger. ) => dagger
     Procedure                     :: diag                 => ks_matrix_diag
     Generic                       :: Operator( * )        => multiply
     Generic                       :: Operator( * )        => post_scale, pre_scale
     Generic                       :: Operator( * )        => pre_mult_diag, post_mult_diag
     Generic                       :: Operator( + )        => add, post_add_diag, pre_add_diag
     Generic                       :: Operator( - )        => subtract
     Generic                       :: Operator( - )        => post_subtract_diag
     Procedure                     :: Choleski             => ks_matrix_Choleski
     Procedure                     :: Solve                => ks_matrix_Solve
     Procedure                     :: set_to_identity      => ks_matrix_set_to_identity
     Generic                       :: set_by_global        => sgr, sgc
     Generic                       :: set_by_local         => slr, slc
     Generic                       :: get_by_global        => ggr, ggc
     Generic                       :: get_by_local         => glr, glc
     Procedure                     :: extract              => ks_matrix_extract
     Procedure                     :: global_to_local      => ks_matrix_g_to_l
     Procedure                     :: local_to_global      => ks_matrix_l_to_g
     Procedure                     :: local_size           => ks_matrix_local_size
     Procedure                     :: size                 => ks_matrix_size
     Procedure                     :: get_comm             => ks_matrix_get_communicator
     ! Private implementations
     Procedure, Private            :: dagger               => ks_matrix_dagger
     Procedure, Private            :: multiply             => ks_matrix_mult
     Procedure, Private            :: post_scale           => ks_matrix_post_scale
     Procedure, Private, Pass( A ) :: pre_scale            => ks_matrix_pre_scale
     Procedure, Private            :: post_mult_diag       => ks_matrix_post_mult_diag
     Procedure, Private, Pass( A ) :: pre_mult_diag        => ks_matrix_pre_mult_diag
     Procedure, Private            :: add                  => ks_matrix_add
     Procedure, Private, Pass( A ) :: pre_add_diag         => ks_matrix_pre_add_diag
     Procedure, Private            :: post_add_diag        => ks_matrix_post_add_diag
     Procedure, Private            :: subtract             => ks_matrix_subtract
     Procedure, Private            :: post_subtract_diag   => ks_matrix_post_subtract_diag
     Procedure, Private            :: sgr                  => set_global_real
     Procedure, Private            :: sgc                  => set_global_complex
     Procedure, Private            :: slr                  => set_local_real
     Procedure, Private            :: slc                  => set_local_complex
     Procedure, Private            :: ggr                  => get_global_real
     Procedure, Private            :: ggc                  => get_global_complex
     Procedure, Private            :: glr                  => get_local_real
     Procedure, Private            :: glc                  => get_local_complex
  End Type ks_matrix

  ! Method for k point //ism
  ! a) Set up base matrix which is 1 matrix across all
  ! b) Provide splits

  Public :: ks_matrix_init
  Public :: ks_matrix_comm_to_base
  Public :: ks_matrix_finalise
  Public :: ks_matrix_remap_data
  
  Private

  Integer, Parameter :: INVALID = -1

Contains

  Subroutine ks_matrix_init
    Call distributed_matrix_init
  End Subroutine ks_matrix_init

  Subroutine ks_matrix_comm_to_base( comm, base_matrix )

    Integer                        , Intent( In    ) :: comm
    Type   ( ks_matrix ), Intent(   Out ) :: base_matrix

    Allocate( k_point_matrix:: base_matrix%k_point )

    base_matrix%k_point%this_spin    = INVALID
    base_matrix%k_point%this_k_point = INVALID

    Allocate( real_distributed_matrix:: base_matrix%k_point%matrix )

    Call distributed_matrix_comm_to_base( comm, base_matrix%k_point%matrix )
    
  End Subroutine ks_matrix_comm_to_base

  Subroutine ks_matrix_finalise

    Call distributed_matrix_finalise
    
  End Subroutine ks_matrix_finalise
  
  Subroutine ks_matrix_create( matrix, is_complex, this_spin, this_k_point, m, n, source_matrix )

    Class  ( ks_matrix ), Intent(   Out ) :: matrix
    Logical                        , Intent( In    ) :: is_complex
    Integer                        , Intent( In    ) :: this_spin
    Integer, Dimension( 1:3 )      , Intent( In    ) :: this_k_point
    Integer                        , Intent( In    ) :: m
    Integer                        , Intent( In    ) :: n
    Type   ( ks_matrix ), Intent( In    ) :: source_matrix

    Allocate( k_point_matrix:: matrix%k_point )

    If( is_complex ) Then
      Allocate( complex_distributed_matrix :: matrix%k_point%matrix ) 
    Else
      Allocate( real_distributed_matrix    :: matrix%k_point%matrix ) 
    End If

    matrix%k_point%this_spin    = this_spin
    matrix%k_point%this_k_point = this_k_point

    Call matrix%k_point%matrix%create( m, n, source_matrix%k_point%matrix )
    
  End Subroutine ks_matrix_create

  Subroutine ks_matrix_diag( A, Q, evals )
    
    Class( ks_matrix )          , Intent( In    ) :: A
    Type ( ks_matrix )          , Intent(   Out ) :: Q
    Real( wp ), Dimension( : ), Allocatable, Intent(   Out ) :: evals

    Type(    real_distributed_matrix ), Allocatable :: Q_real
    Type( complex_distributed_matrix ), Allocatable :: Q_complex
    
    Allocate( Q%k_point )

    Q%k_point%this_spin    = A%k_point%this_spin
    Q%k_point%this_k_point = A%k_point%this_k_point

    Associate( Akm => A%k_point%matrix )
    
      Select Type( Akm )

      Class Default
         Stop "Illegal type in ks_matrix_diag"

      Class is ( real_distributed_matrix )
         Call Akm%diag( Q_real, evals )
         Allocate( Q%k_point%matrix, Source = Q_real )

      Class is ( complex_distributed_matrix )
         Call Akm%diag( Q_complex, evals )
         Allocate( Q%k_point%matrix, Source = Q_complex )

      End Select

    End Associate

  End Subroutine ks_matrix_diag
  
  Function ks_matrix_dagger( A ) Result( tA )
    
    Type( ks_matrix ), Allocatable :: tA

    Class( ks_matrix ), Intent( In ) :: A

    Type(    real_distributed_matrix ) :: tA_real
    Type( complex_distributed_matrix ) :: tA_complex

    Allocate( tA )
    Allocate( tA%k_point )
    
    tA%k_point%this_spin    = A%k_point%this_spin
    tA%k_point%this_k_point = A%k_point%this_k_point

    Associate( Akm => A%k_point%matrix )
      Select Type( Akm )
      Class Default
         Stop "Illegal type in ks_matrix_dagger"
      Type is ( real_distributed_matrix )
         tA_real = .Dagger. Akm
         Allocate( tA%k_point%matrix, Source = tA_real )
      Type is ( complex_distributed_matrix )
         tA_complex = .Dagger. Akm
         Allocate( tA%k_point%matrix, Source = tA_complex )
      End Select
         
    End Associate

  End Function ks_matrix_dagger
    
  Function ks_matrix_Choleski( A ) Result( C )
    
    Type( ks_matrix ), Allocatable :: C

    Class( ks_matrix ), Intent( In ) :: A

    Type(    real_distributed_matrix ) :: C_real
    Type( complex_distributed_matrix ) :: C_complex

    Allocate( C )
    Allocate( C%k_point )
    
    C%k_point%this_spin    = A%k_point%this_spin
    C%k_point%this_k_point = A%k_point%this_k_point

    Associate( Akm => A%k_point%matrix )
      Select Type( Akm )
      Class Default
         Stop "Illegal type in ks_matrix_Choleski"
      Type is ( real_distributed_matrix )
         C_real = Akm%Choleski()
         Allocate( C%k_point%matrix, Source = C_real )
      Type is ( complex_distributed_matrix )
         C_complex = Akm%Choleski()
         Allocate( C%k_point%matrix, Source = C_complex )
      End Select
         
    End Associate

  End Function ks_matrix_Choleski
    
  Function ks_matrix_mult( A, B ) Result( C )
    
    Type( ks_matrix ), Allocatable :: C

    Class( ks_matrix ), Intent( In ) :: A
    Type ( ks_matrix ), Intent( In ) :: B

    Type(    real_distributed_matrix ) :: C_real
    Type( complex_distributed_matrix ) :: C_complex

    Allocate( C )
    Allocate( C%k_point )

    C%k_point%this_spin    = A%k_point%this_spin
    C%k_point%this_k_point = A%k_point%this_k_point

    Associate( Akm => A%k_point%matrix )
    
      Select Type( Akm )

      Class Default
         Stop "Illegal type in ks_matrix_mult"
         
      Type is ( real_distributed_matrix )
         
         Associate( Bkm => B%k_point%matrix )
           Select Type( Bkm )
           Class Default
              Stop "Illegal type in ks_matrix_mult"
           Type is ( real_distributed_matrix )
              C_real = Akm * Bkm
              Allocate( C%k_point%matrix, Source = C_real )
           End Select
         End Associate

      Type is ( complex_distributed_matrix )
         Associate( Bkm => B%k_point%matrix )
           Select Type( Bkm )
           Class Default
              Stop "Illegal type in ks_matrix_mult"
           Type is ( complex_distributed_matrix )
              C_complex = Akm * Bkm
              Allocate( C%k_point%matrix, Source = C_complex )
           End Select
         End Associate
         
      End Select
    End Associate

  End Function ks_matrix_mult

  Function ks_matrix_post_mult_diag( A, d ) Result( B )
    
    Type( ks_matrix ), Allocatable :: B

    Class( ks_matrix ), Intent( In ) :: A
    Real ( wp ), Dimension( : )  , Intent( In ) :: d

    Type(    real_distributed_matrix ) :: B_real
    Type( complex_distributed_matrix ) :: B_complex

    Allocate( B )
    Allocate( B%k_point )
    
    B%k_point%this_spin    = A%k_point%this_spin
    B%k_point%this_k_point = A%k_point%this_k_point

    Associate( Akm => A%k_point%matrix )
      Select Type( Akm )
      Class Default
         Stop "Illegal type in ks_matrix_post_mult_diag"
      Type is ( real_distributed_matrix )
         B_real = Akm * d
         Allocate( B%k_point%matrix, Source = B_real )
      Type is ( complex_distributed_matrix )
         B_complex = Akm * Cmplx( d, Kind = wp )
         Allocate( B%k_point%matrix, Source = B_complex )
      End Select
    End Associate

  End Function ks_matrix_post_mult_diag

  Function ks_matrix_pre_mult_diag( d, A ) Result( B )
    
    Type( ks_matrix ), Allocatable :: B

    Real ( wp ), Dimension( : )  , Intent( In ) :: d
    Class( ks_matrix ), Intent( In ) :: A

    Type(    real_distributed_matrix ) :: B_real
    Type( complex_distributed_matrix ) :: B_complex

    Allocate( B )
    Allocate( B%k_point )
    
    B%k_point%this_spin    = A%k_point%this_spin
    B%k_point%this_k_point = A%k_point%this_k_point

    Associate( Akm => A%k_point%matrix )
      Select Type( Akm )
      Class Default
         Stop "Illegal type in ks_matrix_pre_mult_diag"
      Type is ( real_distributed_matrix )
         B_real = d * Akm 
         Allocate( B%k_point%matrix, Source = B_real )
      Type is ( complex_distributed_matrix )
         B_complex = Cmplx( d, Kind = wp ) * Akm
         Allocate( B%k_point%matrix, Source = B_complex )
      End Select
    End Associate

  End Function ks_matrix_pre_mult_diag

  Function ks_matrix_solve( A, B ) Result( C )
    
    Type( ks_matrix ), Allocatable :: C

    Class( ks_matrix ), Intent( In ) :: A
    Type ( ks_matrix ), Intent( In ) :: B

    Type(    real_distributed_matrix ) :: C_real
    Type( complex_distributed_matrix ) :: C_complex

    Allocate( C )
    Allocate( k_point_matrix :: C%k_point )

    C%k_point%this_spin    = A%k_point%this_spin
    C%k_point%this_k_point = A%k_point%this_k_point

    Associate( Akm => A%k_point%matrix )
    
      Select Type( Akm )

      Class Default
         Stop "Illegal type in ks_matrix_solve"
         
      Type is ( real_distributed_matrix )
         
         Associate( Bkm => B%k_point%matrix )
           Select Type( Bkm )
           Class Default
              Stop "Illegal type in ks_matrix_solve"
           Type is ( real_distributed_matrix )
              C_real = Akm%solve( Bkm )
              Allocate( C%k_point%matrix, Source = C_real )
           End Select
         End Associate

      Type is ( complex_distributed_matrix )
         Associate( Bkm => B%k_point%matrix )
           Select Type( Bkm )
           Class Default
              Stop "Illegal type in ks_matrix_solve"
           Type is ( complex_distributed_matrix )
              C_complex = Akm%solve( Bkm )
              Allocate( C%k_point%matrix, Source = C_complex )
           End Select
         End Associate
         
      End Select
    End Associate

  End Function ks_matrix_solve

  Function ks_matrix_add( A, B ) Result( C )
    
    Type( ks_matrix ), Allocatable :: C

    Class( ks_matrix ), Intent( In ) :: A
    Type ( ks_matrix ), Intent( In ) :: B

    Type(    real_distributed_matrix ) :: C_real
    Type( complex_distributed_matrix ) :: C_complex

    Allocate( C )
    Allocate( C%k_point )

    C%k_point%this_spin    = A%k_point%this_spin
    C%k_point%this_k_point = A%k_point%this_k_point

    Associate( Akm => A%k_point%matrix )
    
      Select Type( Akm )

      Class Default
         Stop "Illegal type in ks_matrix_add"
         
      Type is ( real_distributed_matrix )
         
         Associate( Bkm => B%k_point%matrix )
           Select Type( Bkm )
           Class Default
              Stop "Illegal type in ks_matrix_add"
           Type is ( real_distributed_matrix )
              C_real = Akm + Bkm
              Allocate( C%k_point%matrix, Source = C_real )
           End Select
         End Associate

      Type is ( complex_distributed_matrix )
         Associate( Bkm => B%k_point%matrix )
           Select Type( Bkm )
           Class Default
              Stop "Illegal type in ks_matrix_add"
           Type is ( complex_distributed_matrix )
              C_complex = Akm + Bkm
              Allocate( C%k_point%matrix, Source = C_complex )
           End Select
         End Associate
         
      End Select
    End Associate

  End Function ks_matrix_add

  Function ks_matrix_subtract( A, B ) Result( C )
    
    Type( ks_matrix ), Allocatable :: C

    Class( ks_matrix ), Intent( In ) :: A
    Type ( ks_matrix ), Intent( In ) :: B

    Type(    real_distributed_matrix ) :: C_real
    Type( complex_distributed_matrix ) :: C_complex

    Allocate( C )
    Allocate( C%k_point )

    C%k_point%this_spin    = A%k_point%this_spin
    C%k_point%this_k_point = A%k_point%this_k_point

    Associate( Akm => A%k_point%matrix )
    
      Select Type( Akm )

      Class Default
         Stop "Illegal type in ks_matrix_subtract"
         
      Type is ( real_distributed_matrix )
         
         Associate( Bkm => B%k_point%matrix )
           Select Type( Bkm )
           Class Default
              Stop "Illegal type in ks_matrix_subtract"
           Type is ( real_distributed_matrix )
              C_real = Akm - Bkm
              Allocate( C%k_point%matrix, Source = C_real )
           End Select
         End Associate

      Type is ( complex_distributed_matrix )
         Associate( Bkm => B%k_point%matrix )
           Select Type( Bkm )
           Class Default
              Stop "Illegal type in ks_matrix_subtract"
           Type is ( complex_distributed_matrix )
              C_complex = Akm - Bkm
              Allocate( C%k_point%matrix, Source = C_complex )
           End Select
         End Associate
         
      End Select
    End Associate

  End Function ks_matrix_subtract

  Function ks_matrix_post_subtract_diag( A, d ) Result( B )
    
    Type( ks_matrix ), Allocatable :: B

    Class( ks_matrix ), Intent( In ) :: A
    Real ( wp ), Dimension( : )  , Intent( In ) :: d

    Type(    real_distributed_matrix ) :: B_real
    Type( complex_distributed_matrix ) :: B_complex

    Allocate( B )
    Allocate( B%k_point )
    
    B%k_point%this_spin    = A%k_point%this_spin
    B%k_point%this_k_point = A%k_point%this_k_point

    Associate( Akm => A%k_point%matrix )
      Select Type( Akm )
      Class Default
         Stop "Illegal type in ks_matrix_post_scale"
      Type is ( real_distributed_matrix )
         B_real = Akm - d
         Allocate( B%k_point%matrix, Source = B_real )
      Type is ( complex_distributed_matrix )
         B_complex = Akm - Cmplx( d, Kind = wp )
         Allocate( B%k_point%matrix, Source = B_complex )
      End Select
    End Associate

  End Function ks_matrix_post_subtract_diag

  Subroutine set_global_real( A, m, n, p, q, data )
    
    Class( ks_matrix )   , Intent( InOut ) :: A
    Integer                         , Intent( In    ) :: m
    Integer                         , Intent( In    ) :: n
    Integer                         , Intent( In    ) :: p
    Integer                         , Intent( In    ) :: q
    Real( wp ), Dimension( m:, p: ) , Intent( In    ) :: data

    Associate( Akm => A%k_point%matrix )
    
      Select Type( Akm )

      Class Default
         Stop "Illegal type in set_global_real"
         
      Type is ( real_distributed_matrix )
         Call Akm%set_by_global( m, n, p, q, data )

      End Select
      
    End Associate
    
  End Subroutine set_global_real
    
  Subroutine set_global_complex( A, m, n, p, q, data )
    
    Class( ks_matrix )      , Intent( InOut ) :: A
    Integer                            , Intent( In    ) :: m
    Integer                            , Intent( In    ) :: n
    Integer                            , Intent( In    ) :: p
    Integer                            , Intent( In    ) :: q
    Complex( wp ), Dimension( m:, p: ) , Intent( In    ) :: data

    Associate( Akm => A%k_point%matrix )
    
      Select Type( Akm )

      Class Default
         Stop "Illegal type in set_global_complex"
         
      Type is ( complex_distributed_matrix )
         Call Akm%set_by_global( m, n, p, q, data )

      End Select
      
    End Associate
    
  End Subroutine set_global_complex
    
  Subroutine set_local_real( A, m, n, p, q, data )
    
    Class( ks_matrix )   , Intent( InOut ) :: A
    Integer                         , Intent( In    ) :: m
    Integer                         , Intent( In    ) :: n
    Integer                         , Intent( In    ) :: p
    Integer                         , Intent( In    ) :: q
    Real( wp ), Dimension( m:, p: ) , Intent( In    ) :: data

    Associate( Akm => A%k_point%matrix )
    
      Select Type( Akm )

      Class Default
         Stop "Illegal type in set_local_real"
         
      Type is ( real_distributed_matrix )
         Call Akm%set_by_local( m, n, p, q, data )

      End Select
      
    End Associate
    
  End Subroutine set_local_real
    
  Subroutine set_local_complex( A, m, n, p, q, data )
    
    Class( ks_matrix )      , Intent( InOut ) :: A
    Integer                            , Intent( In    ) :: m
    Integer                            , Intent( In    ) :: n
    Integer                            , Intent( In    ) :: p
    Integer                            , Intent( In    ) :: q
    Complex( wp ), Dimension( m:, p: ) , Intent( In    ) :: data

    Associate( Akm => A%k_point%matrix )
    
      Select Type( Akm )

      Class Default
         Stop "Illegal type in set_local_complex"
         
      Type is ( complex_distributed_matrix )
         Call Akm%set_by_local( m, n, p, q, data )

      End Select
      
    End Associate
    
  End Subroutine set_local_complex
    
  Subroutine get_global_real( A, m, n, p, q, data )
    
    Class( ks_matrix )   , Intent( In    ) :: A
    Integer                         , Intent( In    ) :: m
    Integer                         , Intent( In    ) :: n
    Integer                         , Intent( In    ) :: p
    Integer                         , Intent( In    ) :: q
    Real( wp ), Dimension( m:, p: ) , Intent(   Out ) :: data

    Associate( Akm => A%k_point%matrix )
    
      Select Type( Akm )

      Class Default
         Stop "Illegal type in get_global_real"
         
      Type is ( real_distributed_matrix )
         Call Akm%get_by_global( m, n, p, q, data )

      End Select
      
    End Associate
    
  End Subroutine get_global_real
    
  Subroutine get_global_complex( A, m, n, p, q, data )
    
    Class( ks_matrix )      , Intent( In    ) :: A
    Integer                            , Intent( In    ) :: m
    Integer                            , Intent( In    ) :: n
    Integer                            , Intent( In    ) :: p
    Integer                            , Intent( In    ) :: q
    Complex( wp ), Dimension( m:, p: ) , Intent(   Out ) :: data

    Associate( Akm => A%k_point%matrix )
    
      Select Type( Akm )

      Class Default
         Stop "Illegal type in get_global_complex"
         
      Type is ( complex_distributed_matrix )
         Call Akm%get_by_global( m, n, p, q, data )

      End Select
      
    End Associate
    
  End Subroutine get_global_complex
    
  Subroutine get_local_real( A, m, n, p, q, data )
    
    Class( ks_matrix )   , Intent( In    ) :: A
    Integer                         , Intent( In    ) :: m
    Integer                         , Intent( In    ) :: n
    Integer                         , Intent( In    ) :: p
    Integer                         , Intent( In    ) :: q
    Real( wp ), Dimension( m:, p: ) , Intent(   Out ) :: data

    Associate( Akm => A%k_point%matrix )
    
      Select Type( Akm )

      Class Default
         Stop "Illegal type in get_local_real"
         
      Type is ( real_distributed_matrix )
         Call Akm%get_by_local( m, n, p, q, data )

      End Select
      
    End Associate
    
  End Subroutine get_local_real
    
  Subroutine get_local_complex( A, m, n, p, q, data )
    
    Class( ks_matrix )      , Intent( In    ) :: A
    Integer                            , Intent( In    ) :: m
    Integer                            , Intent( In    ) :: n
    Integer                            , Intent( In    ) :: p
    Integer                            , Intent( In    ) :: q
    Complex( wp ), Dimension( m:, p: ) , Intent(   Out ) :: data

    Associate( Akm => A%k_point%matrix )
    
      Select Type( Akm )

      Class Default
         Stop "Illegal type in get_local_complex"
         
      Type is ( complex_distributed_matrix )
         Call Akm%get_by_local( m, n, p, q, data )

      End Select
      
    End Associate
    
  End Subroutine get_local_complex

  Function ks_matrix_extract( A, r1, r2, c1, c2 ) Result( B )
    
    Type ( ks_matrix ), Allocatable :: B

    Class( ks_matrix ), Intent( In    ) :: A
    Integer                      , Intent( In    ) :: r1 
    Integer                      , Intent( In    ) :: r2
    Integer                      , Intent( In    ) :: c1 
    Integer                      , Intent( In    ) :: c2

    Type(    real_distributed_matrix ) :: B_real
    Type( complex_distributed_matrix ) :: B_complex

    Allocate( B )
    Allocate( B%k_point )
    B%k_point%this_spin    = A%k_point%this_spin
    B%k_point%this_k_point = A%k_point%this_k_point

    Associate( Akm => A%k_point%matrix )
    
      Select Type( Akm )

      Class Default
         Stop "Illegal type in ks_matrix_extract"
         
      Type is ( real_distributed_matrix )
         B_real = Akm%extract( r1, r2, c1, c2 )
         Allocate( B%k_point%matrix, Source = B_real )

      Type is ( complex_distributed_matrix )
         B_complex = Akm%extract( r1, r2, c1, c2 )
         Allocate( B%k_point%matrix, Source = B_complex )

      End Select
    End Associate
         
  End Function ks_matrix_extract
  
  Function ks_matrix_post_scale( A, s ) Result( B )
    
    Type( ks_matrix ), Allocatable :: B

    Class( ks_matrix ), Intent( In ) :: A
    Real ( wp )                  , Intent( In ) :: s

    Type(    real_distributed_matrix ) :: B_real
    Type( complex_distributed_matrix ) :: B_complex

    Allocate( B )
    Allocate( B%k_point )
    
    B%k_point%this_spin    = A%k_point%this_spin
    B%k_point%this_k_point = A%k_point%this_k_point

    Associate( Akm => A%k_point%matrix )
      Select Type( Akm )
      Class Default
         Stop "Illegal type in ks_matrix_post_scale"
      Type is ( real_distributed_matrix )
         B_real = s * Akm
         Allocate( B%k_point%matrix, Source = B_real )
      Type is ( complex_distributed_matrix )
         B_complex = Cmplx( s, Kind = wp ) * Akm
         Allocate( B%k_point%matrix, Source = B_complex )
      End Select
    End Associate

  End Function ks_matrix_post_scale

  Function ks_matrix_pre_scale( s, A ) Result( B )
    
    Type( ks_matrix ), Allocatable :: B

    Real ( wp )                  , Intent( In ) :: s
    Class( ks_matrix ), Intent( In ) :: A

    Type(    real_distributed_matrix ) :: B_real
    Type( complex_distributed_matrix ) :: B_complex

    Allocate( B )
    Allocate( B%k_point )
    
    B%k_point%this_spin    = A%k_point%this_spin
    B%k_point%this_k_point = A%k_point%this_k_point

    Associate( Akm => A%k_point%matrix )
      Select Type( Akm )
      Class Default
         Stop "Illegal type in ks_matrix_pre_scale"
      Type is ( real_distributed_matrix )
         B_real = s * Akm
         Allocate( B%k_point%matrix, Source = B_real )
      Type is ( complex_distributed_matrix )
         B_complex = Cmplx( s, Kind = wp ) * Akm
         Allocate( B%k_point%matrix, Source = B_complex )
      End Select
    End Associate

  End Function ks_matrix_pre_scale

  Function ks_matrix_post_add_diag( A, d ) Result( B )
    
    Type( ks_matrix ), Allocatable :: B

    Class( ks_matrix ), Intent( In ) :: A
    Real ( wp ), Dimension( : )  , Intent( In ) :: d

    Type(    real_distributed_matrix ) :: B_real
    Type( complex_distributed_matrix ) :: B_complex

    Allocate( B )
    Allocate( B%k_point )
    
    B%k_point%this_spin    = A%k_point%this_spin
    B%k_point%this_k_point = A%k_point%this_k_point

    Associate( Akm => A%k_point%matrix )
      Select Type( Akm )
      Class Default
         Stop "Illegal type in ks_matrix_post_add_diag"
      Type is ( real_distributed_matrix )
         B_real = Akm + d
         Allocate( B%k_point%matrix, Source = B_real )
      Type is ( complex_distributed_matrix )
         B_complex = Akm + Cmplx( d, Kind = wp )
         Allocate( B%k_point%matrix, Source = B_complex )
      End Select
    End Associate

  End Function ks_matrix_post_add_diag

  Function ks_matrix_pre_add_diag( d, A ) Result( B )
    
    Type( ks_matrix ), Allocatable :: B

    Real ( wp ), Dimension( : )  , Intent( In ) :: d
    Class( ks_matrix ), Intent( In ) :: A

    Type(    real_distributed_matrix ) :: B_real
    Type( complex_distributed_matrix ) :: B_complex

    Allocate( B )
    Allocate( B%k_point )
    
    B%k_point%this_spin    = A%k_point%this_spin
    B%k_point%this_k_point = A%k_point%this_k_point

    Associate( Akm => A%k_point%matrix )
      Select Type( Akm )
      Class Default
         Stop "Illegal type in ks_matrix_pre_add_diag"
      Type is ( real_distributed_matrix )
         B_real = Akm + d
         Allocate( B%k_point%matrix, Source = B_real )
      Type is ( complex_distributed_matrix )
         B_complex = Akm + Cmplx( d, Kind = wp )
         Allocate( B%k_point%matrix, Source = B_complex )
      End Select
    End Associate

  End Function ks_matrix_pre_add_diag

  Subroutine ks_matrix_remap_data( parent_comm, A, B )

    ! Horrible routine to remap the data part of the matrix
    ! Horrible becuase we know nothing about overlap of the two sets of processors
    ! which hold the matrices.
    ! Parent_comm is a communicator holding the union of processes that hold A and will hold B
    ! A is the matrix from which we are remapping
    ! B is the matrix to which we are remapping
    ! There is a horrible complication that not all processes required to call this neccessarily hold
    ! any part of A or B - thy may hold only one (and the routine requires they hold at least parts of
    ! one). To indicate A this process does not hold part of the matrix the actual argument
    ! should be in a deallocated state

    ! NOTE THIS DOES NO MATRIX SETTING UP - all that happens is the data in the matrix is redistributed
    
    Integer                     ,              Intent( In    ) :: parent_comm
    Type( ks_matrix ), Allocatable, Intent( In    ) :: A
    Type( ks_matrix ), Allocatable, Intent( InOut ) :: B

    Type(    real_distributed_matrix ), Allocatable :: A_mat_real   , B_mat_real
    Type( complex_distributed_matrix ), Allocatable :: A_mat_complex, B_mat_complex

    Logical :: is_real
    Logical :: p_A, p_B

    p_A = Allocated( A )
    p_B = Allocated( B )
    
    If( .Not. p_A .And. .Not. p_B ) Then
       Stop "In ks_matrix_remap_data_complex one of A or B must be supplied"
    End If

    If( p_A .And. p_B ) Then
       If( .Not. same_type_as( A%k_point%matrix, B%k_point%matrix ) ) Then
          Stop "Type mismatch in ks_matrix_remap_data"
       End If
    End If

    If( p_A ) Then
       Associate( Akm => A%k_point%matrix )
         Select Type( Akm )
         Class Default
            Stop "Illegal type in ks_matrix_remap_data"
         Type is ( real_distributed_matrix )
            Allocate( A_mat_real, Source = Akm ) 
         Type is ( complex_distributed_matrix )
            Allocate( A_mat_complex, Source = Akm )
         End Select
       End Associate
    End If

    If( p_B ) Then
       Associate( Bkm => B%k_point%matrix )
         Select Type( Bkm )
         Class Default
            Stop "Illegal type in ks_matrix_remap_data"
         Type is ( real_distributed_matrix )
            Allocate( B_mat_real, Source = Bkm ) 
         Type is ( complex_distributed_matrix )
            Allocate( B_mat_complex, Source = Bkm )
         End Select
       End Associate
    End If

    ! Now we know
    ! 1) If there is both A and B they are of the same type
    ! 2) At least one of them is allocated
    ! Thus can work out easily what typ we are dealing with
    is_real = Allocated( A_mat_real ) .Or. Allocated( B_mat_real )

    If( is_real ) Then
       Call distributed_matrix_remap_data( parent_comm, A_mat_real, B_mat_real )
       If( p_B ) Then
          Deallocate( B%k_point%matrix )
          Allocate  ( B%k_point%matrix, Source = B_mat_real )
       End If
    Else
       Call distributed_matrix_remap_data( parent_comm, A_mat_complex, B_mat_complex )
       If( p_B ) Then
          Deallocate( B%k_point%matrix )
          Allocate  ( B%k_point%matrix, Source = B_mat_complex )
       End If
    End If
       
  End Subroutine ks_matrix_remap_data

  Subroutine ks_matrix_set_to_identity( A )
    
    Class( ks_matrix ), Intent( InOut ) :: A

    Associate( Akm => A%k_point%matrix )
      Select Type( Akm )
      Class Default
         Stop "Illegal type in ks_matrix_set_to_identity"
      Type is ( real_distributed_matrix )
         Call Akm%set_to_identity()
      Type is ( complex_distributed_matrix )
         Call Akm%set_to_identity()
      End Select
    End Associate

  End Subroutine ks_matrix_set_to_identity

  Function ks_matrix_g_to_l( A, what ) Result( gl_indexing )

    Integer, Dimension( : ), Allocatable :: gl_indexing

    Class( ks_matrix ), Intent( In ) :: A
    Character( Len = * )         , Intent( In ) :: what

    gl_indexing = A%k_point%matrix%global_to_local( what )

  End Function ks_matrix_g_to_l

  Function ks_matrix_l_to_g( A, what ) Result( lg_indexing )

    Integer, Dimension( : ), Allocatable :: lg_indexing

    Class( ks_matrix ), Intent( In ) :: A
    Character( Len = * )         , Intent( In ) :: what

    lg_indexing = A%k_point%matrix%local_to_global( what )

  End Function ks_matrix_l_to_g

  Function ks_matrix_size( A, dim ) Result( n )

    Integer :: n

    Class( ks_matrix ), Intent( In )           :: A
    Integer                      , Intent( In ), Optional :: dim

    If( .Not. Present( dim ) ) Then
       n = A%k_point%matrix%size( 1 ) * A%k_point%matrix%size( 2 )
    Else
       n = A%k_point%matrix%size( dim )
    End If

  End Function ks_matrix_size

  Function ks_matrix_local_size( A, dim ) Result( n )

    Integer :: n

    Class( ks_matrix ), Intent( In )           :: A
    Integer                      , Intent( In ), Optional :: dim

    If( .Not. Present( dim ) ) Then
       n = A%k_point%matrix%local_size( 1 ) * A%k_point%matrix%local_size( 2 )
    Else
       n = A%k_point%matrix%local_size( dim )
    End If

  End Function ks_matrix_local_size

  Function ks_matrix_get_communicator( A ) Result( c )

    Integer :: c

    Class( ks_matrix ), Intent( In ) :: A

    c = A%k_point%matrix%get_comm()
    
  End Function ks_matrix_get_communicator

End Module ks_matrix_module

