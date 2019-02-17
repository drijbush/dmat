Module distributed_k_module

  ! Wrapper for dist matrix types so can make opaque array - only 1 k point at this level

  Use numbers_module       , Only : wp
  Use distributed_matrix_module, Only : distributed_matrix, real_distributed_matrix, complex_distributed_matrix, &
       distributed_matrix_init, distributed_matrix_finalise
  Implicit None

  Type, Private :: k_point_matrix
     Integer                                 , Private :: this_spin
     Integer, Dimension( 1:3 )               , Private :: this_k_point
     Class( distributed_matrix ), Allocatable, Private :: matrix
  End Type k_point_matrix

  Type, Extends( k_point_matrix ), Private :: k_wave_function
     Real( wp ), Dimension( : ), Allocatable, Private :: evals
  End Type k_wave_function

  Type, Public :: distributed_k_matrix
     Class( k_point_matrix ), Allocatable, Private :: k_point
   Contains
     Procedure            :: create               => distributed_k_matrix_create
     Procedure            :: dagger               => distributed_k_matrix_dagger
     Generic              :: operator( .Dagger. ) => dagger
     Procedure            :: diag                 => distributed_k_matrix_diag
     Procedure            :: multiply             => distributed_k_matrix_mult
     Generic              :: Operator( * )        => multiply
     Procedure            :: post_scale           => distributed_k_matrix_post_scale
     Procedure, Pass( A ) :: pre_scale            => distributed_k_matrix_pre_scale
     Generic              :: Operator( * )        => post_scale, pre_scale
     Procedure            :: add                  => distributed_k_matrix_add
     Generic              :: Operator( + )        => add
     Procedure            :: subtract             => distributed_k_matrix_subtract
     Generic              :: Operator( - )        => subtract
     Procedure            :: Choleski             => distributed_k_matrix_Choleski
     Procedure            :: Solve                => distributed_k_matrix_Solve
     Procedure            :: set_to_identity      => distributed_k_matrix_set_to_identity
     Procedure, Private   :: sgr                  => set_global_real
     Procedure, Private   :: sgc                  => set_global_complex
     Generic              :: set_by_global        => sgr, sgc
     Procedure, Private   :: slr                  => set_local_real
     Procedure, Private   :: slc                  => set_local_complex
     Generic              :: set_by_local         => slr, slc
     Procedure, Private   :: ggr                  => get_global_real
     Procedure, Private   :: ggc                  => get_global_complex
     Generic              :: get_by_global        => ggr, ggc
     Procedure, Private   :: glr                  => get_local_real
     Procedure, Private   :: glc                  => get_local_complex
     Generic              :: get_by_local         => glr, glc
     Procedure            :: extract_cols         => distributed_k_matrix_extract_cols
     Procedure            :: extract              => distributed_k_matrix_extract
     Procedure            :: global_to_local      => distributed_k_matrix_g_to_l
     Procedure            :: local_to_global      => distributed_k_matrix_l_to_g
     Procedure            :: local_size           => distributed_k_matrix_local_size
  End Type distributed_k_matrix
  
  ! a) Set up base matrix which is 1 matrix across all
  ! b) Provide splits

  Public :: distributed_k_matrix_init
  Public :: distributed_k_matrix_finalise
  
  Private

  Integer, Parameter :: INVALID = -1

Contains

  Subroutine distributed_k_matrix_init( comm, base_matrix )

    Integer                        , Intent( In    ) :: comm
    Type   ( distributed_k_matrix ), Intent(   Out ) :: base_matrix

    Allocate( k_point_matrix:: base_matrix%k_point )

    base_matrix%k_point%this_spin    = INVALID
    base_matrix%k_point%this_k_point = INVALID

    Allocate( real_distributed_matrix:: base_matrix%k_point%matrix )

    Call distributed_matrix_init( comm, base_matrix%k_point%matrix )
    
  End Subroutine distributed_k_matrix_init

  Subroutine distributed_k_matrix_finalise

    Call distributed_matrix_finalise
    
  End Subroutine distributed_k_matrix_finalise
  
  Subroutine distributed_k_matrix_create( matrix, is_complex, this_spin, this_k_point, m, n, source_matrix )

    Class  ( distributed_k_matrix ), Intent(   Out ) :: matrix
    Logical                        , Intent( In    ) :: is_complex
    Integer                        , Intent( In    ) :: this_spin
    Integer, Dimension( 1:3 )      , Intent( In    ) :: this_k_point
    Integer                        , Intent( In    ) :: m
    Integer                        , Intent( In    ) :: n
    Type   ( distributed_k_matrix ), Intent( In    ) :: source_matrix

    Allocate( k_point_matrix:: matrix%k_point )

    If( is_complex ) Then
      Allocate( complex_distributed_matrix :: matrix%k_point%matrix ) 
    Else
      Allocate( real_distributed_matrix    :: matrix%k_point%matrix ) 
    End If

    matrix%k_point%this_spin    = this_spin
    matrix%k_point%this_k_point = this_k_point

    Call matrix%k_point%matrix%create( m, n, source_matrix%k_point%matrix )
    
  End Subroutine distributed_k_matrix_create

  Subroutine distributed_k_matrix_diag( A, Q, evals )
    
    Class( distributed_k_matrix )          , Intent( In    ) :: A
    Type ( distributed_k_matrix )          , Intent(   Out ) :: Q
    Real( wp ), Dimension( : ), Allocatable, Intent(   Out ) :: evals

    Type(    real_distributed_matrix ), Allocatable :: Q_real
    Type( complex_distributed_matrix ), Allocatable :: Q_complex
    
    Allocate( k_wave_function  :: Q%k_point )

    Q%k_point%this_spin    = A%k_point%this_spin
    Q%k_point%this_k_point = A%k_point%this_k_point

    Associate( Akm => A%k_point%matrix )
    
      Select Type( Akm )

      Class Default
         Stop "Illegal type in distributed_k_matrix_diag"

      Class is ( real_distributed_matrix )
         Call Akm%diag( Q_real, evals )
         Allocate( Q%k_point%matrix, Source = Q_real )

      Class is ( complex_distributed_matrix )
         Call Akm%diag( Q_complex, evals )
         Allocate( Q%k_point%matrix, Source = Q_complex )

      End Select

    End Associate

    Associate( Qk => Q%k_point )
      Select Type( Qk )
      Class is ( k_wave_function )
         Qk%evals = evals
      End Select
    End Associate
      
  End Subroutine distributed_k_matrix_diag
  
  Function distributed_k_matrix_dagger( A ) Result( tA )
    
    Type( distributed_k_matrix ), Allocatable :: tA

    Class( distributed_k_matrix ), Intent( In ) :: A

    Type(    real_distributed_matrix ) :: tA_real
    Type( complex_distributed_matrix ) :: tA_complex

    Allocate( tA )
    Associate( Ak => A%k_point )
      Select Type( Ak )
      Class Default
         Stop "Illegal type in distributed_k_matrix_diag"
      Type is ( k_point_matrix )
         Allocate( k_point_matrix :: tA%k_point )
      Type is ( k_wave_function )
         Allocate( k_wave_function :: tA%k_point )
      End Select
    End Associate
    
    tA%k_point%this_spin    = A%k_point%this_spin
    tA%k_point%this_k_point = A%k_point%this_k_point

    Associate( Akm => A%k_point%matrix )
      Select Type( Akm )
      Class Default
         Stop "Illegal type in distributed_k_matrix_diag"
      Type is ( real_distributed_matrix )
         tA_real = .Dagger. Akm
         Allocate( tA%k_point%matrix, Source = tA_real )
      Type is ( complex_distributed_matrix )
         tA_complex = .Dagger. Akm
         Allocate( tA%k_point%matrix, Source = tA_complex )
      End Select
         
    End Associate

  End Function distributed_k_matrix_dagger
    
  Function distributed_k_matrix_Choleski( A ) Result( C )
    
    Type( distributed_k_matrix ), Allocatable :: C

    Class( distributed_k_matrix ), Intent( In ) :: A

    Type(    real_distributed_matrix ) :: C_real
    Type( complex_distributed_matrix ) :: C_complex

    Allocate( C )
    Associate( Ak => A%k_point )
      Select Type( Ak )
      Class Default
         Stop "Illegal type in distributed_k_matrix_diag"
      Type is ( k_point_matrix )
         Allocate( k_point_matrix :: C%k_point )
      Type is ( k_wave_function )
         Allocate( k_wave_function :: C%k_point )
      End Select
    End Associate
    
    C%k_point%this_spin    = A%k_point%this_spin
    C%k_point%this_k_point = A%k_point%this_k_point

    Associate( Akm => A%k_point%matrix )
      Select Type( Akm )
      Class Default
         Stop "Illegal type in distributed_k_matrix_diag"
      Type is ( real_distributed_matrix )
         C_real = Akm%Choleski()
         Allocate( C%k_point%matrix, Source = C_real )
      Type is ( complex_distributed_matrix )
         C_complex = Akm%Choleski()
         Allocate( C%k_point%matrix, Source = C_complex )
      End Select
         
    End Associate

  End Function distributed_k_matrix_Choleski
    
  Function distributed_k_matrix_mult( A, B ) Result( C )
    
    Type( distributed_k_matrix ), Allocatable :: C

    Class( distributed_k_matrix ), Intent( In ) :: A
    Type ( distributed_k_matrix ), Intent( In ) :: B

    Type(    real_distributed_matrix ) :: C_real
    Type( complex_distributed_matrix ) :: C_complex

    Allocate( C )
    Associate( Ak => A%k_point )
      Select Type( Ak )
      Class Default
         Stop "Illegal type in distributed_k_matrix_diag"
      Type is ( k_point_matrix )
         Allocate( k_point_matrix :: C%k_point )
      Type is ( k_wave_function )
         Allocate( k_wave_function :: C%k_point )
      End Select
    End Associate
    C%k_point%this_spin    = A%k_point%this_spin
    C%k_point%this_k_point = A%k_point%this_k_point

    Associate( Akm => A%k_point%matrix )
    
      Select Type( Akm )

      Class Default
         Stop "Illegal type in distributed_k_matrix_diag"
         
      Type is ( real_distributed_matrix )
         
         Associate( Bkm => B%k_point%matrix )
           Select Type( Bkm )
           Class Default
              Stop "Illegal type in distributed_k_matrix_diag"
           Type is ( real_distributed_matrix )
              C_real = Akm * Bkm
              Allocate( C%k_point%matrix, Source = C_real )
           End Select
         End Associate

      Type is ( complex_distributed_matrix )
         Associate( Bkm => B%k_point%matrix )
           Select Type( Bkm )
           Class Default
              Stop "Illegal type in distributed_k_matrix_diag"
           Type is ( complex_distributed_matrix )
              C_complex = Akm * Bkm
              Allocate( C%k_point%matrix, Source = C_complex )
           End Select
         End Associate
         
      End Select
    End Associate

  End Function distributed_k_matrix_mult

  Function distributed_k_matrix_solve( A, B ) Result( C )
    
    Type( distributed_k_matrix ), Allocatable :: C

    Class( distributed_k_matrix ), Intent( In ) :: A
    Type ( distributed_k_matrix ), Intent( In ) :: B

    Type(    real_distributed_matrix ) :: C_real
    Type( complex_distributed_matrix ) :: C_complex

    Allocate( C )
    Associate( Ak => A%k_point )
      Select Type( Ak )
      Class Default
         Stop "Illegal type in distributed_k_matrix_diag"
      Type is ( k_point_matrix )
         Allocate( k_point_matrix :: C%k_point )
      Type is ( k_wave_function )
         Allocate( k_wave_function :: C%k_point )
      End Select
    End Associate
    C%k_point%this_spin    = A%k_point%this_spin
    C%k_point%this_k_point = A%k_point%this_k_point

    Associate( Akm => A%k_point%matrix )
    
      Select Type( Akm )

      Class Default
         Stop "Illegal type in distributed_k_matrix_diag"
         
      Type is ( real_distributed_matrix )
         
         Associate( Bkm => B%k_point%matrix )
           Select Type( Bkm )
           Class Default
              Stop "Illegal type in distributed_k_matrix_diag"
           Type is ( real_distributed_matrix )
              C_real = Akm%solve( Bkm )
              Allocate( C%k_point%matrix, Source = C_real )
           End Select
         End Associate

      Type is ( complex_distributed_matrix )
         Associate( Bkm => B%k_point%matrix )
           Select Type( Bkm )
           Class Default
              Stop "Illegal type in distributed_k_matrix_diag"
           Type is ( complex_distributed_matrix )
              C_complex = Akm%solve( Bkm )
              Allocate( C%k_point%matrix, Source = C_complex )
           End Select
         End Associate
         
      End Select
    End Associate

  End Function distributed_k_matrix_solve

  Function distributed_k_matrix_add( A, B ) Result( C )
    
    Type( distributed_k_matrix ), Allocatable :: C

    Class( distributed_k_matrix ), Intent( In ) :: A
    Type ( distributed_k_matrix ), Intent( In ) :: B

    Type(    real_distributed_matrix ) :: C_real
    Type( complex_distributed_matrix ) :: C_complex

    Allocate( C )
    Associate( Ak => A%k_point )
      Select Type( Ak )
      Class Default
         Stop "Illegal type in distributed_k_matrix_diag"
      Type is ( k_point_matrix )
         Allocate( k_point_matrix :: C%k_point )
      Type is ( k_wave_function )
         Allocate( k_wave_function :: C%k_point )
      End Select
    End Associate
    C%k_point%this_spin    = A%k_point%this_spin
    C%k_point%this_k_point = A%k_point%this_k_point

    Associate( Akm => A%k_point%matrix )
    
      Select Type( Akm )

      Class Default
         Stop "Illegal type in distributed_k_matrix_diag"
         
      Type is ( real_distributed_matrix )
         
         Associate( Bkm => B%k_point%matrix )
           Select Type( Bkm )
           Class Default
              Stop "Illegal type in distributed_k_matrix_diag"
           Type is ( real_distributed_matrix )
              C_real = Akm + Bkm
              Allocate( C%k_point%matrix, Source = C_real )
           End Select
         End Associate

      Type is ( complex_distributed_matrix )
         Associate( Bkm => B%k_point%matrix )
           Select Type( Bkm )
           Class Default
              Stop "Illegal type in distributed_k_matrix_diag"
           Type is ( complex_distributed_matrix )
              C_complex = Akm + Bkm
              Allocate( C%k_point%matrix, Source = C_complex )
           End Select
         End Associate
         
      End Select
    End Associate

  End Function distributed_k_matrix_add

  Function distributed_k_matrix_subtract( A, B ) Result( C )
    
    Type( distributed_k_matrix ), Allocatable :: C

    Class( distributed_k_matrix ), Intent( In ) :: A
    Type ( distributed_k_matrix ), Intent( In ) :: B

    Type(    real_distributed_matrix ) :: C_real
    Type( complex_distributed_matrix ) :: C_complex

    Allocate( C )
    Associate( Ak => A%k_point )
      Select Type( Ak )
      Class Default
         Stop "Illegal type in distributed_k_matrix_diag"
      Type is ( k_point_matrix )
         Allocate( k_point_matrix :: C%k_point )
      Type is ( k_wave_function )
         Allocate( k_wave_function :: C%k_point )
      End Select
    End Associate
    C%k_point%this_spin    = A%k_point%this_spin
    C%k_point%this_k_point = A%k_point%this_k_point

    Associate( Akm => A%k_point%matrix )
    
      Select Type( Akm )

      Class Default
         Stop "Illegal type in distributed_k_matrix_diag"
         
      Type is ( real_distributed_matrix )
         
         Associate( Bkm => B%k_point%matrix )
           Select Type( Bkm )
           Class Default
              Stop "Illegal type in distributed_k_matrix_diag"
           Type is ( real_distributed_matrix )
              C_real = Akm - Bkm
              Allocate( C%k_point%matrix, Source = C_real )
           End Select
         End Associate

      Type is ( complex_distributed_matrix )
         Associate( Bkm => B%k_point%matrix )
           Select Type( Bkm )
           Class Default
              Stop "Illegal type in distributed_k_matrix_diag"
           Type is ( complex_distributed_matrix )
              C_complex = Akm - Bkm
              Allocate( C%k_point%matrix, Source = C_complex )
           End Select
         End Associate
         
      End Select
    End Associate

  End Function distributed_k_matrix_subtract

  Subroutine set_global_real( A, m, n, p, q, data )
    
    Class( distributed_k_matrix )   , Intent( InOut ) :: A
    Integer                         , Intent( In    ) :: m
    Integer                         , Intent( In    ) :: n
    Integer                         , Intent( In    ) :: p
    Integer                         , Intent( In    ) :: q
    Real( wp ), Dimension( m:, p: ) , Intent( In    ) :: data

    Associate( Akm => A%k_point%matrix )
    
      Select Type( Akm )

      Class Default
         Stop "Illegal type in distributed_k_matrix_diag"
         
      Type is ( real_distributed_matrix )
         Call Akm%set_by_global( m, n, p, q, data )

      End Select
      
    End Associate
    
  End Subroutine set_global_real
    
  Subroutine set_global_complex( A, m, n, p, q, data )
    
    Class( distributed_k_matrix )      , Intent( InOut ) :: A
    Integer                            , Intent( In    ) :: m
    Integer                            , Intent( In    ) :: n
    Integer                            , Intent( In    ) :: p
    Integer                            , Intent( In    ) :: q
    Complex( wp ), Dimension( m:, p: ) , Intent( In    ) :: data

    Associate( Akm => A%k_point%matrix )
    
      Select Type( Akm )

      Class Default
         Stop "Illegal type in distributed_k_matrix_diag"
         
      Type is ( complex_distributed_matrix )
         Call Akm%set_by_global( m, n, p, q, data )

      End Select
      
    End Associate
    
  End Subroutine set_global_complex
    
  Subroutine set_local_real( A, m, n, p, q, data )
    
    Class( distributed_k_matrix )   , Intent( InOut ) :: A
    Integer                         , Intent( In    ) :: m
    Integer                         , Intent( In    ) :: n
    Integer                         , Intent( In    ) :: p
    Integer                         , Intent( In    ) :: q
    Real( wp ), Dimension( m:, p: ) , Intent( In    ) :: data

    Associate( Akm => A%k_point%matrix )
    
      Select Type( Akm )

      Class Default
         Stop "Illegal type in distributed_k_matrix_diag"
         
      Type is ( real_distributed_matrix )
         Call Akm%set_by_local( m, n, p, q, data )

      End Select
      
    End Associate
    
  End Subroutine set_local_real
    
  Subroutine set_local_complex( A, m, n, p, q, data )
    
    Class( distributed_k_matrix )      , Intent( InOut ) :: A
    Integer                            , Intent( In    ) :: m
    Integer                            , Intent( In    ) :: n
    Integer                            , Intent( In    ) :: p
    Integer                            , Intent( In    ) :: q
    Complex( wp ), Dimension( m:, p: ) , Intent( In    ) :: data

    Associate( Akm => A%k_point%matrix )
    
      Select Type( Akm )

      Class Default
         Stop "Illegal type in distributed_k_matrix_diag"
         
      Type is ( complex_distributed_matrix )
         Call Akm%set_by_local( m, n, p, q, data )

      End Select
      
    End Associate
    
  End Subroutine set_local_complex
    
  Subroutine get_global_real( A, m, n, p, q, data )
    
    Class( distributed_k_matrix )   , Intent( In    ) :: A
    Integer                         , Intent( In    ) :: m
    Integer                         , Intent( In    ) :: n
    Integer                         , Intent( In    ) :: p
    Integer                         , Intent( In    ) :: q
    Real( wp ), Dimension( m:, p: ) , Intent(   Out ) :: data

    Associate( Akm => A%k_point%matrix )
    
      Select Type( Akm )

      Class Default
         Stop "Illegal type in distributed_k_matrix_diag"
         
      Type is ( real_distributed_matrix )
         Call Akm%get_by_global( m, n, p, q, data )

      End Select
      
    End Associate
    
  End Subroutine get_global_real
    
  Subroutine get_global_complex( A, m, n, p, q, data )
    
    Class( distributed_k_matrix )      , Intent( In    ) :: A
    Integer                            , Intent( In    ) :: m
    Integer                            , Intent( In    ) :: n
    Integer                            , Intent( In    ) :: p
    Integer                            , Intent( In    ) :: q
    Complex( wp ), Dimension( m:, p: ) , Intent(   Out ) :: data

    Associate( Akm => A%k_point%matrix )
    
      Select Type( Akm )

      Class Default
         Stop "Illegal type in distributed_k_matrix_diag"
         
      Type is ( complex_distributed_matrix )
         Call Akm%get_by_global( m, n, p, q, data )

      End Select
      
    End Associate
    
  End Subroutine get_global_complex
    
  Subroutine get_local_real( A, m, n, p, q, data )
    
    Class( distributed_k_matrix )   , Intent( In    ) :: A
    Integer                         , Intent( In    ) :: m
    Integer                         , Intent( In    ) :: n
    Integer                         , Intent( In    ) :: p
    Integer                         , Intent( In    ) :: q
    Real( wp ), Dimension( m:, p: ) , Intent(   Out ) :: data

    Associate( Akm => A%k_point%matrix )
    
      Select Type( Akm )

      Class Default
         Stop "Illegal type in distributed_k_matrix_diag"
         
      Type is ( real_distributed_matrix )
         Call Akm%get_by_local( m, n, p, q, data )

      End Select
      
    End Associate
    
  End Subroutine get_local_real
    
  Subroutine get_local_complex( A, m, n, p, q, data )
    
    Class( distributed_k_matrix )      , Intent( In    ) :: A
    Integer                            , Intent( In    ) :: m
    Integer                            , Intent( In    ) :: n
    Integer                            , Intent( In    ) :: p
    Integer                            , Intent( In    ) :: q
    Complex( wp ), Dimension( m:, p: ) , Intent(   Out ) :: data

    Associate( Akm => A%k_point%matrix )
    
      Select Type( Akm )

      Class Default
         Stop "Illegal type in distributed_k_matrix_diag"
         
      Type is ( complex_distributed_matrix )
         Call Akm%get_by_local( m, n, p, q, data )

      End Select
      
    End Associate
    
  End Subroutine get_local_complex

  Subroutine distributed_k_matrix_extract_cols( A, c1, c2, B )
    
    Class( distributed_k_matrix ), Intent( In    ) :: A
    Integer                      , Intent( In    ) :: c1 
    Integer                      , Intent( In    ) :: c2
    Type ( distributed_k_matrix ), Intent(   Out ) :: B

    Type(    real_distributed_matrix ) :: B_real
    Type( complex_distributed_matrix ) :: B_complex
    
    Associate( Ak => A%k_point )
      Select Type( Ak )
      Class Default
         Stop "Illegal type in distributed_k_matrix_extract_cols"
      Type is ( k_point_matrix )
         Allocate( k_point_matrix :: B%k_point )
      Type is ( k_wave_function )
         Allocate( k_wave_function :: B%k_point )
         Associate( Bk => B%k_point )
           Select Type( Bk )
           Type is ( k_wave_function )
              Bk%evals = Ak%evals( c1:c2 )
           End Select
         End Associate
      End Select
    End Associate
    B%k_point%this_spin    = A%k_point%this_spin
    B%k_point%this_k_point = A%k_point%this_k_point

    Associate( Akm => A%k_point%matrix )
    
      Select Type( Akm )

      Class Default
         Stop "Illegal type in distributed_k_matrix_extract_cols"
         
      Type is ( real_distributed_matrix )
         Call Akm%extract_cols( c1, c2, B_real )
         Allocate( B%k_point%matrix, Source = B_real )

      Type is ( complex_distributed_matrix )
         Call Akm%extract_cols( c1, c2, B_complex )
         Allocate( B%k_point%matrix, Source = B_complex )

      End Select
    End Associate
         
  End Subroutine distributed_k_matrix_extract_cols

  Function distributed_k_matrix_extract( A, r1, r2, c1, c2 ) Result( B )
    
    Type ( distributed_k_matrix ), Allocatable :: B

    Class( distributed_k_matrix ), Intent( In    ) :: A
    Integer                      , Intent( In    ) :: r1 
    Integer                      , Intent( In    ) :: r2
    Integer                      , Intent( In    ) :: c1 
    Integer                      , Intent( In    ) :: c2

    Type(    real_distributed_matrix ) :: B_real
    Type( complex_distributed_matrix ) :: B_complex

    Allocate( B )
    Associate( Ak => A%k_point )
      Select Type( Ak )
      Class Default
         Stop "Illegal type in distributed_k_matrix_extract"
      Type is ( k_point_matrix )
         Allocate( k_point_matrix :: B%k_point )
      Type is ( k_wave_function )
         Allocate( k_wave_function :: B%k_point )
         Associate( Bk => B%k_point )
           Select Type( Bk )
           Type is ( k_wave_function )
              Bk%evals = Ak%evals( c1:c2 )
           End Select
         End Associate
      End Select
    End Associate
    B%k_point%this_spin    = A%k_point%this_spin
    B%k_point%this_k_point = A%k_point%this_k_point

    Associate( Akm => A%k_point%matrix )
    
      Select Type( Akm )

      Class Default
         Stop "Illegal type in distributed_k_matrix_extract"
         
      Type is ( real_distributed_matrix )
         B_real = Akm%extract( r1, r2, c1, c2 )
         Allocate( B%k_point%matrix, Source = B_real )

      Type is ( complex_distributed_matrix )
         B_complex = Akm%extract( r1, r2, c1, c2 )
         Allocate( B%k_point%matrix, Source = B_complex )

      End Select
    End Associate
         
  End Function distributed_k_matrix_extract
  
  Function distributed_k_matrix_post_scale( A, s ) Result( B )
    
    Type( distributed_k_matrix ), Allocatable :: B

    Class( distributed_k_matrix ), Intent( In ) :: A
    Real ( wp )                  , Intent( In ) :: s

    Type(    real_distributed_matrix ) :: B_real
    Type( complex_distributed_matrix ) :: B_complex

    Allocate( B )
    Associate( Ak => A%k_point )
      Select Type( Ak )
      Class Default
         Stop "Illegal type in distributed_k_matrix_post_scale"
      Type is ( k_point_matrix )
         Allocate( k_point_matrix :: B%k_point )
      Type is ( k_wave_function )
         Allocate( k_wave_function :: B%k_point )
      End Select
    End Associate
    
    B%k_point%this_spin    = A%k_point%this_spin
    B%k_point%this_k_point = A%k_point%this_k_point

    Associate( Akm => A%k_point%matrix )
      Select Type( Akm )
      Class Default
         Stop "Illegal type in distributed_k_matrix_post_scale"
      Type is ( real_distributed_matrix )
         B_real = s * Akm
         Allocate( B%k_point%matrix, Source = B_real )
      Type is ( complex_distributed_matrix )
         B_complex = Cmplx( s, Kind = wp ) * Akm
         Allocate( B%k_point%matrix, Source = B_complex )
      End Select
    End Associate

  End Function distributed_k_matrix_post_scale

  Function distributed_k_matrix_pre_scale( s, A ) Result( B )
    
    Type( distributed_k_matrix ), Allocatable :: B

    Real ( wp )                  , Intent( In ) :: s
    Class( distributed_k_matrix ), Intent( In ) :: A

    Type(    real_distributed_matrix ) :: B_real
    Type( complex_distributed_matrix ) :: B_complex

    Allocate( B )
    Associate( Ak => A%k_point )
      Select Type( Ak )
      Class Default
         Stop "Illegal type in distributed_k_matrix_pre_scale"
      Type is ( k_point_matrix )
         Allocate( k_point_matrix :: B%k_point )
      Type is ( k_wave_function )
         Allocate( k_wave_function :: B%k_point )
      End Select
    End Associate
    
    B%k_point%this_spin    = A%k_point%this_spin
    B%k_point%this_k_point = A%k_point%this_k_point

    Associate( Akm => A%k_point%matrix )
      Select Type( Akm )
      Class Default
         Stop "Illegal type in distributed_k_matrix_pre_scale"
      Type is ( real_distributed_matrix )
         B_real = s * Akm
         Allocate( B%k_point%matrix, Source = B_real )
      Type is ( complex_distributed_matrix )
         B_complex = Cmplx( s, Kind = wp ) * Akm
         Allocate( B%k_point%matrix, Source = B_complex )
      End Select
    End Associate

  End Function distributed_k_matrix_pre_scale

  Subroutine distributed_k_matrix_set_to_identity( A )
    
    Class( distributed_k_matrix ), Intent( InOut ) :: A

    Associate( Akm => A%k_point%matrix )
      Select Type( Akm )
      Class Default
         Stop "Illegal type in distributed_k_matrix_pre_scale"
      Type is ( real_distributed_matrix )
         Call Akm%set_to_identity()
      Type is ( complex_distributed_matrix )
         Call Akm%set_to_identity()
      End Select
    End Associate

  End Subroutine distributed_k_matrix_set_to_identity

  Function distributed_k_matrix_g_to_l( A, what ) Result( gl_indexing )

    Integer, Dimension( : ), Allocatable :: gl_indexing

    Class( distributed_k_matrix ), Intent( In ) :: A
    Character( Len = * )         , Intent( In ) :: what

    gl_indexing = A%k_point%matrix%global_to_local( what )

  End Function distributed_k_matrix_g_to_l

  Function distributed_k_matrix_l_to_g( A, what ) Result( lg_indexing )

    Integer, Dimension( : ), Allocatable :: lg_indexing

    Class( distributed_k_matrix ), Intent( In ) :: A
    Character( Len = * )         , Intent( In ) :: what

    lg_indexing = A%k_point%matrix%local_to_global( what )

  End Function distributed_k_matrix_l_to_g

  Function distributed_k_matrix_local_size( A, dim ) Result( n )

    Integer :: n

    Class( distributed_k_matrix ), Intent( In )           :: A
    Integer                      , Intent( In ), Optional :: dim

    If( .Not. Present( dim ) ) Then
       n = A%k_point%matrix%local_size( 1 ) * A%k_point%matrix%local_size( 2 )
    Else
       n = A%k_point%matrix%local_size( dim )
    End If

  End Function distributed_k_matrix_local_size

End Module distributed_k_module
