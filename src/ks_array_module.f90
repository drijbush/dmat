Module ks_array_module

  ! IRREPS NEED MORE THOUGHT
  !-------------------------
  ! Currently assume 1 - i.e. nosymada
  
  Use numbers_module  , Only : wp
  Use ks_matrix_module, Only : ks_matrix_init, ks_matrix_comm_to_base, &
       ks_matrix_finalise, ks_matrix_remap_data, ks_matrix

  Implicit None

  Integer, Public, Parameter :: K_POINT_REAL = 0
  Integer, Public, Parameter :: K_POINT_COMPLEX = 1

  Integer, Private, Parameter :: INVALID = -1
  Integer, Private, Parameter :: NOT_ME  = -2

  Type k_point_info
     Integer                              :: k_type
     Integer                              :: spin
     Integer, Dimension( : ), Allocatable :: k_indices
  End type k_point_info

  Type, Private :: k_point_matrices
     ! External label - e.g. irrep number
     Integer                      :: label = INVALID
     Type( ks_matrix ) :: matrix
  End type k_point_matrices

  Type, Private :: k_point
     Type( k_point_info     )                              :: info
     ! This splitting allows irreps within this k point
     Type( k_point_matrices ), Dimension( : ), Allocatable :: data
     ! Want to hide eventually
     Integer                                               :: communicator = INVALID
  End type k_point

  Type, Public :: ks_array
     Type( k_point_info ), Dimension( : ), Allocatable :: all_k_point_info
     ! this splitting allows multiple k points on this process
     Type( k_point      ), Dimension( : ), Allocatable :: my_k_points
     ! Want to hide eventually
     Integer                                           :: parent_communicator = INVALID
   Contains
     Procedure                     :: create               => ks_array_create
     Procedure                     :: split_ks             => ks_array_split_ks
     Procedure                     :: diag                 => ks_array_diag
     Procedure                     :: print_info           => ks_array_print_info
     Procedure, Private            :: get_all_ks_index
     Procedure, Private            :: get_my_ks_index
     Procedure, Private            :: get_ks
     Procedure, Private            :: set_by_global_r      => ks_array_set_global_real
     Procedure, Private            :: set_by_global_c      => ks_array_set_global_complex
     Generic                       :: set_by_global        => set_by_global_r, set_by_global_c
     Procedure, Private            :: get_by_global_r      => ks_array_get_global_real
     Procedure, Private            :: get_by_global_c      => ks_array_get_global_complex
     Generic                       :: get_by_global        => get_by_global_r, get_by_global_c
     Procedure, Private            :: multiply             => ks_array_mult
     Procedure, Private, Pass( A ) :: pre_scale            => ks_array_pre_scale
     Procedure, Private            :: post_scale           => ks_array_post_scale
     Procedure, Private, Pass( A ) :: pre_mult_diag        => ks_array_pre_mult_diag
     Procedure, Private            :: post_mult_diag       => ks_array_post_mult_diag
     Generic                       :: Operator( * )        => multiply, pre_scale, post_scale, &
                                                              pre_mult_diag, post_mult_diag
     Procedure, Private            :: dagger               => ks_array_dagger
     Generic                       :: Operator( .Dagger. ) => dagger
     Procedure, Private            :: add                  => ks_array_add
     Generic                       :: Operator( + )        => add
  End type ks_array
  
  Type, Public :: eval_storage
     Integer                                 :: spin
     Integer   , Dimension( : ), Allocatable :: k_indices
     Real( wp ), Dimension( : ), Allocatable :: evals
  End type eval_storage

  Public :: ks_array_init
  Public :: ks_array_comm_to_base
  Public :: ks_array_finalise
  
  Private

Contains

  Subroutine ks_array_init
    Call ks_matrix_init
  End Subroutine ks_array_init

  Subroutine ks_array_comm_to_base( comm, n_spin, k_point_type, k_points, base_ks_array )

    ! Want to think about what kind of thing is returned here ...
    
    Integer                   , Intent( In    ) :: comm
    Integer                   , Intent( In    ) :: n_spin
    Integer, Dimension( :    ), Intent( In    ) :: k_point_type
    Integer, Dimension( :, : ), Intent( In    ) :: k_points
    Type( ks_array )          , Intent(   Out ) :: base_ks_array

    Type( ks_matrix ) :: base_matrix
    
    Integer :: n_k_points
    Integer :: s, k, ks
    
    n_k_points = Size( k_point_type)

    ! Set up the all k point data structure
    Allocate( base_ks_array%all_k_point_info( 1:n_spin * n_k_points ) )
    ks = 0
    Do s = 1, n_spin
       Do k = 1, n_k_points
          ks = ks + 1
          base_ks_array%all_k_point_info( ks )%k_type    = k_point_type( k )
          base_ks_array%all_k_point_info( ks )%spin      = s
          base_ks_array%all_k_point_info( ks )%k_indices = k_points( :, k )
       End Do
    End Do

    Call ks_matrix_comm_to_base( comm, base_matrix ) 

    ! Now my k points
    Allocate( base_ks_array%my_k_points( 1:n_spin * n_k_points ) )
    ks = 0
    Do s = 1, n_spin
       Do k = 1, n_k_points
          ks = ks + 1
          base_ks_array%my_k_points( ks )%info%k_type    = k_point_type( k )
          base_ks_array%my_k_points( ks )%info%spin      = s
          base_ks_array%my_k_points( ks )%info%k_indices = k_points( :, k )
          Allocate( base_ks_array%my_k_points( ks )%data( 1:1 ) )
          base_ks_array%my_k_points( ks )%data( 1 )%label = 1
          base_ks_array%my_k_points( ks )%data( 1 )%matrix = base_matrix
          base_ks_array%my_k_points( ks )%communicator = base_matrix%get_comm()
       End Do
    End Do

    base_ks_array%parent_communicator = base_matrix%get_comm()

  End Subroutine ks_array_comm_to_base

  Subroutine ks_array_finalise

    Call ks_matrix_finalise
    
  End Subroutine ks_array_finalise

  Subroutine ks_array_create( A, m, n, source )

    ! M and N should be arrays!!!!!!!
    
    ! Create a matrix in all k point mode with no irreps parallelism
    
    ! Also want to think about what kind of object should be used as a source
    
    Class( ks_array ), Intent(   Out ) :: A
    Integer          , Intent( In    ) :: m
    Integer          , Intent( In    ) :: n
    Type ( ks_array ), Intent( In    ) :: source

    Integer :: n_all_ks, n_my_ks
    Integer :: ks
    
    ! Set up the all k point data structure
    n_all_ks = Size( source%all_k_point_info )
    Allocate( A%all_k_point_info( 1:n_all_ks ) )
    A%all_k_point_info = source%all_k_point_info

    ! Now my k points
    n_my_ks = Size( source%my_k_points )
    Allocate( A%my_k_points( 1:n_my_ks ) )
    Do ks = 1, n_my_ks
       A%my_k_points( ks )%info = source%my_k_points( ks )%info
          Allocate( A%my_k_points( ks )%data( 1:1 ) )
          A%my_k_points( ks )%data( 1 )%label = 1
          Call A%my_k_points( ks )%data( 1 )%matrix%create( &
               A%my_k_points( ks )%info%k_type == K_POINT_COMPLEX, &
               A%my_k_points( ks )%info%spin,                      &
               A%my_k_points( ks )%info%k_indices,                 &
               m, n, source%my_k_points( ks )%data( 1 )%matrix )
          A%my_k_points( ks )%communicator = source%my_k_points( ks )%communicator 
    End Do

    A%parent_communicator = source%parent_communicator
    
  End Subroutine ks_array_create

  Subroutine ks_array_print_info( A, name, verbosity )

    Use mpi

    Class( ks_array     ),           Intent( In ) :: A
    Character( Len = * ) , Optional, Intent( In ) :: name
    Integer              , Optional, Intent( In ) :: verbosity

    Integer, Dimension( : ), Allocatable :: this_comm_ranks
    Integer, Dimension( : ), Allocatable :: packed_ranks
    
    Integer :: me, nproc
    Integer :: error
    Integer :: n_ks
    Integer :: n_packed_ranks
    Integer :: ks, my_ks
    Integer :: i

    Logical :: in_range

    Call MPI_Comm_rank( A%parent_communicator, me   , error )
    Call MPI_Comm_size( A%parent_communicator, nproc, error )

    n_ks = Size( A%all_k_point_info )
    
    If( me == 0 ) Then
       Write( *, * )
       If( Present( name ) ) Then
          Write( *, '( a, a )' ) 'Information on ks_array ', name
       Else
          Write( *, '( a )' ) 'Information on a ks_array'
       End If
       Write( *, '( a )' ) 'This ks array holds the following spins and k points'
       Write( *, '( a, t10, a, t30, a )' ) 'Spin', 'Indices', 'Data Type'
       Do ks = 1, n_ks
          Write( *, '( i0, t10, "( ", i2, ", ", i2, ", ", i2, " )", t30, a )' )               &
               A%all_k_point_info( ks )%spin, A%all_k_point_info( ks )%k_indices,             &
               Merge( 'Real   ', 'Complex', A%all_k_point_info( ks )%k_type == K_POINT_REAL )
       End Do
    End If

    ! Print Which procs hold which k point
    If( Present( verbosity ) ) Then
       If( verbosity > 99 ) Then
          If( me == 0 ) Then
             Write( *, '( a )' ) 'The ranks within the parent communicator map onto the k points as below:'
          End If
          ! For each k point in turn work out who owns it
          Do ks = 1, n_ks
             Allocate( this_comm_ranks( 0:nproc - 1 ) )
             this_comm_ranks = 0
             my_ks = A%get_my_ks_index( ks )
             If( my_ks /= NOT_ME ) Then
                this_comm_ranks( me ) = me + 1 ! Avoid problems when me == 0
             End If
             Call MPI_Allreduce( MPI_IN_PLACE, this_comm_ranks, Size( this_comm_ranks ), &
                  MPI_INTEGER, MPI_SUM, A%parent_communicator, error )
             If( me == 0 ) Then
                n_packed_ranks = Size( Pack( this_comm_ranks - 1, this_comm_ranks /= 0 ) )
                Allocate( packed_ranks( 0:n_packed_ranks  -1 ) )
                packed_ranks = Pack( this_comm_ranks - 1, this_comm_ranks /= 0 )
                Write( *, '( i0, t10, "( ", i2, ", ", i2, ", ", i2, " )", 3x )', Advance = 'No' )    &
                     A%all_k_point_info( ks )%spin, A%all_k_point_info( ks )%k_indices
                Write( *, '( i0 )', Advance = 'No' ) packed_ranks( 0 )
                If( Size( packed_ranks ) > 1 ) Then
                   in_range = .False.
                   Do i = Lbound( packed_ranks, Dim = 1 ) + 1, Ubound( packed_ranks, Dim = 1 ) - 1
                      If( packed_ranks( i ) - packed_ranks( i - 1 ) == 1 ) Then
                         If( .Not. in_range ) Then
                            Write( *, '( "-" )', Advance = 'No' )
                         End If
                         in_range = .True.
                      Else
                         If( in_range ) Then
                            Write( *, '( i0, ",", i0 )', Advance = 'No' ) &
                                 packed_ranks( i - 1 ), packed_ranks( i )
                            in_range = .False.
                         Else
                            Write( *, '( ",", i0, "," )', Advance = 'No' ) packed_ranks( i - 1 )
                         End If
                      End If
                   End Do
                   If( in_range ) Then
                      Write( *, '( i0 )' ) packed_ranks( i )
                   Else
                      Write( *, '( ",", i0 )' ) packed_ranks( i )
                   End If
                Else
                   Write( *, * )
                End If
                Deallocate( packed_ranks )
             End If
             Deallocate( this_comm_ranks )
          End Do
          If( me == 0 ) Then
             Write( *, * )
          End If
       End If
    End If
    
  End Subroutine ks_array_print_info

  Subroutine ks_array_split_ks( A, complex_weight, split_A, redistribute )

    ! Split a ks_array A so the resulting ks_array is k point distributed

    Use mpi
    
    Class( ks_array     ), Intent( In    ) :: A
    Real ( wp           ), Intent( In    ) :: complex_weight
    Type ( ks_array     ), Intent(   Out ) :: split_A
    Logical, Optional    , Intent( In    ) :: redistribute

    Type( ks_matrix ), Allocatable :: this_ks_matrix
    Type( ks_matrix ), Allocatable :: split_ks_matrix

    Type( ks_matrix ) :: base_matrix

    Real( wp ), Dimension( : ), Allocatable :: weights

    Integer, Dimension( : ), Allocatable :: n_procs_ks
    Integer, Dimension( : ), Allocatable :: my_ks
    Integer, Dimension( : ), Allocatable :: this_k_indices

    Integer :: m, n
    Integer :: n_ks
    Integer :: n_procs_parent, me_parent, my_colour, k_comm, n_my_ks
    Integer :: this_k_type, this_s
    Integer :: top_rank
    Integer :: cost
    Integer :: error
    Integer :: ks, this_ks

    Logical :: loc_redist

    If( Present( redistribute ) ) Then
       loc_redist = redistribute
    Else
       loc_redist =.True.
    End If

    ! Set up stuff relevant to all k points
    split_A%all_k_point_info    = A%all_k_point_info
    split_A%parent_communicator = A%parent_communicator

    ! Now split the k points and return the new structure in split_A

    ! First split the communicator
    ! Generate an array for the weights
    n_ks = Size( split_A%all_k_point_info )
    Allocate( weights( 1:n_ks ) )
    Do ks = 1, n_ks
       weights( ks ) = Merge( 1.0_wp, complex_weight, split_A%all_k_point_info( ks )%k_type == K_POINT_REAL )
    End Do
    ! Makes a difference if somehow somebody has all k points complex
    weights = weights / Minval( weights )

    !THIS WHOLE SPLITTING STRATEGY PROBABLY NEEDS MORE THOUGHT
    Call MPI_Comm_size( split_A%parent_communicator, n_procs_parent, error )
    Call MPI_Comm_rank( split_A%parent_communicator,      me_parent, error )
    cost = Nint( Sum( weights ) )
    ! Scale up weights so fits inside the parent comm
    Allocate( n_procs_ks( 1:n_ks ) )
    n_procs_ks = Nint( weights * ( n_procs_parent / cost ) )

    ! Two possible cases

    k_split_strategy: If( cost <= n_procs_parent ) Then
       
       ! 1) There are sufficent processors for each k point to have its own, separate set if processors which work on it
       n_my_ks  = 1
       Allocate( my_ks( 1:n_my_ks ) ) 

       ! Decide which group ( if any ) I am in can probably write this more neatly but lets keep
       ! it explicit as I'm getting a little confused
       If( me_parent > Sum( n_procs_ks ) - 1 ) Then
          my_colour  = MPI_UNDEFINED
          my_ks( 1 ) = INVALID
          n_my_ks    = 0
       Else
          top_rank = 0
          Do ks = 1, n_ks
             top_rank = top_rank + n_procs_ks( ks )
             If( me_parent < top_rank ) Then
                ! Colour for comm spliting
                my_colour = ks
                ! As in this strategy we only have 1 k point store which one it is
                my_ks( 1 ) = ks
                Exit
             End If
          End Do
       End If

    Else

       ! 2) Not enough procs to make distibuting matrices wortwhile. Just use a simple
       ! round robin assignment and do operations in serial

       ! First work out how many ks points I will hold
       n_my_ks = n_ks / n_procs_parent
       If( n_my_ks * n_procs_parent < n_ks ) Then
          n_my_ks = n_my_ks + 1
       End If
       Allocate( my_ks( 1:n_my_ks ) )

       ! Now assign them
       this_ks = 0
       Do ks = me_parent + 1, n_ks, n_procs_parent
          this_ks = this_ks + 1
          my_ks( this_ks ) = ks
       End Do

       ! And colour my communicators, an mpi_comm_split is next
       my_colour = Merge( me_parent, MPI_UNDEFINED, me_parent < n_ks )

    End If k_split_strategy

    ! Can now split the communicator
    Call MPI_Comm_split( split_A%parent_communicator, my_colour, 0, k_comm, error )

    ! Now start setting up the k points held by this set of processes (if any!)
    Allocate( split_A%my_k_points( 1:n_my_ks ) )
    Allocate( this_k_indices( 1:Size( A%all_k_point_info( 1 )%k_indices ) ) )
    Do ks = 1, n_my_ks
       split_A%my_k_points( ks )%info = split_A%all_k_point_info( my_ks( ks ) )
       ! Irreps not split yet hence no split at this level
       Allocate( split_A%my_k_points( ks )%data( 1:1 ) )
       split_A%my_k_points( ks )%data( 1 )%label = 1
       this_k_type    = split_A%all_k_point_info( my_ks( ks ) )%k_type
       this_s         = split_A%all_k_point_info( my_ks( ks ) )%spin
       this_k_indices = split_A%all_k_point_info( my_ks( ks ) )%k_indices
       ! Now need to generate a source matrix from the communicator - precisely what the init routine does!!
       Call ks_matrix_comm_to_base( k_comm, base_matrix )
       ! Need to get sizes for creation
       m = A%my_k_points( my_ks( ks ) )%data( 1 )%matrix%size( 1 )
       n = A%my_k_points( my_ks( ks ) )%data( 1 )%matrix%size( 2 )
       Call split_A%my_k_points( ks )%data( 1 )%matrix%create( this_k_type == K_POINT_COMPLEX, &
            this_s, this_k_indices, m, n, base_matrix )
       split_A%my_k_points( ks )%communicator = base_matrix%get_comm()
    End Do

    ! Finally if required redistribute the data from the all procs distribution into the split distribution
    ! Note all procs in parent comm must be involved as all hold data for the unsplit matrx
    If( loc_redist ) Then
       Allocate( this_ks_matrix )
       Do ks = 1, n_ks
          this_ks_matrix = A%my_k_points( ks )%data( 1 )%matrix
          ! Find out if I hold this in the split distribution
          ! Note indicate to the remap routine that we hod no data by having an unallocated array
          ! c.f. Null pointer in C type languages
          this_ks = split_A%get_my_ks_index( ks )
          If( this_ks /= NOT_ME ) Then
             Allocate( split_ks_matrix )
             split_ks_matrix = split_A%my_k_points( this_ks )%data( 1 )%matrix
          End If
          Call ks_matrix_remap_data( A%parent_communicator, this_ks_matrix, split_ks_matrix )
          If( this_ks /= NOT_ME ) Then
             split_A%my_k_points( this_ks )%data( 1 )%matrix = split_ks_matrix
             Deallocate( split_ks_matrix ) ! Important as an deallocated matrix indicates no data on this process in the remap routine
          End If
       End Do
       Deallocate( this_ks_matrix )
    End If

  End Subroutine ks_array_split_ks

  Subroutine ks_array_diag( A, Q, E )

    Use mpi

    Class( ks_array     ),                 Intent( In    ) :: A
    Type ( ks_array     ),                 Intent(   Out ) :: Q
    Type ( eval_storage ), Dimension( : ), Intent(   Out ) :: E

    Real( wp ) :: rdum

    Integer, Dimension( 1:2 ) :: buff_send, buff_recv
    
    Integer :: me, me_parent
    Integer :: ks_root, nb
    Integer :: error
    Integer :: request
    Integer :: rsize, handle
    Integer :: my_ks, ks
    Integer :: my_irrep

    Logical :: sending_data
    
    ! Make Q have the same set up as A
    Q = A

    Do my_ks = 1, Size( A%my_k_points )
       ! Irreps will need more thought - worrk currenly as burnt into as 1
       Do my_irrep = 1, Size( A%my_k_points( my_ks )%data )
          ks = A%get_all_ks_index( my_ks )
          Associate( Aks => A%my_k_points( my_ks )%data( my_irrep )%matrix, &
                     Qks => Q%my_k_points( my_ks )%data( my_irrep )%matrix )
            E( ks )%spin      = A%all_k_point_info( ks )%spin
            E( ks )%k_indices = A%all_k_point_info( ks )%k_indices
            Call Aks%diag( Qks, E( ks )%evals )
          End Associate
       End Do
    End Do

    ! Replicate evals
    ! Again needs thought for ireps
    Call mpi_comm_rank( A%parent_communicator, me_parent, error )
    Do ks = 1, Size( A%all_k_point_info )
       my_ks = A%get_my_ks_index( ks )
       ! Work out who holds this set of evals and send how many there are
       ! the root node of the communicator holding them back to the root node of the parent communicator
       sending_data = .False.
       If( my_ks /= NOT_ME ) Then
          Call mpi_comm_rank( A%my_k_points( my_ks )%communicator, me, error )
          sending_data = me == 0
          If( sending_data ) then
             buff_send( 1 ) = me_parent
             buff_send( 2 ) = Size( E( ks )%evals )
             Call mpi_isend( buff_send, 2, MPI_INTEGER, 0, ks, A%parent_communicator, request, error )
          End If
       End If
       If( me_parent == 0 ) Then
          Call mpi_recv( buff_recv, 2, MPI_INTEGER, MPI_ANY_SOURCE, ks, A%parent_communicator, MPI_STATUS_IGNORE, error )
       End If
       If( sending_data ) Then
          Call mpi_wait( request, MPI_STATUS_IGNORE, error )
       End If
       ! Now on root of parent bcast back to all
       Call mpi_bcast( buff_recv, 2, MPI_INTEGER, 0, A%parent_communicator, error )
       ks_root = buff_recv( 1 )
       nb      = buff_recv( 2 )
       ! Now know how many evals we will recv - allocate memory if haven't done so already because I don't 'own' this k point
       If( .Not. Allocated( E( ks )%evals ) ) Then
          Allocate( E( ks )%evals( 1:nb ) )
       End If
       ! And finally bcast out the values from the root node for this set of evals
       Call mpi_sizeof( rdum, rsize, error )
       Call mpi_type_match_size( MPI_TYPECLASS_REAL, rsize, handle, error )
       Call mpi_bcast( E( ks )%evals, Size( E( ks )%evals ), handle, ks_root, A%parent_communicator, error )
    End Do
    
  End Subroutine ks_array_diag

  Function ks_array_dagger( A ) Result( tA )

    Type( ks_array ), Allocatable :: tA

    Class( ks_array ), Intent( In ) :: A

    Integer :: my_ks, my_irrep

    Allocate( tA )
    tA = A
    
    Do my_ks = 1, Size( A%my_k_points )
       ! Irreps will need more thought - worrk currenly as burnt into as 1
       Do my_irrep = 1, Size( A%my_k_points( my_ks )%data )
          Associate( Aks  =>  A%my_k_points( my_ks )%data( my_irrep )%matrix, &
                     tAks => tA%my_k_points( my_ks )%data( my_irrep )%matrix )
            tAks = .Dagger. Aks
          End Associate
       End Do
    End Do

  End Function ks_array_dagger

  Function ks_array_mult( A, B ) Result( C )

    Type( ks_array ), Allocatable :: C

    Class( ks_array ), Intent( In ) :: A
    Type ( ks_array ), Intent( In ) :: B

    Integer :: my_ks, my_irrep

    Allocate( C )
    C = A
    
    Do my_ks = 1, Size( A%my_k_points )
       ! Irreps will need more thought - work currenly as burnt into as 1
       Do my_irrep = 1, Size( A%my_k_points( my_ks )%data )
          Associate( Aks => A%my_k_points( my_ks )%data( my_irrep )%matrix, &
                     Bks => B%my_k_points( my_ks )%data( my_irrep )%matrix, &
                     Cks => C%my_k_points( my_ks )%data( my_irrep )%matrix )
            Cks = Aks * Bks
          End Associate
       End Do
    End Do

  End Function ks_array_mult

  Function ks_array_pre_scale( s, A ) Result( C )

    Type( ks_array ), Allocatable :: C

    Real( wp )       , Intent( In ) :: s
    Class( ks_array ), Intent( In ) :: A

    Integer :: my_ks, my_irrep

    Allocate( C )
    C = A
    
    Do my_ks = 1, Size( A%my_k_points )
       ! Irreps will need more thought - work currenly as burnt into as 1
       Do my_irrep = 1, Size( A%my_k_points( my_ks )%data )
          Associate( Aks => A%my_k_points( my_ks )%data( my_irrep )%matrix, &
                     Cks => C%my_k_points( my_ks )%data( my_irrep )%matrix )
            Cks = s * Aks
          End Associate
       End Do
    End Do

  End Function ks_array_pre_scale

  Function ks_array_post_scale( A, s ) Result( C )

    Type( ks_array ), Allocatable :: C

    Class( ks_array ), Intent( In ) :: A
    Real( wp )       , Intent( In ) :: s

    Integer :: my_ks, my_irrep

    Allocate( C )
    C = A
    
    Do my_ks = 1, Size( A%my_k_points )
       ! Irreps will need more thought - work currenly as burnt into as 1
       Do my_irrep = 1, Size( A%my_k_points( my_ks )%data )
          Associate( Aks => A%my_k_points( my_ks )%data( my_irrep )%matrix, &
                     Cks => C%my_k_points( my_ks )%data( my_irrep )%matrix )
            Cks = Aks * s
          End Associate
       End Do
    End Do

  End Function ks_array_post_scale

  Function ks_array_pre_mult_diag( s, A ) Result( C )

    Type( ks_array ), Allocatable :: C

    Real( wp )       , Dimension( : ), Intent( In ) :: s
    Class( ks_array ),                Intent( In ) :: A

    Integer :: my_ks, my_irrep

    Allocate( C )
    C = A
    
    Do my_ks = 1, Size( A%my_k_points )
       ! Irreps will need more thought - work currenly as burnt into as 1
       Do my_irrep = 1, Size( A%my_k_points( my_ks )%data )
          Associate( Aks => A%my_k_points( my_ks )%data( my_irrep )%matrix, &
                     Cks => C%my_k_points( my_ks )%data( my_irrep )%matrix )
            Cks = s * Aks
          End Associate
       End Do
    End Do

  End Function ks_array_pre_mult_diag

  Function ks_array_post_mult_diag( A, s ) Result( C )

    Type( ks_array ), Allocatable :: C

    Class( ks_array ),                Intent( In ) :: A
    Real( wp )       , Dimension( : ), Intent( In ) :: s

    Integer :: my_ks, my_irrep

    Allocate( C )
    C = A
    
    Do my_ks = 1, Size( A%my_k_points )
       ! Irreps will need more thought - work currenly as burnt into as 1
       Do my_irrep = 1, Size( A%my_k_points( my_ks )%data )
          Associate( Aks => A%my_k_points( my_ks )%data( my_irrep )%matrix, &
                     Cks => C%my_k_points( my_ks )%data( my_irrep )%matrix )
            Cks = Aks * s
          End Associate
       End Do
    End Do

  End Function ks_array_post_mult_diag

  Function ks_array_add( A, B ) Result( C )

    Type( ks_array ), Allocatable :: C

    Class( ks_array ), Intent( In ) :: A
    Type ( ks_array ), Intent( In ) :: B

    Integer :: my_ks, my_irrep

    Allocate( C )
    C = A
    
    Do my_ks = 1, Size( A%my_k_points )
       ! Irreps will need more thought - work currenly as burnt into as 1
       Do my_irrep = 1, Size( A%my_k_points( my_ks )%data )
          Associate( Aks => A%my_k_points( my_ks )%data( my_irrep )%matrix, &
                     Bks => B%my_k_points( my_ks )%data( my_irrep )%matrix, &
                     Cks => C%my_k_points( my_ks )%data( my_irrep )%matrix )
            Cks = Aks + Bks
          End Associate
       End Do
    End Do

  End Function ks_array_add

  Pure Function get_all_ks_index( A, my_ks ) Result( ks )

    Integer :: ks

    Class( ks_array ), Intent( In ) :: A
    Integer          , Intent( In ) :: my_ks

    Do ks = 1, Size( A%all_k_point_info )
       If( A%all_k_point_info( ks )%spin == A%my_k_points( my_ks )%info%spin ) Then
          If( All( A%all_k_point_info( ks )%k_indices == A%my_k_points( my_ks )%info%k_indices ) ) Then
             Exit
          End If
       End If
    End Do
    
  End Function get_all_ks_index

  Pure Function get_my_ks_index( A, ks ) Result( my_ks )

    Integer :: my_ks

    Class( ks_array ), Intent( In ) :: A
    Integer          , Intent( In ) :: ks

    Integer :: k
    
    my_ks = NOT_ME

    Do k = 1, Size( A%my_k_points )
       If( A%all_k_point_info( ks )%spin == A%my_k_points( k )%info%spin ) Then
          If( All( A%all_k_point_info( ks )%k_indices == A%my_k_points( k )%info%k_indices ) ) Then
             my_ks = k
             Exit
          End If
       End If
    End Do

  End Function get_my_ks_index

  Pure Function get_ks( A, k, s ) Result( ks )

    Integer :: ks

    Class( ks_array )        , Intent( In ) :: A
    Integer, Dimension( 1:3 ), Intent( In ) :: k
    Integer                  , Intent( In ) :: s

    Do ks = 1, Size( A%all_k_point_info )
       If( A%all_k_point_info( ks )%spin == s ) Then
          If( All( A%all_k_point_info( ks )%k_indices == k ) ) Then
             Exit
          End If
       End If
    End Do
    
  End Function get_ks

  Subroutine ks_array_set_global_real( A, s, k, m, n, p, q, data )

    ! Need to overload for irreps

    Class( ks_array )              , Intent( InOut ) :: A
    Integer                        , Intent( In    ) :: s
    Integer    , Dimension( : )    , Intent( In    ) :: k
    Integer                        , Intent( In    ) :: m
    Integer                        , Intent( In    ) :: n
    Integer                        , Intent( In    ) :: p
    Integer                        , Intent( In    ) :: q
    Real( wp ), Dimension( m:, p: ), Intent( In    ) :: data

    Integer :: ks, my_ks

    ks = A%get_ks( k, s )
    
    my_ks = A%get_my_ks_index( ks )

    If( my_ks /= NOT_ME ) Then
       Call A%my_k_points( my_ks )%data( 1 )%matrix%set_by_global( m, n, p, q, data )
    End If

  End Subroutine ks_array_set_global_real

  Subroutine ks_array_set_global_complex( A, s, k, m, n, p, q, data )

    ! Need to overload for irreps

    Class( ks_array )                 , Intent( InOut ) :: A
    Integer                           , Intent( In    ) :: s
    Integer    , Dimension( : )       , Intent( In    ) :: k
    Integer                           , Intent( In    ) :: m
    Integer                           , Intent( In    ) :: n
    Integer                           , Intent( In    ) :: p
    Integer                           , Intent( In    ) :: q
    Complex( wp ), Dimension( m:, p: ), Intent( In    ) :: data

    Integer :: ks, my_ks

    ks = A%get_ks( k, s )
    
    my_ks = A%get_my_ks_index( ks )

    If( my_ks /= NOT_ME ) Then
       Call A%my_k_points( my_ks )%data( 1 )%matrix%set_by_global( m, n, p, q, data )
    End If

  End Subroutine ks_array_set_global_complex  

  Subroutine ks_array_get_global_real( A, s, k, m, n, p, q, data )

    Use mpi

    ! Need to overload for irreps

    Class( ks_array )              , Intent( In    ) :: A
    Integer                        , Intent( In    ) :: s
    Integer    , Dimension( : )    , Intent( In    ) :: k
    Integer                        , Intent( In    ) :: m
    Integer                        , Intent( In    ) :: n
    Integer                        , Intent( In    ) :: p
    Integer                        , Intent( In    ) :: q
    Real( wp ), Dimension( m:, p: ), Intent(   Out ) :: data

    Real( wp ) :: rdum

    Integer :: me_parent, me
    Integer :: buff_send, buff_recv
    Integer :: ks, my_ks
    Integer :: request
    Integer :: ks_root
    Integer :: rsize, handle
    Integer :: error

    Logical :: sending_data
    
    ks = A%get_ks( k, s )
    
    my_ks = A%get_my_ks_index( ks )

    If( my_ks /= NOT_ME ) Then
       Call A%my_k_points( my_ks )%data( 1 )%matrix%get_by_global( m, n, p, q, data )
    End If

    ! Need to replicate data over parent communicator if split ks points
    
    ! Work out who holds this set of evals and send how many there are
    ! the root node of the communicator holding them back to the root node of the parent communicator
    Call mpi_comm_rank( A%parent_communicator, me_parent, error )
    sending_data = .False.
    If( my_ks /= NOT_ME ) Then
       Call mpi_comm_rank( A%my_k_points( my_ks )%communicator, me, error )
       sending_data = me == 0
       If( sending_data ) then
          buff_send = me_parent
          Call mpi_isend( buff_send, 1, MPI_INTEGER, 0, ks, A%parent_communicator, request, error )
       End If
    End If
    If( me_parent == 0 ) Then
       Call mpi_recv( buff_recv, 1, MPI_INTEGER, MPI_ANY_SOURCE, ks, A%parent_communicator, MPI_STATUS_IGNORE, error )
    End If
    If( sending_data ) Then
       Call mpi_wait( request, MPI_STATUS_IGNORE, error )
    End If
    ! Now on root of parent bcast back to all
    Call mpi_bcast( buff_recv, 1, MPI_INTEGER, 0, A%parent_communicator, error )
    ks_root = buff_recv 
    Call mpi_sizeof( rdum, rsize, error )
    Call mpi_type_match_size( MPI_TYPECLASS_REAL, rsize, handle, error )
    Call mpi_bcast( data, Size( data ), handle, ks_root, A%parent_communicator, error )    

  End Subroutine ks_array_get_global_real

  Subroutine ks_array_get_global_complex( A, s, k, m, n, p, q, data )

    Use mpi

    ! Need to overload for irreps

    Class( ks_array )                 , Intent( In    ) :: A
    Integer                           , Intent( In    ) :: s
    Integer    , Dimension( : )       , Intent( In    ) :: k
    Integer                           , Intent( In    ) :: m
    Integer                           , Intent( In    ) :: n
    Integer                           , Intent( In    ) :: p
    Integer                           , Intent( In    ) :: q
    Complex( wp ), Dimension( m:, p: ), Intent(   Out ) :: data

    Complex( wp ) :: cdum

    Integer :: me_parent, me
    Integer :: buff_send, buff_recv
    Integer :: ks, my_ks
    Integer :: request
    Integer :: ks_root
    Integer :: rsize, handle
    Integer :: error

    Logical :: sending_data
    
    ks = A%get_ks( k, s )
    
    my_ks = A%get_my_ks_index( ks )

    If( my_ks /= NOT_ME ) Then
       Call A%my_k_points( my_ks )%data( 1 )%matrix%get_by_global( m, n, p, q, data )
    End If

    ! Need to replicate data over parent communicator if split ks points
    
    ! Work out who holds this set of evals and send how many there are
    ! the root node of the communicator holding them back to the root node of the parent communicator
    Call mpi_comm_rank( A%parent_communicator, me_parent, error )
    sending_data = .False.
    If( my_ks /= NOT_ME ) Then
       Call mpi_comm_rank( A%my_k_points( my_ks )%communicator, me, error )
       sending_data = me == 0
       If( sending_data ) then
          buff_send = me_parent
          Call mpi_isend( buff_send, 1, MPI_INTEGER, 0, ks, A%parent_communicator, request, error )
       End If
    End If
    If( me_parent == 0 ) Then
       Call mpi_recv( buff_recv, 1, MPI_INTEGER, MPI_ANY_SOURCE, ks, A%parent_communicator, MPI_STATUS_IGNORE, error )
    End If
    If( sending_data ) Then
       Call mpi_wait( request, MPI_STATUS_IGNORE, error )
    End If
    ! Now on root of parent bcast back to all
    Call mpi_bcast( buff_recv, 1, MPI_INTEGER, 0, A%parent_communicator, error )
    ks_root = buff_recv 
    Call mpi_sizeof( cdum, rsize, error )
    Call mpi_type_match_size( MPI_TYPECLASS_COMPLEX, rsize, handle, error )
    Call mpi_bcast( data, Size( data ), handle, ks_root, A%parent_communicator, error )    

  End Subroutine ks_array_get_global_complex

End Module ks_array_module
