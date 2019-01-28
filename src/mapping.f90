Module mapping_module

  Use mpi

  Implicit None

  Integer, Parameter, Private :: INVALID = -1

  Type, Public :: mapping
     Character( Len = 128 )   , Private :: name
     Integer                  , Private :: communicator
!!$     Integer, Dimension( 1:9 ), Private :: descriptor
     Integer                  , Private :: parent_communicator
!!$     Integer, Dimension( 1:9 ), Private :: parent_descriptor
   Contains
     Procedure :: set   => set_mapping
     Procedure :: print => print_mapping
     Procedure :: split => split_mapping
  End Type mapping

  Type( mapping ), Public :: mapping_base_map = mapping( name = 'BASE_MAP', &
       communicator = MPI_COMM_NULL, parent_communicator = MPI_COMM_NULL )      

!!$  Integer, Parameter, Public :: mapping_get_global_n_row = 3
!!$  Integer, Parameter, Public :: mapping_get_global_n_col = 4

  Public :: mapping_init
  Public :: mapping_finalise
!!$  Public :: mapping_get_data
  
  Private

Contains

  Subroutine mapping_init( comm )

    Integer, Intent( In ) :: comm

!!$    Integer, Dimension( 1:9 ) :: descriptor
!!$    Integer, Dimension( 1:9 ) :: parent_descriptor

    Integer :: parent_communicator
    
!!$    descriptor = INVALID

    parent_communicator = MPI_COMM_NULL
!!$    parent_descriptor   = INVALID

    Call mapping_base_map%set( 'BASE_MAP', comm, parent_communicator )
    
  End Subroutine mapping_init

  Subroutine mapping_finalise
    
    Integer :: communicator
    Integer :: parent_communicator

    communicator = MPI_COMM_NULL

    parent_communicator = MPI_COMM_NULL

    Call mapping_base_map%set( 'BASE_MAP', communicator, parent_communicator )

  End Subroutine mapping_finalise

  Subroutine print_mapping( map )

    Use mpi
    
    Class( mapping ), Intent( In ) :: map

    Integer :: rank, nproc, nproc_parent
    Integer :: error

    Call mpi_comm_size( map%communicator, nproc, error )
    Call mpi_comm_rank( map%communicator, rank , error )
    If( rank == 0 ) Then
       Write( *, '( a, a, a, i0 )' ) 'Size of mapping ', Trim( Adjustl( map%name ) ), ' is ', nproc
       ! Catch base map - can probably do better
       If( map%parent_communicator /= MPI_COMM_NULL ) Then
          Call mpi_comm_size( map%communicator, nproc_parent, error )
          Write( *, '( a, a, a, i0 )' ) 'Size of parent of mapping ', Trim( Adjustl( map%name ) ), ' is ', nproc_parent
       End If
    End If
    
  End Subroutine print_mapping

  Subroutine split_mapping( map, weights, split_name, split_map, i_hold )

    ! Totally hacked together at the moment - need to think of a better splitting algorithm -
    ! maybe sort weights and then grab proc in decreasing size and use that to grab procs

    Use mpi

    Class( mapping )                                    , Intent( In    ) :: map
    Integer, Dimension( : )                             , Intent( In    ) :: weights
    Character( Len = * )                                , Intent( In    ) :: split_name
    Class( mapping )       , Dimension( : ), Allocatable, Intent(   Out ) :: split_map
    Integer, Dimension( : ),                 Allocatable, Intent(   Out ) :: i_hold

    Integer :: base_cost
    Integer :: rank, nproc
    Integer :: parent_rank, parent_nproc
    Integer :: n_split
    Integer :: split_comm
    Integer :: colour
    Integer :: error
    Integer :: i, j

    Call mpi_comm_size( map%communicator, parent_nproc, error )
    Call mpi_comm_rank( map%communicator, parent_rank , error )
    
    base_cost = Sum( weights )

    n_split = parent_nproc / base_cost

    ! This needs to be generalised!!!
    If( parent_rank < base_cost * n_split ) Then 
       colour = 0
       rank   = 0 
       Outer: Do i = 1, Size( weights )
          colour = colour + 1
          Do j = 1, n_split * weights( i )
             If( rank == parent_rank ) Then
                Exit Outer
             End If
             rank = rank + 1
          End Do
       End Do Outer
       Allocate( i_hold( 1:1 ) )
       Allocate( split_map( 1:1 ) )
       i_hold( 1 ) = i
    Else
       colour = MPI_UNDEFINED
       Allocate( i_hold( 1:0 ) )
       Allocate( split_map( 1:0 ) )
    End If

    Call mpi_comm_split( map%communicator, colour, 1, split_comm, error )
    If( split_comm /= MPI_COMM_NULL ) Then
       Call mpi_comm_size( split_comm, nproc, error )
       Call mpi_comm_rank( split_comm, rank , error )
    Else
       rank  = INVALID
       nproc = INVALID
    End If

    Write( *, * ) parent_nproc, parent_rank, nproc, rank, colour
    
    Call split_map( 1 )%set( split_name, split_comm, map%communicator)

  End Subroutine split_mapping

  Subroutine set_mapping( map, name, communicator, parent_communicator )

    Class( mapping )       , Intent(   Out ) :: map
    Character( Len = * )   , Intent( In    ) :: name
    Integer                , Intent( In    ) :: communicator
    Integer                , Intent( In    ) :: parent_communicator

    map%name                = name
    map%communicator        = communicator
    map%parent_communicator = parent_communicator

  End Subroutine set_mapping
  
End Module mapping_module
