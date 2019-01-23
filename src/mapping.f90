Module mapping_module

  Implicit None

  Integer, Parameter, Private :: INVALID = -1

  Type, Public :: mapping
     Character( Len = 128 )   , Private :: name
     Integer                  , Private :: communicator
     Integer, Dimension( 1:9 ), Private :: descriptor
     Integer                  , Private :: parent_communicator
     Integer, Dimension( 1:9 ), Private :: parent_descriptor
   Contains
     Procedure :: print 
     Procedure :: split
  End Type mapping

  Type( mapping ), Public :: mapping_base_map = mapping( name = 'BASE_MAP', communicator = INVALID, descriptor = INVALID, &
       parent_communicator = INVALID, parent_descriptor = INVALID )      

  Integer, Parameter, Public :: mapping_get_global_n_row = 3
  Integer, Parameter, Public :: mapping_get_global_n_col = 4

  Public :: mapping_init
  Public :: mapping_finalise
  Public :: mapping_get_data
  
  Private

Contains

  Subroutine mapping_init( comm )

    Integer, Intent( In ) :: comm

    Integer :: base_context
    
    mapping_base_map%communicator = comm

    ! Get base context here ...
    base_context = INVALID - 1
    mapping_base_map%descriptor( 2 ) = base_context
    
  End Subroutine mapping_init

  Subroutine mapping_get_data( map, what, data )

    ! get out of BLACS descriptor - only so far for
    ! dense in core matrix, but easily modifiable for other type

    Type( mapping ), Intent( In    ) :: map
    Integer        , Intent( In    ) :: what
    Integer        , Intent(   Out ) :: data

    Select Case( what )

    Case( mapping_get_global_n_row )
       data = map%descriptor( 3 )
    
    Case( mapping_get_global_n_col )
       data = map%descriptor( 4 )

    Case Default
       Stop 'Illegal WHAT in get_mapping_data'

    End Select
    
  End Subroutine mapping_get_data

  Subroutine mapping_finalise
    
    mapping_base_map = mapping( name = 'BASE_MAP', communicator = INVALID, descriptor = INVALID, &
         parent_communicator = INVALID, parent_descriptor = INVALID )      

  End Subroutine mapping_finalise

  Subroutine print( map )

    Use mpi
    
    Class( mapping ), Intent( In ) :: map

    Integer :: rank, nproc, nproc_parent
    Integer :: error

    Call mpi_comm_size( map%communicator, nproc, error )
    Call mpi_comm_rank( map%communicator, rank , error )
    If( rank == 0 ) Then
       Write( *, '( a, a, a, i0 )' ) 'Size of mapping ', Trim( Adjustl( map%name ) ), ' is ', nproc
       ! Catch base map - can probably do better
       If( map%parent_communicator /= INVALID ) Then
          Call mpi_comm_size( map%communicator, nproc_parent, error )
          Write( *, '( a, a, a, i0 )' ) 'Size of parent of mapping ', Trim( Adjustl( map%name ) ), ' is ', nproc_parent
       End If
    End If
    
  End Subroutine print

  Subroutine split( map, weights, split_map, i_hold )

    Use mpi

    Class( mapping )                    , Intent( In    ) :: map
    Integer, Dimension( : )             , Intent( In    ) :: weights
    Class( mapping )                    , Intent(   Out ) :: split_map
    Integer, Dimension( : ), Allocatable, Intent(   Out ) :: i_hold

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
    Else
       colour = MPI_UNDEFINED
    End If

    Call mpi_comm_split( map%communicator, colour, 1, split_comm, error )
    If( colour /= MPI_UNDEFINED ) Then
       Call mpi_comm_size( split_comm, nproc, error )
       Call mpi_comm_rank( split_comm, rank , error )
    Else
       rank  = INVALID
       nproc = INVALID
    End If

    Write( *, * ) parent_nproc, parent_rank, nproc, rank, colour

    
  End Subroutine split

End Module mapping_module
