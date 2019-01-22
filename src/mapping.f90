Module mapping_module

  Implicit None

  Type, Public :: mapping
     Integer                             , Private :: communicator
     Integer, Dimension( : ), Allocatable, Private :: descriptor
  End Type mapping

  Integer, Parameter, Public :: mapping_get_global_n_row = 3
  Integer, Parameter, Public :: mapping_get_global_n_col = 4

  Public :: mapping_init
  Public :: mapping_get_data
  
  Private

Contains

  Subroutine mapping_init( n, bfac, comm, map )

    Integer        , Intent( In    ) :: n
    Integer        , Intent( In    ) :: bfac
    Integer        , Intent( In    ) :: comm
    Type( mapping ), Intent(   Out ) :: map

    map%communicator = comm
    
    Allocate( map%descriptor( 1:9 ) )
    ! At some point make consistent with parameters above
    map%descriptor( 1 ) = 1
    map%descriptor( 3 ) = n
    map%descriptor( 4 ) = n
    map%descriptor( 5 ) = bfac
    map%descriptor( 6 ) = bfac

    
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

End Module mapping_module
