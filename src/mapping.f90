Module mapping_module

  Implicit None

  Integer, Parameter, Private :: INVALID = -1

  Type, Public :: mapping
     Integer                  , Private :: communicator
     Integer, Dimension( 1:9 ), Private :: descriptor
  End Type mapping

  Type( mapping ), Public :: mapping_base_map = mapping( communicator = INVALID, descriptor = INVALID )

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
    mapping_base_map = mapping( communicator = INVALID, descriptor = INVALID )
  End Subroutine mapping_finalise

End Module mapping_module
