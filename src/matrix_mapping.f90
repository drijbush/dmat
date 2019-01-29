Module matrix_mapping_module

  Use numbers_module     , Only : wp
  Use mpi
  Use proc_mapping_module, Only : proc_mapping, proc_mapping_init, proc_mapping_finalise, proc_mapping_base, proc_mapping_base_start
  
  Implicit None

  Type, Public, Extends( proc_mapping )  :: matrix_mapping
     Integer, Dimension( 1:9 ), Private :: descriptor
   Contains
     Procedure, Private :: set_matrix_mapping
     Procedure, Private :: split_matrix_mapping
     Procedure, Public  :: print => print_matrix_mapping
     Generic  , Public  :: set   => set_matrix_mapping
     Generic  , Public  :: split => split_matrix_mapping
  End type matrix_mapping

  Integer, Parameter, Private :: INVALID = -1

  Type( matrix_mapping ), Private, Parameter :: matrix_mapping_base_start = &
       matrix_mapping( descriptor = [ INVALID, INVALID, INVALID,     &
                                      INVALID, INVALID, INVALID,     &
                                      INVALID, INVALID, INVALID ],   &
                                      proc_mapping =  proc_mapping_base_start )
  Type( matrix_mapping ), Public, Protected :: matrix_mapping_base = matrix_mapping_base_start

  Public :: matrix_mapping_init
  Public :: matrix_mapping_finalise
  
  Private

  ! What the elements of the descriptor mean
  Integer, Parameter, Private :: dtype_a = 1 ! descriptor type
  Integer, Parameter, Private :: ctxt_a  = 2 ! context
  Integer, Parameter, Private :: m_a     = 3 ! number of global rows
  Integer, Parameter, Private :: n_a     = 4 ! number of global columns
  Integer, Parameter, Private :: mb_a    = 5 ! row blocking factor
  Integer, Parameter, Private :: nb_a    = 6 ! column blocking factor
  Integer, Parameter, Private :: rsrc_a  = 7 ! first process row which holds a
  Integer, Parameter, Private :: csrc_a  = 8 ! first process col which golds a
  Integer, Parameter, Private :: lld_a   = 9 ! leading dimension of LOCAL a

Contains

  Subroutine matrix_mapping_init( comm )

    Integer, Intent( In ) :: comm

    Call proc_mapping_init( comm )

    matrix_mapping_base%descriptor   = INVALID
    matrix_mapping_base%proc_mapping = proc_mapping_base
    
  End Subroutine matrix_mapping_init

  Subroutine matrix_mapping_finalise

    matrix_mapping_base = matrix_mapping_base_start
    
    Call proc_mapping_finalise
    
  End Subroutine matrix_mapping_finalise

  Subroutine print_matrix_mapping( map )

    Class( matrix_mapping ), Intent( In ) :: map

    Integer :: rank
    Integer :: error
    
    Call mpi_comm_rank( map%get_comm(), rank, error )
    If( rank == 0 ) Then
       Write( *, '( a, 9( i0, 1x ) )' ) 'Descriptor for matrix mapping ', map%descriptor
       Call map%proc_mapping%print
    End If
    
  End Subroutine print_matrix_mapping

  Subroutine set_matrix_mapping( map, proc_map, m, n, mb, nb, rsrc, csrc )

    Class( matrix_mapping ), Intent(   Out ) :: map
    Type ( proc_mapping   ), Intent( In    ) :: proc_map
    Integer                , Intent( In    ) :: m
    Integer                , Intent( In    ) :: n
    Integer                , Intent( In    ) :: mb
    Integer                , Intent( In    ) :: nb
    Integer                , Intent( In    ) :: rsrc
    Integer                , Intent( In    ) :: csrc

    Integer :: nproc, nprow, npcol
    Integer :: error

    map%proc_mapping = proc_map

    map%descriptor = INVALID
    
    map%descriptor( dtype_a ) = 1 ! Dense matrix

    Call mpi_comm_size( map%get_comm(), nproc, error )
    Call factor( nproc, nprow, npcol )
    Call blacs_gridinit( map%get_comm(), 'C', nprow, npcol, map%descriptor( ctxt_a ) )

    map%descriptor( m_a     ) = m ! Global rows
    map%descriptor( n_a     ) = n ! Global cols
    
    map%descriptor( mb_a    ) = mb ! row block fac
    map%descriptor( nb_a    ) = nb ! col block fac

    map%descriptor( rsrc_a  ) = rsrc ! first proc row
    map%descriptor( csrc_a  ) = csrc ! first proc col
    
  End Subroutine set_matrix_mapping

  Subroutine split_matrix_mapping( map, weights, split_name, split_map, i_hold )

    Class( matrix_mapping )                             , Intent( In    ) :: map
    Integer, Dimension( : )                             , Intent( In    ) :: weights
    Character( Len = * )                                , Intent( In    ) :: split_name
    Type ( matrix_mapping ), Dimension( : ), Allocatable, Intent(   Out ) :: split_map
    Integer,                 Dimension( : ), Allocatable, Intent(   Out ) :: i_hold

    Type( proc_mapping ), Dimension( : ), Allocatable :: proc_map

    Integer :: i
    
    Call map%proc_mapping%split( weights, split_name, proc_map, i_hold )

    Allocate( split_map( 1:Size( proc_map ) ) )

    Do i = 1, Size( split_map )
       Call split_map( i )%set( proc_map( i ) ,                 &
            map%descriptor( m_a    ), map%descriptor( n_a    ), &
            map%descriptor( mb_a   ), map%descriptor( nb_a   ), &
            map%descriptor( rsrc_a ), map%descriptor( csrc_a ) )
    End Do

  End Subroutine split_matrix_mapping

  Subroutine factor( a, b, c )

    Integer, Intent( In    ) :: a
    Integer, Intent(   Out ) :: b
    Integer, Intent(   Out ) :: c

    b = Nint( Sqrt( Real( a ) ) )

    Do
       c = a / b
       If( b * c == a ) Then
          Exit
       End If
       b = b - 1 ! Loop will terminate when b = 1
    End Do

  End Subroutine factor

End Module matrix_mapping_module
 
