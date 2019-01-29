Module matrix_mapping_module

  Use numbers_module     , Only : wp
  Use mpi
  Use proc_mapping_module, Only : proc_mapping
  
  Implicit None

  Type, Public  :: matrix_mapping
     Type( proc_mapping )     , Private :: proc_map
     Integer, Dimension( 1:9 ), Private :: descriptor
   Contains
     Procedure :: print => print_matrix_mapping
     Procedure :: set   => set_matrix_mapping
!!$     Procedure :: split => split_matrix_mapping
  End type matrix_mapping
  
  Private

  ! What the elements of the descriptor mean
  Integer, Parameter :: dtype_a = 1 ! descriptor type
  Integer, Parameter :: ctxt_a  = 2 ! context
  Integer, Parameter :: m_a     = 3 ! number of global rows
  Integer, Parameter :: n_a     = 4 ! number of global columns
  Integer, Parameter :: mb_a    = 5 ! row blocking factor
  Integer, Parameter :: nb_a    = 6 ! column blocking factor
  Integer, Parameter :: rsrc_a  = 7 ! first process row which holds a
  Integer, Parameter :: csrc_a  = 8 ! first process col which golds a
  Integer, Parameter :: lld_a   = 9 ! leading dimension of LOCAL a

Contains

  Subroutine print_matrix_mapping( matrix_map )

    Class( matrix_mapping ), Intent( In ) :: matrix_map

    Integer :: rank
    Integer :: error
    
    Call mpi_comm_rank( matrix_map%proc_map%get_comm(), rank, error )
    If( rank == 0 ) Then
       Write( *, '( a, 9( i0, 1x ) )' ) 'Descriptor for matrix mapping'
       Call matrix_map%proc_map%print
    End If
    
  End Subroutine print_matrix_mapping

  Subroutine set_matrix_mapping( matrix_map, proc_map, m, n, mb, nb, rsrc, csrc )

    Class( matrix_mapping ), Intent(   Out ) :: matrix_map
    Type ( proc_mapping   ), Intent( In    ) :: proc_map
    Integer                , Intent( In    ) :: m
    Integer                , Intent( In    ) :: n
    Integer                , Intent( In    ) :: mb
    Integer                , Intent( In    ) :: nb
    Integer                , Intent( In    ) :: rsrc
    Integer                , Intent( In    ) :: csrc

    Integer :: nproc, nprow, npcol
    Integer :: error

    matrix_map%proc_map = proc_map

    matrix_map%descriptor( dtype_a ) = 1 ! Dense matrix

    Call mpi_comm_size( proc_map%get_comm(), nproc, error )
    Call factor( nproc, nprow, npcol )
    Call blacs_gridinit( proc_map%get_comm(), 'C', nprow, npcol, matrix_map%descriptor( ctxt_a ) )

    matrix_map%descriptor( m_a     ) = m ! Global rows
    matrix_map%descriptor( n_a     ) = n ! Global cols
    
    matrix_map%descriptor( mb_a    ) = mb ! row block fac
    matrix_map%descriptor( nb_a    ) = nb ! col block fac

    matrix_map%descriptor( rsrc_a  ) = rsrc ! first proc row
    matrix_map%descriptor( csrc_a  ) = csrc ! first proc col
    
  End Subroutine set_matrix_mapping

!!$  Subroutine split_matrix_mapping( matrix_map, weights, split_name, split_matrix_map, i_hold )
!!$
!!$    Class( matrix_mapping )                             , Intent( In    ) :: matrix_map
!!$    Integer, Dimension( : )                             , Intent( In    ) :: weights
!!$    Character( Len = * )                                , Intent( In    ) :: split_name
!!$    Class( matrix_mapping ), Dimension( : ), Allocatable, Intent(   Out ) :: split_matrix_map
!!$    Integer,                 Dimension( : ), Allocatable, Intent(   Out ) :: i_hold
!!$
!!$    Call matrix_map%proc_map%split( weights, split_name, split_matrix_map, i_hold )
!!$
!!$  End Subroutine split_matrix_mapping

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
 
