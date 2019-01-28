Module matrix_mapping_module

  Use numbers_module     , Only : wp 
  Use proc_mapping_module, Only : proc_mapping
  
  Implicit None

  Type, Public  :: matrix_mapping
     Type( proc_mapping )     , Private :: proc_map
     Integer, Dimension( 1:9 ), Private :: descriptor
   Contains
     Procedure :: set => set_matrix_mapping
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

  Subroutine set_matrix_mapping( matrix_map, proc_map, m, n, mb, nb, rsrc, csrc )

    Class( matrix_mapping ), Intent(   Out ) :: matrix_map
    Type ( proc_mapping   ), Intent( In    ) :: proc_map
    Integer                , Intent( In    ) :: m
    Integer                , Intent( In    ) :: n
    Integer                , Intent( In    ) :: mb
    Integer                , Intent( In    ) :: nb
    Integer                , Intent( In    ) :: rsrc
    Integer                , Intent( In    ) :: csrc

    matrix_map%proc_map = proc_map

    matrix_map%descriptor( dtype_a ) = 1 ! Dense matrix

    matrix_map%descriptor( m_a     ) = m ! Global rows
    matrix_map%descriptor( n_a     ) = n ! Global cols
    
    matrix_map%descriptor( mb_a    ) = mb ! row block fac
    matrix_map%descriptor( nb_a    ) = nb ! col block fac

    matrix_map%descriptor( rsrc_a  ) = rsrc ! first proc row
    matrix_map%descriptor( csrc_a  ) = csrc ! first proc col
    
  End Subroutine set_matrix_mapping

End Module matrix_mapping_module
 
