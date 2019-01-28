Module matrix_mapping_module

  Use numbers_module     , Only : wp 
  Use proc_mapping_module, Only : proc_mapping
  
  Implicit None

  Type, Extends( proc_mapping ), Public  :: matrix_mapping
     Integer, Dimension( 1:9 ) :: descriptor
  End type matrix_mapping
  
  Private

!!$  Integer

End Module matrix_mapping_module
 
