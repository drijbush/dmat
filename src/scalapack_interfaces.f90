Module scalapack_interfaces
 
  Use numbers_module, Only : wp

  Interface 
 
     Pure Function numroc( n, nb, iproc, isrcproc, nprocs  ) Result( loc_n )
       Implicit None
       Integer               :: loc_n
       Integer, Intent( In ) :: n
       Integer, Intent( In ) :: nb
       Integer, Intent( In ) :: iproc
       Integer, Intent( In ) :: isrcproc
       Integer, Intent( In ) :: nprocs
     End Function numroc

     Subroutine pdsyevd( jobz, uplo, n, a, ia, ja, desca, w, z, iz, jz, descz, work, lwork, iwork, liwork, info )
       Import :: wp
       Implicit None
       Character                    , Intent( In    ) :: jobz
       Character                    , Intent( In    ) :: uplo
       Integer                      , Intent( In    ) :: n
       Real( wp ), Dimension( :, : ), Intent( InOut ) :: a
       Integer                      , Intent( In    ) :: ia
       Integer                      , Intent( In    ) :: ja
       Integer   , Dimension( :    ), Intent( In    ) :: desca
       Real( wp ), Dimension( :    ), Intent(   Out ) :: w
       Real( wp ), Dimension( :, : ), Intent(   Out ) :: z
       Integer                      , Intent( In    ) :: iz
       Integer                      , Intent( In    ) :: jz
       Integer   , Dimension( :    ), Intent( In    ) :: descz
       Real( wp ), Dimension( :    ), Intent(   Out ) :: work
       Integer                      , Intent( In    ) :: lwork
       Integer   , Dimension( :    ), Intent(   Out ) :: iwork
       Integer   ,                    Intent( In    ) :: liwork
       Integer   ,                    Intent( In    ) :: info
     End Subroutine pdsyevd

     Subroutine pzheevd( jobz, uplo, n, a, ia, ja, desca, w, z, iz, jz, descz, cwork, lcwork, rwork, lrwork, iwork, liwork, info )
       Import :: wp
       Implicit None
       Character                       , Intent( In    ) :: jobz
       Character                       , Intent( In    ) :: uplo
       Integer                         , Intent( In    ) :: n
       Complex( wp ), Dimension( :, : ), Intent( InOut ) :: a
       Integer                         , Intent( In    ) :: ia
       Integer                         , Intent( In    ) :: ja
       Integer      , Dimension( :    ), Intent( In    ) :: desca
       Real   ( wp ), Dimension( :    ), Intent(   Out ) :: w
       Complex( wp ), Dimension( :, : ), Intent(   Out ) :: z
       Integer                         , Intent( In    ) :: iz
       Integer                         , Intent( In    ) :: jz
       Integer      , Dimension( :    ), Intent( In    ) :: descz
       Complex( wp ), Dimension( :    ), Intent(   Out ) :: cwork
       Integer                         , Intent( In    ) :: lcwork
       Real( wp )   , Dimension( :    ), Intent(   Out ) :: rwork
       Integer                         , Intent( In    ) :: lrwork
       Integer      , Dimension( :    ), Intent(   Out ) :: iwork
       Integer      ,                    Intent( In    ) :: liwork
       Integer      ,                    Intent( In    ) :: info
     End Subroutine pzheevd
 
  End Interface

End Module scalapack_interfaces
