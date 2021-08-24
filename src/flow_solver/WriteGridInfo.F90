!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: WriteGridInfo.F90                              !
!    CONTAINS: subroutine WriteGridInfo                   !
!                                                         ! 
!    PURPOSE: Write the grid information in               !
!     cordin_info.h5                                      !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine WriteGridInfo
      use mpih
      use param
      use hdf5

      IMPLICIT none

      character*70 namfile
      character*30 :: dsetname


      if (ismaster) then 
       namfile='outputdir/cordin_info.h5'
       call HdfCreateBlankFile(namfile)

       dsetname = trim('xm')
       call HdfSerialWriteReal1D(dsetname,namfile,xm,nxm)
       dsetname = trim('xc')
       call HdfSerialWriteReal1D(dsetname,namfile,xc,nx)
       dsetname = trim('ym')
       call HdfSerialWriteReal1D(dsetname,namfile,ym(1:nym),nym)
       dsetname = trim('zm')
       call HdfSerialWriteReal1D(dsetname,namfile,zm(1:nzm),nzm)

       if (multires) then
        dsetname = trim('xmr')
        call HdfSerialWriteReal1D(dsetname,namfile,xmr(1:nxmr),nxmr)
        dsetname = trim('xcr')
        call HdfSerialWriteReal1D(dsetname,namfile,xcr,nxr)
        dsetname = trim('ymr')
        call HdfSerialWriteReal1D(dsetname,namfile,ymr(1:nymr),nymr)
        dsetname = trim('zmr')
        call HdfSerialWriteReal1D(dsetname,namfile,zmr(1:nzmr),nzmr)
       end if

      endif

      return
      end


