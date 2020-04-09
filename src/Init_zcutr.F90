!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: Init_zcutr.F90                                 !
!    CONTAINS: subroutine Initgrid for zcut               !
!                                                         ! 
!    PURPOSE: Initgrid zcut                               !
!     cordin_zcutr_info.h5                                !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine Init_zcutr
      use mpih
      use param
      use hdf5

      IMPLICIT none

      character*70 namfile
      character*30 :: dsetname

      real xmov(nxr,nymr),ymov(nxr,nymr)
      integer i,j

      integer hdf_error,ndims
      integer(HID_T) :: file_id
      integer(HID_T) :: cordin_dset_rid
      integer(HID_T) :: cordin_dset_zid
      integer(HID_T) :: cordin_dspace_id
      integer(HSIZE_T) :: dims(2)

      do j=1,nymr
        do i=1,nxr
          xmov(i,j) = xcr(i)
          ymov(i,j) = ymr(j)
        end do
      end do

      ndims = 2
      dims(1)=nxr
      dims(2)=nymr
      if (ismaster) then 
       namfile='./outputdir/flowmov/cordin_zcutr_info.h5'
       call h5fcreate_f(namfile, H5F_ACC_TRUNC_F, file_id, hdf_error)
       call h5screate_simple_f(ndims, dims, cordin_dspace_id, hdf_error)
       call h5dcreate_f(file_id, 'x', H5T_NATIVE_DOUBLE,cordin_dspace_id,cordin_dset_rid, hdf_error)
       call h5dwrite_f(cordin_dset_rid, H5T_NATIVE_DOUBLE,xmov, dims, hdf_error)
       call h5dcreate_f(file_id, 'y', H5T_NATIVE_DOUBLE, cordin_dspace_id,cordin_dset_zid, hdf_error)
       call h5dwrite_f(cordin_dset_zid, H5T_NATIVE_DOUBLE,ymov, dims, hdf_error)
       call h5dclose_f(cordin_dset_rid, hdf_error)
       call h5dclose_f(cordin_dset_zid, hdf_error)
       call h5sclose_f(cordin_dspace_id, hdf_error)
       call h5fclose_f(file_id, hdf_error)

      endif

      return
      end


