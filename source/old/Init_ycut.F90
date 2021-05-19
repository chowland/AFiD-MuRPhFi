!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: Init_ycut.F90                                  !
!    CONTAINS: subroutine Initgrid for ycut               !
!                                                         ! 
!    PURPOSE: Initgrid ycut                               !
!     cordin_ycut_info.h5                                 !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine Init_ycut
      use mpih
      use param
      use hdf5

      IMPLICIT none

      character*70 namfile
      character*30 :: dsetname

      real xmov(nx,nzm),zmov(nx,nzm)
      integer i,j

      integer hdf_error,ndims
      integer(HID_T) :: file_id
      integer(HID_T) :: cordin_dset_rid
      integer(HID_T) :: cordin_dset_zid
      integer(HID_T) :: cordin_dspace_id
      integer(HSIZE_T) :: dims(2)

      do j=1,nzm
        do i=1,nx
          xmov(i,j) = xc(i)
          zmov(i,j) = zm(j)
        end do
      end do

      ndims = 2
      dims(1)=nx
      dims(2)=nzm
      if (ismaster) then 
       namfile='./outputdir/flowmov/cordin_ycut_info.h5'
       call h5fcreate_f(namfile, H5F_ACC_TRUNC_F, file_id, hdf_error)
       call h5screate_simple_f(ndims, dims, cordin_dspace_id, hdf_error)
       call h5dcreate_f(file_id, 'x', H5T_NATIVE_DOUBLE, cordin_dspace_id, cordin_dset_rid, hdf_error)
       call h5dwrite_f(cordin_dset_rid, H5T_NATIVE_DOUBLE, xmov, dims, hdf_error)
       call h5dcreate_f(file_id, 'z', H5T_NATIVE_DOUBLE, cordin_dspace_id, cordin_dset_zid, hdf_error)
       call h5dwrite_f(cordin_dset_zid, H5T_NATIVE_DOUBLE, zmov, dims, hdf_error)
       call h5dclose_f(cordin_dset_rid, hdf_error)
       call h5dclose_f(cordin_dset_zid, hdf_error)
       call h5sclose_f(cordin_dspace_id, hdf_error)
       call h5fclose_f(file_id, hdf_error)

      endif

      return
      end

      subroutine Init_ycutr
      use mpih
      use param
      use hdf5

      IMPLICIT none

      character*70 namfile
      character*30 :: dsetname

      real xmov(nxr,nzmr),zmov(nxr,nzmr)
      integer i,j

      integer hdf_error,ndims
      integer(HID_T) :: file_id
      integer(HID_T) :: cordin_dset_rid
      integer(HID_T) :: cordin_dset_zid
      integer(HID_T) :: cordin_dspace_id
      integer(HSIZE_T) :: dims(2)

      do j=1,nzmr
        do i=1,nxr
          xmov(i,j) = xcr(i)
          zmov(i,j) = zmr(j)
        end do
      end do

      ndims = 2
      dims(1)=nxr
      dims(2)=nzmr
      if (ismaster) then 
       namfile='./outputdir/flowmov/cordin_ycutr_info.h5'
       call h5fcreate_f(namfile, H5F_ACC_TRUNC_F, file_id, hdf_error)
       call h5screate_simple_f(ndims, dims, cordin_dspace_id, hdf_error)
       call h5dcreate_f(file_id, 'x', H5T_NATIVE_DOUBLE, cordin_dspace_id, cordin_dset_rid, hdf_error)
       call h5dwrite_f(cordin_dset_rid, H5T_NATIVE_DOUBLE, xmov, dims, hdf_error)
       call h5dcreate_f(file_id, 'z', H5T_NATIVE_DOUBLE, cordin_dspace_id, cordin_dset_zid, hdf_error)
       call h5dwrite_f(cordin_dset_zid, H5T_NATIVE_DOUBLE, zmov, dims, hdf_error)
       call h5dclose_f(cordin_dset_rid, hdf_error)
       call h5dclose_f(cordin_dset_zid, hdf_error)
       call h5sclose_f(cordin_dspace_id, hdf_error)
       call h5fclose_f(file_id, hdf_error)

      endif

      return
      end

