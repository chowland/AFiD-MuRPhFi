!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: HdfRoutines.F90                                !
!    CONTAINS: subroutine hdf_read_serial_1d              !
!                                                         ! 
!    PURPOSE: I/O routines.                               !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      subroutine HdfCreateBlankFile(filename)
      use hdf5
      implicit none
      character*30,intent(in) :: filename
      integer(HID_T) :: file_id
      integer :: hdf_error

      call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, hdf_error)

      call h5fclose_f(file_id, hdf_error)
      
      end subroutine HdfCreateBlankFile

!====================================================================
      subroutine HdfSerialWriteRealScalar(dsetname,filename,n)
      use hdf5
      implicit none
      character*30,intent(in) :: dsetname,filename
      real, intent(in) :: n
      integer(HID_T) :: file_id
      integer(HID_T) :: dset, filespace
      integer :: hdf_error
      integer(HSIZE_T) :: dims(1)
      logical :: fileexists

      call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, hdf_error)

      dims(1)=1

      call h5lexists_f(file_id,dsetname,fileexists,hdf_error)

      if(fileexists) call h5ldelete_f(file_id,dsetname,hdf_error)

      call h5screate_simple_f(1, dims, &
     &                        filespace, hdf_error)

      call h5dcreate_f(file_id, dsetname, H5T_NATIVE_DOUBLE, &
     &                filespace, dset, hdf_error)

       call h5dwrite_f(dset, H5T_NATIVE_DOUBLE, &
     &   n, dims, hdf_error)

      call h5dclose_f(dset, hdf_error)

      call h5sclose_f(filespace, hdf_error)

      call h5fclose_f(file_id, hdf_error)
      
      end subroutine HdfSerialWriteRealScalar

!====================================================================
      subroutine HdfSerialWriteIntScalar(dsetname,filename,n)
      use hdf5
      implicit none
      character*30,intent(in) :: dsetname,filename
      integer, intent(in) :: n
      integer(HID_T) :: file_id
      integer(HID_T) :: dset, filespace
      integer :: hdf_error
      integer(HSIZE_T) :: dims(1)
      logical :: fileexists

      call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, hdf_error)

      dims(1)=1

      call h5lexists_f(file_id,dsetname,fileexists,hdf_error)

      if(fileexists) call h5ldelete_f(file_id,dsetname,hdf_error)

      call h5screate_simple_f(1, dims, &
     &                        filespace, hdf_error)

      call h5dcreate_f(file_id, dsetname, H5T_NATIVE_INTEGER, &
     &                filespace, dset, hdf_error)

       call h5dwrite_f(dset, H5T_NATIVE_INTEGER, &
     &   n, dims, hdf_error)

      call h5dclose_f(dset, hdf_error)

      call h5sclose_f(filespace, hdf_error)

      call h5fclose_f(file_id, hdf_error)
      
      end subroutine HdfSerialWriteIntScalar

!====================================================================
      subroutine HdfSerialWriteInt1D(dsetname,filename,var,sz)
            use hdf5
            implicit none
            character*30,intent(in) :: dsetname,filename
            integer, intent(in) :: sz
            real, dimension(sz), intent(in) :: var
            integer(HID_T) :: file_id
            integer(HID_T) :: dset, filespace
            integer :: hdf_error
            integer(HSIZE_T) :: dims(1)
            logical :: fileexists
      
      
            call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, hdf_error)
      
            dims(1)=sz
      
            call h5screate_simple_f(1, dims, &
           &                        filespace, hdf_error)
      
            call h5lexists_f(file_id,dsetname,fileexists,hdf_error)
      
            if(fileexists) call h5ldelete_f(file_id,dsetname,hdf_error)
      
            call h5dcreate_f(file_id, dsetname, H5T_NATIVE_INTEGER, &
           &                filespace, dset, hdf_error)
      
      
             call h5dwrite_f(dset, H5T_NATIVE_INTEGER, &
           &   var(1:sz), dims, hdf_error)
      
      
            call h5dclose_f(dset, hdf_error)
      
            call h5sclose_f(filespace, hdf_error)
      
            call h5fclose_f(file_id, hdf_error)
            
            end subroutine HdfSerialWriteInt1D
!====================================================================
      subroutine HdfSerialWriteReal1D(dsetname,filename,var,sz)
      use hdf5
      implicit none
      character*30,intent(in) :: dsetname,filename
      integer, intent(in) :: sz
      real, dimension(sz), intent(in) :: var
      integer(HID_T) :: file_id
      integer(HID_T) :: dset, filespace
      integer :: hdf_error
      integer(HSIZE_T) :: dims(1)
      logical :: fileexists


      call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, hdf_error)

      dims(1)=sz

      call h5screate_simple_f(1, dims, &
     &                        filespace, hdf_error)

      call h5lexists_f(file_id,dsetname,fileexists,hdf_error)

      if(fileexists) call h5ldelete_f(file_id,dsetname,hdf_error)

      call h5dcreate_f(file_id, dsetname, H5T_NATIVE_DOUBLE, &
     &                filespace, dset, hdf_error)


       call h5dwrite_f(dset, H5T_NATIVE_DOUBLE, &
     &   var(1:sz), dims, hdf_error)


      call h5dclose_f(dset, hdf_error)

      call h5sclose_f(filespace, hdf_error)

      call h5fclose_f(file_id, hdf_error)
      
      end subroutine HdfSerialWriteReal1D
!====================================================================
      subroutine HdfWriteReal2D(dsetname,filename,var)
      use mpih
      use param
      use hdf5
      use decomp_2d, only: xstart,xend
      implicit none
      character*30,intent(in) :: dsetname,filename
      real, intent(in) :: var(xstart(2):xend(2) &
     &                  ,xstart(3):xend(3))

      integer(HID_T) :: file_id
      integer(HID_T) :: filespace
      integer(HID_T) :: slabspace
      integer(HID_T) :: memspace
      integer(HID_T) :: dset
      integer(HID_T) :: plist_id
      integer(HSIZE_T) :: dims(2)
      integer(HSIZE_T), dimension(2) :: data_count  
      integer(HSSIZE_T), dimension(2) :: data_offset 
      integer :: hdf_error, ndims
      integer :: comm, info

!RO   Sort out MPI definitions

      comm = MPI_COMM_WORLD
      info = MPI_INFO_NULL

      ndims = 2
      dims(1)=nym
      dims(2)=nzm

      call h5screate_simple_f(ndims, dims,  &
     &                        filespace, hdf_error)

      data_count(1) = xend(2)-xstart(2)+1
      data_count(2) = xend(3)-xstart(3)+1

      data_offset(1) = xstart(2)-1
      data_offset(2) = xstart(3)-1

      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, &
     &    hdf_error)
 
      call h5pset_fapl_mpio_f(plist_id, comm, info, &
     &  hdf_error)

      call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, &
     & hdf_error, access_prp=plist_id)

      call h5pclose_f(plist_id, hdf_error)

      call h5dcreate_f(file_id, dsetname, H5T_NATIVE_DOUBLE, &
     &                filespace, dset, hdf_error)
      call h5screate_simple_f(ndims, data_count, memspace, hdf_error) 
      call h5dget_space_f(dset, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
     &                      data_offset, data_count, hdf_error)

      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
     &                        hdf_error)

       call h5dwrite_f(dset, H5T_NATIVE_DOUBLE, &
     &   var(xstart(2):xend(2),xstart(3):xend(3)), dims,  &
     &   hdf_error, file_space_id = slabspace, mem_space_id = memspace,  &
     &   xfer_prp = plist_id)

      call h5pclose_f(plist_id, hdf_error)

      call h5dclose_f(dset, hdf_error)

      call h5sclose_f(memspace, hdf_error)
      call h5fclose_f(file_id, hdf_error)

      
      end subroutine HdfWriteReal2D

!====================================================================
      subroutine HdfSerialReadRealScalar(dsetname,filename,n)
      use hdf5
      implicit none
      character*30,intent(in) :: dsetname,filename
      real, intent(out) :: n
      integer(HID_T) :: file_id
      integer(HID_T) :: dset
      integer :: hdf_error
      integer(HSIZE_T) :: dims(1)


      call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, hdf_error)

      dims(1)=1

      call h5dopen_f(file_id, dsetname, dset, hdf_error)

       call h5dread_f(dset, H5T_NATIVE_DOUBLE, &
     &   n, dims, hdf_error)


      call h5dclose_f(dset, hdf_error)

      call h5fclose_f(file_id, hdf_error)
      
      end subroutine HdfSerialReadRealScalar

!====================================================================
      subroutine HdfSerialReadIntScalar(dsetname,filename,n)
      use hdf5
      implicit none
      character*30,intent(in) :: dsetname,filename
      integer, intent(out) :: n
      integer(HID_T) :: file_id
      integer(HID_T) :: dset
      integer :: hdf_error
      integer(HSIZE_T) :: dims(1)


      call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, hdf_error)

      dims(1)=1

      call h5dopen_f(file_id, dsetname, dset, hdf_error)

      call h5dread_f(dset, H5T_NATIVE_INTEGER, &
     &   n, dims, hdf_error)

      call h5dclose_f(dset, hdf_error)

      call h5fclose_f(file_id, hdf_error)
      
      end subroutine HdfSerialReadIntScalar
!====================================================================
      subroutine HdfSerialReadReal1D(dsetname,filename,var,sz)
      use hdf5
      implicit none
      character*30,intent(in) :: dsetname,filename
      integer, intent(in) :: sz
      real, dimension(sz), intent(out) :: var
      integer(HID_T) :: file_id
      integer(HID_T) :: dset
      integer :: hdf_error
      integer(HSIZE_T) :: dims(1)


      call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, hdf_error)

      dims(1)=sz

      call h5dopen_f(file_id, dsetname, dset, hdf_error)

      call h5dread_f(dset, H5T_NATIVE_DOUBLE, &
     &   var, dims, hdf_error)

      call h5dclose_f(dset, hdf_error)

      call h5fclose_f(file_id, hdf_error)
      
      end subroutine HdfSerialReadReal1D

!====================================================================
      subroutine HdfReadRealHalo3D(filnam1,qua,st2,en2,st3,en3,n1o,n2o,n3o)
      use param
      use mpih
      use hdf5
      
      implicit none

      integer, intent(in) :: n1o,n2o,n3o,st2,en2,st3,en3
      real, intent(out), dimension(1:n3o,st2-lvlhalo:en2+lvlhalo, &
         st3-lvlhalo:en3+lvlhalo)::qua

      integer hdf_error

      integer(HID_T) :: file_id
      integer(HID_T) :: slabspace
      integer(HID_T) :: memspace

      integer(HID_T) :: dset_qua

      integer(HSIZE_T) :: dims(3)

      integer(HID_T) :: plist_id
      integer(HSIZE_T), dimension(3) :: data_count  
      integer(HSSIZE_T), dimension(3) :: data_offset 

      integer :: comm, info
      integer :: ndims

!      real, intent(in), dimension(1:n3,xstart(2)-1:xend(2)+1, &
!     & xstart(3)-1:xend(3)+1)::qua

      character*30,intent(in) :: filnam1

!RO   Sort out MPI definitions

      comm = MPI_COMM_WORLD
      info = MPI_INFO_NULL

!RO   Set offsets and element counts
   
      ndims = 3

      dims(1)=n3o
      dims(2)=n2o-1
      dims(3)=n1o-1

      data_count(1) = n3o
      data_count(2) = en2-st2+1
      data_count(3) = en3-st3+1

      data_offset(1) = 0
      data_offset(2) = st2-1
      data_offset(3) = st3-1


      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, &
          hdf_error)

      call h5pset_fapl_mpio_f(plist_id, comm, info, &
     &  hdf_error)

      call h5fopen_f(filnam1, H5F_ACC_RDONLY_F, file_id, &
     & hdf_error, access_prp=plist_id)
      call h5pclose_f(plist_id, hdf_error)

      call h5dopen_f(file_id, 'var', &
     &                dset_qua, hdf_error)

      call h5screate_simple_f(ndims, data_count, memspace, hdf_error) 

      call h5dget_space_f(dset_qua, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
     &                      data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
     &                        hdf_error)
       call h5dread_f(dset_qua, H5T_NATIVE_DOUBLE, &
     &   qua(1:n3o,st2:en2,st3:en3), dims,  &
     &   hdf_error, file_space_id = slabspace, mem_space_id = memspace,  &
     &   xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)

      call h5dclose_f(dset_qua, hdf_error)

      call h5sclose_f(memspace, hdf_error)
      call h5fclose_f(file_id, hdf_error)

      end subroutine HdfReadRealHalo3D
!====================================================================
      subroutine HdfWriteRealHalo3D(filnam1,qua)
      use param
      use mpih
      use hdf5
      use decomp_2d, only: xstart,xend
      
      implicit none

      integer hdf_error

      integer(HID_T) :: file_id
      integer(HID_T) :: filespace
      integer(HID_T) :: slabspace
      integer(HID_T) :: memspace

      integer(HID_T) :: dset

      integer(HSIZE_T) :: dims(3)

      integer(HID_T) :: plist_id
      integer(HSIZE_T), dimension(3) :: data_count  
      integer(HSSIZE_T), dimension(3) :: data_offset 

      integer :: comm, info
      integer :: ndims

      real, intent(in), dimension(1:nx, &
       xstart(2)-lvlhalo:xend(2)+lvlhalo, &
       xstart(3)-lvlhalo:xend(3)+lvlhalo)::qua

      character*30,intent(in) :: filnam1

!RO   Sort out MPI definitions

      comm = MPI_COMM_WORLD
      info = MPI_INFO_NULL

!RO   Set offsets and element counts
   
      ndims = 3

      dims(1)=nx
      dims(2)=nym
      dims(3)=nzm

      call h5screate_simple_f(ndims, dims,  &
     &                        filespace, hdf_error)

      data_count(1) = nx
      data_count(2) = xend(2)-xstart(2)+1
      data_count(3) = xend(3)-xstart(3)+1

      data_offset(1) = 0
      data_offset(2) = xstart(2)-1
      data_offset(3) = xstart(3)-1

      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, &
          hdf_error)

      call h5pset_fapl_mpio_f(plist_id, comm, info, &
     &  hdf_error)

      call h5fcreate_f(filnam1, H5F_ACC_TRUNC_F, file_id, &
     & hdf_error, access_prp=plist_id)

      call h5pclose_f(plist_id, hdf_error)

      ! call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, hdf_error)
      ! call h5pset_chunk_f(plist_id, ndims, data_count, hdf_error)
      call h5dcreate_f(file_id, 'var', H5T_NATIVE_DOUBLE, &
     &                filespace, dset, hdf_error)
      ! call h5pclose_f(plist_id, hdf_error)

      call h5screate_simple_f(ndims, data_count, memspace, hdf_error) 

      call h5dget_space_f(dset, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
     &                      data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
     &                        hdf_error)
       call h5dwrite_f(dset, H5T_NATIVE_DOUBLE, &
     &   qua(1:nx,xstart(2):xend(2),xstart(3):xend(3)), dims,  &
     &   hdf_error, file_space_id = slabspace, mem_space_id = memspace,  &
     &   xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)

      call h5dclose_f(dset, hdf_error)

      call h5sclose_f(memspace, hdf_error)
      call h5fclose_f(file_id, hdf_error)

      end subroutine HdfWriteRealHalo3D
      
!====================================================================
      subroutine HdfWriteRealHalo3DR(filnam1,qua)
      use param
      use mpih
      use hdf5
      use decomp_2d, only: xstartr,xendr
      
      implicit none

      integer hdf_error

      integer(HID_T) :: file_id
      integer(HID_T) :: filespace
      integer(HID_T) :: slabspace
      integer(HID_T) :: memspace

      integer(HID_T) :: dset

      integer(HSIZE_T) :: dims(3)

      integer(HID_T) :: plist_id
      integer(HSIZE_T), dimension(3) :: data_count  
      integer(HSSIZE_T), dimension(3) :: data_offset 

      integer :: comm, info
      integer :: ndims

      real, intent(in), dimension(1:nxr, &
       xstartr(2)-lvlhalo:xendr(2)+lvlhalo, &
       xstartr(3)-lvlhalo:xendr(3)+lvlhalo)::qua

      character*30,intent(in) :: filnam1

!RO   Sort out MPI definitions

      comm = MPI_COMM_WORLD
      info = MPI_INFO_NULL

!RO   Set offsets and element counts
   
      ndims = 3

      dims(1)=nxr
      dims(2)=nymr
      dims(3)=nzmr

      call h5screate_simple_f(ndims, dims,  &
     &                        filespace, hdf_error)

      data_count(1) = nxr
      data_count(2) = xendr(2)-xstartr(2)+1
      data_count(3) = xendr(3)-xstartr(3)+1

      data_offset(1) = 0
      data_offset(2) = xstartr(2)-1
      data_offset(3) = xstartr(3)-1

      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, &
          hdf_error)

      call h5pset_fapl_mpio_f(plist_id, comm, info, &
     &  hdf_error)

      call h5fcreate_f(filnam1, H5F_ACC_TRUNC_F, file_id, &
     & hdf_error, access_prp=plist_id)

      call h5pclose_f(plist_id, hdf_error)

      call h5dcreate_f(file_id, 'var', H5T_NATIVE_DOUBLE, &
     &                filespace, dset, hdf_error)

      call h5screate_simple_f(ndims, data_count, memspace, hdf_error) 

      call h5dget_space_f(dset, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
     &                      data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
     &                        hdf_error)
       call h5dwrite_f(dset, H5T_NATIVE_DOUBLE, &
     &   qua(1:nxr,xstartr(2):xendr(2),xstartr(3):xendr(3)), dims,  &
     &   hdf_error, file_space_id = slabspace, mem_space_id = memspace,  &
     &   xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)

      call h5dclose_f(dset, hdf_error)

      call h5sclose_f(memspace, hdf_error)
      call h5fclose_f(file_id, hdf_error)

      end subroutine HdfWriteRealHalo3DR
      
!====================================================================

      subroutine HdfWriteSingleHalo3D(filnam1,qua)
      use param
      use mpih
      use hdf5
      use decomp_2d, only: xstart,xend
      
      implicit none

      integer hdf_error

      integer(HID_T) :: file_id
      integer(HID_T) :: filespace
      integer(HID_T) :: slabspace
      integer(HID_T) :: memspace

      integer(HID_T) :: dset

      integer(HSIZE_T) :: dims(3)

      integer(HID_T) :: plist_id
      integer(HSIZE_T), dimension(3) :: data_count  
      integer(HSSIZE_T), dimension(3) :: data_offset 

      integer :: comm, info
      integer :: ndims

      integer, parameter :: sing = selected_real_kind(6, 37)

      real(sing), intent(in), dimension(1:nx, &
       xstart(2)-lvlhalo:xend(2)+lvlhalo, &
       xstart(3)-lvlhalo:xend(3)+lvlhalo)::qua

      character*30,intent(in) :: filnam1

!RO   Sort out MPI definitions

      comm = MPI_COMM_WORLD
      info = MPI_INFO_NULL

!RO   Set offsets and element counts
   
      ndims = 3

      dims(1)=nx !CS For half slab: (nx-2)/2+2
      dims(2)=nym
      dims(3)=nzm

      call h5screate_simple_f(ndims, dims,  &
     &                        filespace, hdf_error)

      data_count(1) = nx !CS For half slab: (nx-2)/2+2
      data_count(2) = xend(2)-xstart(2)+1
      data_count(3) = xend(3)-xstart(3)+1

      data_offset(1) = 0
      data_offset(2) = xstart(2)-1
      data_offset(3) = xstart(3)-1

      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, &
          hdf_error)

      call h5pset_fapl_mpio_f(plist_id, comm, info, &
     &  hdf_error)

      call h5fcreate_f(filnam1, H5F_ACC_TRUNC_F, file_id, &
     & hdf_error, access_prp=plist_id)

      call h5pclose_f(plist_id, hdf_error)

      call h5dcreate_f(file_id, 'var', H5T_NATIVE_REAL, &
     &                filespace, dset, hdf_error)

      call h5screate_simple_f(ndims, data_count, memspace, hdf_error) 

      call h5dget_space_f(dset, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
     &                      data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
     &                        hdf_error)
       call h5dwrite_f(dset, H5T_NATIVE_REAL, &
     &   qua(1:nx,xstart(2):xend(2),xstart(3):xend(3)), dims,  &
     &   hdf_error, file_space_id = slabspace, mem_space_id = memspace,  &
     &   xfer_prp = plist_id)
!CS     &   qua((nx-2)/2+1:nx,xstart(2):xend(2),xstart(3):xend(3)), dims,  &    ! For half slab
      call h5pclose_f(plist_id, hdf_error)

      call h5dclose_f(dset, hdf_error)

      call h5sclose_f(memspace, hdf_error)
      call h5fclose_f(file_id, hdf_error)

      end subroutine HdfWriteSingleHalo3D

!====================================================================

      subroutine HdfStart
      use hdf5
      implicit none
      integer :: hdf_error

      call h5open_f(hdf_error)

      end subroutine HdfStart

!====================================================================
      subroutine HdfClose
      use hdf5
      implicit none
      integer :: hdf_error

      call h5close_f(hdf_error)

      end subroutine HdfClose
