!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: Mkmov_zcutr.F90                                !
!    CONTAINS: subroutine Mkmov_zcut                      !
!                                                         ! 
!    PURPOSE: Write down zcut snapshot for                !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine Mkmov_zcutr
      use param
      use mpih
      use hdf5
      use decomp_2d, only: xstartr,xendr,DECOMP_2D_COMM_CART_X
      use mgrd_arrays, only: vzr,vyr,vxr!,temp,pr
      implicit none
      character(70) namfile,xdmnam
      character(30) :: dsetname
      character(5) :: ipfi

      real tprfi
      integer itime
      integer ic

      integer(HID_T) :: file_id
      integer(HID_T) :: filespace
      integer(HID_T) :: slabspace
      integer(HID_T) :: memspace
      integer(HID_T) :: dset_vx,dset_vy,dset_vz
      integer(HID_T) :: plist_id
      integer(HSIZE_T) :: dims(2)
      integer(HSIZE_T), dimension(2) :: data_count  
      integer(HSSIZE_T), dimension(2) :: data_offset 
      integer :: hdf_error, ndims
      integer :: comm, info, ierror


      !-- file name      
      tprfi = 1.d0/tframe
      itime=nint(time*tprfi)
      write(ipfi,'(i5.5)')itime

      namfile='outputdir/flowmov/frame_'//ipfi//'_zcutr.h5'
      xdmnam='outputdir/flowmov/frame_'//ipfi//'_zcutr.xmf'

      !-- select mid plane
      ic=nzmr/2+1

      !-- MPI
      call MPI_CART_SUB(DECOMP_2D_COMM_CART_X,(/.true.,.false./),comm,ierror)
      info = MPI_INFO_NULL

      if(ic.le.xendr(3) .and. ic.ge.xstartr(3)) then

         !-- Define dimension
         ndims = 2
         dims(1)=nxr
         dims(2)=nymr

         data_count(1) = nxr
         data_count(2) = xendr(2)-xstartr(2)+1

         data_offset(1) = 0
         data_offset(2) = xstartr(2)-1


         call h5screate_simple_f(ndims, dims,  filespace, hdf_error)

         !-- Begin HDF5
         call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf_error)
         call h5pset_fapl_mpio_f(plist_id, comm, info, hdf_error)
         call h5fcreate_f(namfile, H5F_ACC_TRUNC_F, file_id, hdf_error,access_prp=plist_id)
         call h5pclose_f(plist_id, hdf_error)

         call h5dcreate_f(file_id, 'vxr', H5T_NATIVE_DOUBLE,filespace,dset_vx, hdf_error)
         !call h5dcreate_f(file_id, 'vyr', H5T_NATIVE_DOUBLE,filespace,dset_vy, hdf_error)
         !call h5dcreate_f(file_id, 'vzr', H5T_NATIVE_DOUBLE,filespace,dset_vz, hdf_error)
         ! call h5dcreate_f(file_id, 'fi', H5T_NATIVE_DOUBLE,filespace,dset_fi, hdf_error)

         call h5screate_simple_f(ndims, data_count, memspace, hdf_error) 

            !-- vx
            call h5dget_space_f(dset_vx, filespace, hdf_error)
            call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F,data_offset, data_count, hdf_error)
            call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
            call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, hdf_error)
            call h5dwrite_f(dset_vx, H5T_NATIVE_DOUBLE,vxr(1:nxr,xstartr(2):xendr(2),ic), dims,  hdf_error, file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)
            call h5pclose_f(plist_id, hdf_error)

            !!-- vy
            !call h5dget_space_f(dset_vy, filespace, hdf_error)
            !call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F,data_offset, data_count, hdf_error)
            !call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
            !call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, hdf_error)
            !call h5dwrite_f(dset_vy, H5T_NATIVE_DOUBLE,vyr(1:nxr,xstartr(2):xendr(2),ic), dims,  hdf_error, file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)
            !call h5pclose_f(plist_id, hdf_error)

            !!-- vz
            !call h5dget_space_f(dset_vz, filespace, hdf_error)
            !call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F,data_offset, data_count, hdf_error)
            !call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
            !call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, hdf_error)
            !call h5dwrite_f(dset_vz, H5T_NATIVE_DOUBLE,vzr(1:nxr,xstartr(2):xendr(2),ic), dims,  hdf_error, file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)
            !call h5pclose_f(plist_id, hdf_error)

            !!-- fi
            !call h5dget_space_f(dset_fi, filespace, hdf_error)
            !call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F,data_offset, data_count, hdf_error)
            !call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
            !call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, hdf_error)
            !call h5dwrite_f(dset_fi, H5T_NATIVE_DOUBLE,fi(1:nx,xstart(2):xend(2),ic), dims,  hdf_error, file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)
            !call h5pclose_f(plist_id, hdf_error)

         call h5dclose_f(dset_vx, hdf_error)
         !call h5dclose_f(dset_vy, hdf_error)
         !call h5dclose_f(dset_vz, hdf_error)
         !call h5dclose_f(dset_fi, hdf_error)
         call h5sclose_f(memspace, hdf_error)
         call h5fclose_f(file_id, hdf_error)

      endif


      !-- Write XMF
      if (ismaster) then 
      open(45,file=xdmnam,status='unknown')
      rewind(45)
      write(45,'("<?xml version=""1.0"" ?>")')
      write(45,'("<!DOCTYPE Xdmf SYSTEM ""Xdmf.dtd"" []>")')
      write(45,'("<Xdmf Version=""2.0"">")')
      write(45,'("<Domain>")')
      write(45,'("<Grid Name=""zcutr"" GridType=""Uniform"">")')

      write(45,'("<Topology TopologyType=""2DSMesh"" NumberOfElements=""",i4," ",i4,"""/>")') nym,nx

      write(45,'("<Geometry GeometryType=""X_Y"">")')
      write(45,'("<DataItem Dimensions=""",i4," ",i4,""" NumberType=""Float"" Precision=""4"" Format=""HDF"">")') nym,nx
      write(45,'("cordin_zcutr_info.h5:/x")')
      write(45,'("</DataItem>")')
      write(45,'("<DataItem Dimensions=""",i4," ",i4,""" NumberType=""Float"" Precision=""4"" Format=""HDF"">")') nym,nx
      write(45,'("cordin_zcutr_info.h5:/y")')
      write(45,'("</DataItem>")')
      write(45,'("</Geometry>")')

      write(45,'("<Attribute Name=""vxr"" AttributeType=""Scalar"" Center=""Node"">")')
      write(45,'("<DataItem Dimensions=""",i4," ",i4,""" NumberType=""Float"" Precision=""4"" Format=""HDF"">")')nym,nx
      write(45,'("frame_",i5.5,"_zcutr.h5:/vxr")') itime
      write(45,'("</DataItem>")')
      write(45,'("</Attribute>")')

      !write(45,'("<Attribute Name=""vyr"" AttributeType=""Scalar"" Center=""Node"">")')
      !write(45,'("<DataItem Dimensions=""",i4," ",i4,""" NumberType=""Float"" Precision=""4"" Format=""HDF"">")')nym,nx
      !write(45,'("frame_",i5.5,"_zcutr.h5:/vyr")') itime
      !write(45,'("</DataItem>")')
      !write(45,'("</Attribute>")')

      !write(45,'("<Attribute Name=""vzr"" AttributeType=""Scalar"" Center=""Node"">")')
      !write(45,'("<DataItem Dimensions=""",i4," ",i4,""" NumberType=""Float"" Precision=""4"" Format=""HDF"">")')nym,nx
      !write(45,'("frame_",i5.5,"_zcutr.h5:/vzr")') itime
      !write(45,'("</DataItem>")')
      !write(45,'("</Attribute>")')

      write(45,'("<Time Value=""",e12.5""" />")')time

      write(45,'("</Grid>")')
      write(45,'("</Domain>")')
      write(45,'("</Xdmf>")')
      close(45) 
      endif

      end subroutine Mkmov_zcutr
