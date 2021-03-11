!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: Mkmov_xcut.F90                                 !
!    CONTAINS: subroutine Mkmov_xcut                      !
!                                                         ! 
!    PURPOSE: Write down xcut snapshot for flow vars      !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine Mkmov_xcut
    use param
    use mpih
    use hdf5
    use decomp_2d, only: xstart,xend,DECOMP_2D_COMM_CART_X
    use local_arrays, only: vz,vy,vx,temp,pr
    implicit none
    character(70) namfile
    character(30) :: dsetname
    character(5) :: ipfi

    real tprfi
    integer itime
    integer ic

    integer(HID_T) :: file_id
    integer(HID_T) :: filespace
    integer(HID_T) :: slabspace
    integer(HID_T) :: memspace,tspace,aid
    integer(HID_T) :: dset_vx,dset_vy,dset_vz,dset_temp
    integer(HID_T) :: plist_id
    integer(HSIZE_T) :: dims(2)
    integer(HSIZE_T), dimension(2) :: data_count  
    integer(HSSIZE_T), dimension(2) :: data_offset
    integer :: hdf_error, ndims
    integer :: comm, info, ierror

    ! File name
    tprfi = 1.d0/tframe
    itime = nint(time*tprfi)
    write(ipfi,"(i5.5)")itime

    namfile="outputdir/flowmov/frame_"//ipfi//"_xcut.h5"

    ! Select plane - plane against lower wall
    ic = 1

    ! MPI
    call MPI_CART_SUB(DECOMP_2D_COMM_CART_X,(/.true.,.true./),comm,ierror)
    info = MPI_INFO_NULL

    ! Define dimension
    ndims = 2
    dims(1) = nym
    dims(2) = nzm

    data_count(1) = xend(2)-xstart(2)+1
    data_count(2) = xend(3)-xstart(3)+1

    data_offset(1) = xstart(2)-1
    data_offset(2) = xstart(3)-1

    call h5screate_simple_f(ndims, dims, filespace, hdf_error)

    ! Begin HDF5
    call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf_error)
    call h5pset_fapl_mpio_f(plist_id, comm, info, hdf_error)
    call h5fcreate_f(namfile, H5F_ACC_TRUNC_F, file_id, hdf_error, access_prp=plist_id)
    call h5pclose_f(plist_id, hdf_error)

    ! Write time as attribute to file
    if (ismaster) then
        call h5screate_f(H5S_SCALAR_F,tspace,hdf_error)
        call h5acreate_f(file_id,"Time",H5T_IEEE_F64LE,tspace,aid,hdf_error)
        call h5awrite_f(aid, H5T_NATIVE_DOUBLE, time, 1, hdf_error)
        call h5aclose_f(aid, hdf_error)
        call h5sclose_f(tspace, hdf_error)
    end if

    call h5dcreate_f(file_id, "vx", H5T_NATIVE_DOUBLE, filespace, dset_vx, hdf_error)
    call h5dcreate_f(file_id, "vy", H5T_NATIVE_DOUBLE, filespace, dset_vy, hdf_error)
    call h5dcreate_f(file_id, "vz", H5T_NATIVE_DOUBLE, filespace, dset_vz, hdf_error)
    call h5dcreate_f(file_id, "temp", H5T_NATIVE_DOUBLE, filespace, dset_temp, hdf_error)

    call h5screate_simple_f(ndims, data_count, memspace, hdf_error)

    ! vx
    call h5dget_space_f(dset_vx, filespace, hdf_error)
    call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, data_offset, data_count, hdf_error)
    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error)
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, hdf_error)
    call h5dwrite_f(&
            dset_vx, H5T_NATIVE_DOUBLE,&
            vx(ic,xstart(2):xend(2),xstart(3):xend(3)),&
            data_count, hdf_error, file_space_id = filespace,&
            mem_space_id = memspace, xfer_prp = plist_id&
            )
    call h5pclose_f(plist_id, hdf_error)

    ! vy
    call h5dget_space_f(dset_vy, filespace, hdf_error)
    call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, data_offset, data_count, hdf_error)
    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error)
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, hdf_error)
    call h5dwrite_f(&
            dset_vy, H5T_NATIVE_DOUBLE,&
            vy(ic,xstart(2):xend(2),xstart(3):xend(3)),&
            data_count, hdf_error, file_space_id = filespace,&
            mem_space_id = memspace, xfer_prp = plist_id&
            )
    call h5pclose_f(plist_id, hdf_error)

    ! vz
    call h5dget_space_f(dset_vz, filespace, hdf_error)
    call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, data_offset, data_count, hdf_error)
    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error)
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, hdf_error)
    call h5dwrite_f(&
            dset_vz, H5T_NATIVE_DOUBLE,&
            vz(ic,xstart(2):xend(2),xstart(3):xend(3)),&
            data_count, hdf_error, file_space_id = filespace,&
            mem_space_id = memspace, xfer_prp = plist_id&
            )
    call h5pclose_f(plist_id, hdf_error)

    ! temp
    call h5dget_space_f(dset_temp, filespace, hdf_error)
    call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, data_offset, data_count, hdf_error)
    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error)
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, hdf_error)
    call h5dwrite_f(&
            dset_temp, H5T_NATIVE_DOUBLE,&
            temp(ic,xstart(2):xend(2),xstart(3):xend(3)),&
            data_count, hdf_error, file_space_id = filespace,&
            mem_space_id = memspace, xfer_prp = plist_id&
            )
    call h5pclose_f(plist_id, hdf_error)

    call h5dclose_f(dset_vx, hdf_error)
    call h5dclose_f(dset_vy, hdf_error)
    call h5dclose_f(dset_vz, hdf_error)
    call h5dclose_f(dset_temp, hdf_error)
    call h5sclose_f(memspace, hdf_error)
    call h5fclose_f(file_id, hdf_error)

end subroutine Mkmov_xcut

!=======================================================

subroutine Mkmov_xcutr
    use param
    use mpih
    use hdf5
    use decomp_2d, only: xstartr,xendr,DECOMP_2D_COMM_CART_X
    ! use mgrd_arrays, only: vzr,vyr,vxr ! ,temp,pr
    use local_arrays, only: sal
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
    integer(HID_T) :: memspace,tspace,aid
    integer(HID_T) :: dset_sal!,dset_vx,dset_vy,dset_vz
    integer(HID_T) :: plist_id
    integer(HSIZE_T) :: dims(2)
    integer(HSIZE_T), dimension(2) :: data_count  
    integer(HSSIZE_T), dimension(2) :: data_offset 
    integer :: hdf_error, ndims
    integer :: comm, info, ierror

    ! File name
    tprfi = 1.d0/tframe
    itime = nint(time*tprfi)
    write(ipfi,"(i5.5)")itime

    namfile="outputdir/flowmov/frame_"//ipfi//"_xcutr.h5"

    ! Select plane - plane against lower wall
    ic = 1

    ! MPI
    call MPI_CART_SUB(DECOMP_2D_COMM_CART_X,(/.true.,.true./),comm,ierror)
    info = MPI_INFO_NULL

    ! Define dimension
    ndims = 2
    dims(1) = nymr
    dims(2) = nzmr

    data_count(1) = xendr(2)-xstartr(2)+1
    data_count(2) = xendr(3)-xstartr(3)+1

    data_offset(1) = xstartr(2)-1
    data_offset(2) = xstartr(3)-1

    call h5screate_simple_f(ndims, dims, filespace, hdf_error)

    ! Begin HDF5
    call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf_error)
    call h5pset_fapl_mpio_f(plist_id, comm, info, hdf_error)
    call h5fcreate_f(namfile, H5F_ACC_TRUNC_F, file_id, hdf_error, access_prp=plist_id)
    call h5pclose_f(plist_id, hdf_error)

    ! Write time as attribute to file
    if (ismaster) then
        call h5screate_f(H5S_SCALAR_F,tspace,hdf_error)
        call h5acreate_f(file_id,"Time",H5T_IEEE_F64LE,tspace,aid,hdf_error)
        call h5awrite_f(aid, H5T_NATIVE_DOUBLE, time, 1, hdf_error)
        call h5aclose_f(aid, hdf_error)
        call h5sclose_f(tspace, hdf_error)
    end if

    ! call h5dcreate_f(file_id, "vxr", H5T_NATIVE_DOUBLE, filespace, dset_vx, hdf_error)
    ! call h5dcreate_f(file_id, "vyr", H5T_NATIVE_DOUBLE, filespace, dset_vy, hdf_error)
    ! call h5dcreate_f(file_id, "vzr", H5T_NATIVE_DOUBLE, filespace, dset_vz, hdf_error)
    call h5dcreate_f(file_id, "sal", H5T_NATIVE_DOUBLE, filespace, dset_sal, hdf_error)

    call h5screate_simple_f(ndims, data_count, memspace, hdf_error)

    ! vx
    ! call h5dget_space_f(dset_vx, filespace, hdf_error)
    ! call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, data_offset, data_count, hdf_error)
    ! call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error)
    ! call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, hdf_error)
    ! call h5dwrite_f(&
    !         dset_vx, H5T_NATIVE_DOUBLE,&
    !         vxr(ic,xstartr(2):xendr(2),xstartr(3):xendr(3)),&
    !         data_count, hdf_error, file_space_id = filespace,&
    !         mem_space_id = memspace, xfer_prp = plist_id&
    !         )
    ! call h5pclose_f(plist_id, hdf_error)

    ! vy
    ! call h5dget_space_f(dset_vy, filespace, hdf_error)
    ! call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, data_offset, data_count, hdf_error)
    ! call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error)
    ! call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, hdf_error)
    ! call h5dwrite_f(&
    !         dset_vy, H5T_NATIVE_DOUBLE,&
    !         vyr(ic,xstartr(2):xendr(2),xstartr(3):xendr(3)),&
    !         data_count, hdf_error, file_space_id = filespace,&
    !         mem_space_id = memspace, xfer_prp = plist_id&
    !         )
    ! call h5pclose_f(plist_id, hdf_error)

    ! vz
    ! call h5dget_space_f(dset_vz, filespace, hdf_error)
    ! call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, data_offset, data_count, hdf_error)
    ! call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error)
    ! call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, hdf_error)
    ! call h5dwrite_f(&
    !         dset_vz, H5T_NATIVE_DOUBLE,&
    !         vzr(ic,xstartr(2):xendr(2),xstartr(3):xendr(3)),&
    !         data_count, hdf_error, file_space_id = filespace,&
    !         mem_space_id = memspace, xfer_prp = plist_id&
    !         )
    ! call h5pclose_f(plist_id, hdf_error)

    ! sal
    call h5dget_space_f(dset_sal, filespace, hdf_error)
    call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, data_offset, data_count, hdf_error)
    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error)
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, hdf_error)
    call h5dwrite_f(&
            dset_sal, H5T_NATIVE_DOUBLE,&
            sal(ic,xstartr(2):xendr(2),xstartr(3):xendr(3)),&
            data_count, hdf_error, file_space_id = filespace,&
            mem_space_id = memspace, xfer_prp = plist_id&
            )
    call h5pclose_f(plist_id, hdf_error)

    ! call h5dclose_f(dset_vx, hdf_error)
    ! call h5dclose_f(dset_vy, hdf_error)
    ! call h5dclose_f(dset_vz, hdf_error)
    call h5dclose_f(dset_sal, hdf_error)
    call h5sclose_f(memspace, hdf_error)
    call h5fclose_f(file_id, hdf_error)

end subroutine Mkmov_xcutr
