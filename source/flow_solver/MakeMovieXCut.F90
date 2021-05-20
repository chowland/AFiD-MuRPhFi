!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: MakeMovieXCut.F90                              !
!    CONTAINS: subroutine Mkmov_xcut                      !
!                                                         ! 
!    PURPOSE: Write down xcut snapshot for flow vars      !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine Mkmov_xcut
    use param
    use mpih
    use hdf5
    use decomp_2d, only: xstart,xend,xstartr,xendr,DECOMP_2D_COMM_CART_X
    use local_arrays, only: vz,vy,vx,temp
    use mgrd_arrays, only: sal
    implicit none
    character(70) :: filename
    character(30) :: dsetname
    character( 5) :: frame

    integer(HID_T) :: file_id, group_id, plist_id
    integer(HID_T) :: filespace, slabspace, memspace
    integer(HID_T) :: dset_vx, dset_vy, dset_vz, dset_temp, dset_sal
    integer(HSIZE_T), dimension(2) :: dims, data_count, data_offset
    integer :: ic, hdf_error, ndims, comm, info, ierror

    logical :: fexist

    ! Select plane - plane next to lower wall
    ic = 1

    ! Record frame number and filename as strings
    write(frame,"(i5.5)")nint(time/tframe)
    filename = trim("outputdir/flowmov/movie_xcut.h5")

    ! Check if movie file already exists, and create it if not
    inquire(file=filename,exist=fexist)
    if (.not.fexist) then
        if (ismaster) then
        end if
    end if

    ! MPI
    call MPI_CART_SUB(DECOMP_2D_COMM_CART_X, (/.true.,.true./), comm, ierror)
    info = MPI_INFO_NULL

    ! Define dimension
    ndims = 2
    dims(1) = nym
    dims(2) = nzm
    
    data_count(1) = xend(2) - xstart(2) + 1
    data_count(2) = xend(3) - xstart(3) + 1

    data_offset(1) = xstart(2) - 1
    data_offset(2) = xstart(3) - 1

    call h5screate_simple_f(ndims, dims, filespace, hdf_error)

    ! Begin HDF5
    call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf_error)
    call h5pset_fapl_mpio_f(plist_id, comm, info, hdf_error)
    
    ! Check if movie file already exists, and create it if not
    inquire(file=filename,exist=fexist)
    if (fexist) then
        call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, hdf_error, access_prp=plist_id)
    else
        call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, hdf_error,access_prp=plist_id)
        call h5gcreate_f(file_id, "vx", group_id, hdf_error)
        call h5gclose_f(group_id, hdf_error)
        call h5gcreate_f(file_id, "vy", group_id, hdf_error)
        call h5gclose_f(group_id, hdf_error)
        call h5gcreate_f(file_id, "vz", group_id, hdf_error)
        call h5gclose_f(group_id, hdf_error)
        call h5gcreate_f(file_id, "temp", group_id, hdf_error)
        call h5gclose_f(group_id, hdf_error)
        call h5gcreate_f(file_id, "sal", group_id, hdf_error)
        call h5gclose_f(group_id, hdf_error)
    end if

    call h5pclose_f(plist_id, hdf_error)

    ! Create dataspace
    call h5screate_simple_f(ndims, data_count, memspace, hdf_error)

    !! Create datasets and write data to file
    ! vx
    dsetname = trim("vx/"//frame)
    call h5dcreate_f(file_id, dsetname, H5T_NATIVE_DOUBLE, filespace, dset_vx, hdf_error)
    call h5dget_space_f(dset_vx, filespace, hdf_error)
    call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F,data_offset, data_count, hdf_error)
    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, hdf_error)
    call h5dwrite_f(&
            dset_vx, H5T_NATIVE_DOUBLE,&
            vx(ic,xstart(2):xend(2),xstart(3):xend(3)),&
            data_count,  hdf_error, file_space_id = filespace,&
            mem_space_id = memspace, xfer_prp = plist_id)
    call h5pclose_f(plist_id, hdf_error)
    call h5dclose_f(dset_vx, hdf_error)
    
    ! vy
    dsetname = trim("vy/"//frame)
    call h5dcreate_f(file_id, dsetname, H5T_NATIVE_DOUBLE, filespace, dset_vy, hdf_error)
    call h5dget_space_f(dset_vy, filespace, hdf_error)
    call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F,data_offset, data_count, hdf_error)
    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, hdf_error)
    call h5dwrite_f(&
            dset_vy, H5T_NATIVE_DOUBLE,&
            vy(ic,xstart(2):xend(2),xstart(3):xend(3)),&
            data_count,  hdf_error, file_space_id = filespace,&
            mem_space_id = memspace, xfer_prp = plist_id)
    call h5pclose_f(plist_id, hdf_error)
    call h5dclose_f(dset_vy, hdf_error)
    
    ! vz
    dsetname = trim("vz/"//frame)
    call h5dcreate_f(file_id, dsetname, H5T_NATIVE_DOUBLE, filespace, dset_vz, hdf_error)
    call h5dget_space_f(dset_vz, filespace, hdf_error)
    call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F,data_offset, data_count, hdf_error)
    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, hdf_error)
    call h5dwrite_f(&
            dset_vz, H5T_NATIVE_DOUBLE,&
            vz(ic,xstart(2):xend(2),xstart(3):xend(3)),&
            data_count,  hdf_error, file_space_id = filespace,&
            mem_space_id = memspace, xfer_prp = plist_id)
    call h5pclose_f(plist_id, hdf_error)
    call h5dclose_f(dset_vz, hdf_error)
    
    ! temp
    dsetname = trim("temp/"//frame)
    call h5dcreate_f(file_id, dsetname, H5T_NATIVE_DOUBLE, filespace, dset_temp, hdf_error)
    call h5dget_space_f(dset_temp, filespace, hdf_error)
    call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F,data_offset, data_count, hdf_error)
    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, hdf_error)
    call h5dwrite_f(&
            dset_temp, H5T_NATIVE_DOUBLE,&
            temp(ic,xstart(2):xend(2),xstart(3):xend(3)),&
            data_count,  hdf_error, file_space_id = filespace,&
            mem_space_id = memspace, xfer_prp = plist_id)
    call h5pclose_f(plist_id, hdf_error)
    call h5dclose_f(dset_temp, hdf_error)

    call h5sclose_f(memspace, hdf_error)
    call h5fclose_f(file_id, hdf_error)

    !! Repeat on refined grid to save salinity
    ! Select midplane index for refined grid
    ic = 1

    ! Define dimension
    ndims = 2
    dims(1) = nymr
    dims(2) = nzmr
    
    data_count(1) = xendr(2) - xstartr(2) + 1
    data_count(2) = xendr(3) - xstartr(3) + 1

    data_offset(1) = xstartr(2) - 1
    data_offset(2) = xstartr(3) - 1

    call h5screate_simple_f(ndims, dims, filespace, hdf_error)

    ! Begin HDF5
    call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf_error)
    call h5pset_fapl_mpio_f(plist_id, comm, info, hdf_error)
    call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, hdf_error, access_prp=plist_id)
    call h5pclose_f(plist_id, hdf_error)

    ! Create dataspace
    call h5screate_simple_f(ndims, data_count, memspace, hdf_error)

    ! Create and write salinity data
    dsetname = trim("sal/"//frame)
    call h5dcreate_f(file_id, dsetname, H5T_NATIVE_DOUBLE, filespace, dset_sal, hdf_error)
    call h5dget_space_f(dset_sal, filespace, hdf_error)
    call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F,data_offset, data_count, hdf_error)
    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, hdf_error)
    call h5dwrite_f(&
            dset_sal, H5T_NATIVE_DOUBLE,&
            sal(ic,xstartr(2):xendr(2),xstartr(3):xendr(3)),&
            data_count,  hdf_error, file_space_id = filespace,&
            mem_space_id = memspace, xfer_prp = plist_id)
    call h5pclose_f(plist_id, hdf_error)
    call h5dclose_f(dset_sal, hdf_error)
    
    call h5sclose_f(memspace, hdf_error)
    call h5fclose_f(file_id, hdf_error)

end subroutine Mkmov_xcut