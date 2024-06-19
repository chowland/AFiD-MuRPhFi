module h5_tools
    use HDF5
    use mpih, only: MPI_INFO_NULL, MPI_COMM_WORLD
    use param, only: salinity, phasefield, nx, nxm, nym, nzm, time, tframe, nxmr, nymr, nzmr, IBM, moist,multiRes_Temp
    use decomp_2d, only: xstart, xstartr, xend
    implicit none

    integer(HID_T) :: group_id, plist_id
    integer :: info, hdf_error

contains

subroutine h5_open_or_create(file_id, filename, comm, fexist)
    integer(HID_T), intent(out) :: file_id
    integer, intent(in) :: comm
    character(70), intent(in) :: filename
    logical, intent(out) :: fexist

    info = MPI_INFO_NULL

    call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf_error)
    call h5pset_fapl_mpio_f(plist_id, comm, info, hdf_error)

    inquire(file=filename, exist=fexist)
    if (fexist) then
        call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, hdf_error, access_prp=plist_id)
    else
        call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, hdf_error, access_prp=plist_id)
    end if

    call h5pclose_f(plist_id, hdf_error)
end subroutine

subroutine h5_add_slice_groups(file_id)
    use param, only:multiRes_Temp
    integer(HID_T), intent(inout) :: file_id
     
    integer :: i
    character(len=15), dimension(4) :: var_list
    

   
    var_list = [character(len=15) :: "vx", "vy", "vz", "temp"]
    do i=1,size(var_list)
        call h5gcreate_f(file_id, trim(var_list(i)), group_id, hdf_error)
        call h5gclose_f(group_id, hdf_error)
    end do


    if (salinity) then
        call h5gcreate_f(file_id, "sal", group_id, hdf_error)
        call h5gclose_f(group_id, hdf_error)
        if (IBM) then
            call h5gcreate_f(file_id, "sal2", group_id, hdf_error)
            call h5gclose_f(group_id, hdf_error)
        end if
    end if

    if (phasefield) then
        call h5gcreate_f(file_id, "phi", group_id, hdf_error)
        call h5gclose_f(group_id, hdf_error)
    end if

    if (moist) then
        call h5gcreate_f(file_id, "qhum", group_id, hdf_error)
        call h5gclose_f(group_id, hdf_error)
    end if
end subroutine

subroutine write_H5_plane(file_id, varname, var, axis)
    use param,only:multiRes_Temp
    integer(HID_T), intent(inout) :: file_id
    character(len=15), intent(in) :: varname
    real, intent(in) :: var(:,:)
    character, intent(in) :: axis

    integer, parameter :: ndims = 2
    integer(HID_T) :: filespace, memspace, dset
    integer(HSIZE_T), dimension(2) :: dims, data_count, data_offset
    character(len=5) :: frame
    character(len=20) :: dsetname
    logical :: dexist
    if (axis=='x') then
        dims = [nym, nzm]
        data_offset = [xstart(2) - 1, xstart(3) - 1]
    else
        if (axis=='y') then
            dims(2) = nzm
            data_offset = [0, xstart(3) - 1]
        else
            dims(2) = nym
            data_offset = [0, xstart(2) - 1]
        end if
        if (trim(varname)=='vx') then
            dims(1) = nx
        else
            dims(1) = nxm
        end if
    end if
    if (index('phisal2', trim(varname)) /= 0 .or. (varname== 'temp'.and. multiRes_Temp)) then
        dims = [nxmr, nzmr]
        data_offset = [0, xstartr(3) - 1]
        if (axis=='x') then
            dims(1) = nymr
            data_offset(1) = xstartr(2) - 1
        elseif (axis=='z') then
            dims(2) = nymr
            data_offset(2) = xstartr(2) - 1
        end if
    end if

    data_count = [size(var, dim=1), size(var, dim=2)]

    ! Create dataspaces
    call h5screate_simple_f(ndims, dims, filespace, hdf_error)
    call h5screate_simple_f(ndims, data_count, memspace, hdf_error)


    write(frame,"(i5.5)")nint(time/tframe)
    dsetname = trim(trim(varname)//"/"//frame)
    call h5lexists_f(file_id, dsetname, dexist, hdf_error)
    if (dexist) call h5ldelete_f(file_id, dsetname, hdf_error)
    call h5dcreate_f(file_id, dsetname, H5T_NATIVE_DOUBLE, filespace, dset, hdf_error)
    call h5dget_space_f(dset, filespace, hdf_error)

    call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F,data_offset, data_count, hdf_error)
    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, hdf_error)
    call h5dwrite_f(&
        dset, H5T_NATIVE_DOUBLE,&
        var,&
        data_count,  hdf_error, file_space_id = filespace,&
        mem_space_id = memspace, xfer_prp = plist_id)
    call h5pclose_f(plist_id, hdf_error)
    call h5dclose_f(dset, hdf_error)
    call h5sclose_f(filespace, hdf_error)
    call h5sclose_f(memspace, hdf_error)

end subroutine

!> Initialise the MPI communicators used for writing 2D slices
subroutine InitSliceCommunicators
    use mpih, only: comm_xy, comm_xz, comm_yz
    use decomp_2d, only: DECOMP_2D_COMM_CART_X
    integer :: ierr

    !! yz-plane (cut of constant x)
    call MPI_CART_SUB(DECOMP_2D_COMM_CART_X, (/.true., .true./), comm_yz, ierr)
    !! xy-plane (cut of constant z)
    call MPI_CART_SUB(DECOMP_2D_COMM_CART_X, (/.true.,.false./), comm_xy, ierr)
    !! xz-plane (cut of constant y)
    call MPI_CART_SUB(DECOMP_2D_COMM_CART_X, (/.false.,.true./), comm_xz, ierr)

end subroutine InitSliceCommunicators

!> Write out a 3D array that does not have a halo
subroutine write_3D_array(filename, arr)
    character(len=30), intent(in) :: filename
    real, dimension(1:nx,xstart(2):xend(2),xstart(3):xend(3)), intent(in) :: arr

    integer(hid_t) :: file_id, filespace, slabspace, memspace, dset
    integer(hsize_t), dimension(3) :: dims, data_count, data_offset
    integer :: comm, ndims

    comm = MPI_COMM_WORLD
    info = MPI_INFO_NULL

    ndims = 3

    ! Set the size of the global array
    dims(1)=nx
    dims(2)=nym
    dims(3)=nzm

    ! Create a dataspace for the (global) array in the file
    call h5screate_simple_f(ndims, dims, filespace, hdf_error)

    ! Set the size of the data local to this process
    data_count(1) = nx
    data_count(2) = xend(2)-xstart(2)+1
    data_count(3) = xend(3)-xstart(3)+1

    ! Set the offset position of the local data in the global array
    data_offset(1) = 0
    data_offset(2) = xstart(2)-1
    data_offset(3) = xstart(3)-1

    ! Create a property list and set it up with the MPI comm
    call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf_error)
    call h5pset_fapl_mpio_f(plist_id, comm, info, hdf_error)

    ! Create the HDF5 file
    call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, &
            hdf_error, access_prp=plist_id)
    ! Close the property list used to set the file access property
    call h5pclose_f(plist_id, hdf_error)

    ! Create a dataset with name 'var'
    call h5dcreate_f(file_id, 'var', H5T_NATIVE_DOUBLE, &
            filespace, dset, hdf_error)

    ! Create a dataspace for the array in the (local) memory
    call h5screate_simple_f(ndims, data_count, memspace, hdf_error) 

    ! Get dataspace of file and select slab (subsection) from offsets
    call h5dget_space_f(dset, slabspace, hdf_error)
    call h5sselect_hyperslab_f(slabspace, H5S_SELECT_SET_F, &
            data_offset, data_count, hdf_error)

    ! Create the property defining the collective write
    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
            hdf_error)
    
    ! Write the dataset
    call h5dwrite_f(dset, H5T_NATIVE_DOUBLE, &
        arr, dims, &
        hdf_error, file_space_id = slabspace, mem_space_id = memspace, &
        xfer_prp = plist_id)
    
    ! Close the property list, dataset, dataspace, and file
    call h5pclose_f(plist_id, hdf_error)
    call h5dclose_f(dset, hdf_error)
    call h5sclose_f(filespace, hdf_error)
    call h5sclose_f(memspace, hdf_error)
    call h5sclose_f(slabspace, hdf_error)
    call h5fclose_f(file_id, hdf_error)

end subroutine write_3D_array

!> Write out a 3D spectrum decomposed using spectral space
subroutine write_3D_spectrum(filename, arr)
    use decomp_2d_fft
    character(len=30), intent(in) :: filename
    real, dimension(1:nxm,sp%xst(2):sp%xen(2),sp%xst(3):sp%xen(3)), intent(in) :: arr

    integer(hid_t) :: file_id, filespace, slabspace, memspace, dset
    integer(hsize_t), dimension(3) :: dims, data_count, data_offset
    integer :: comm, ndims

    comm = MPI_COMM_WORLD
    info = MPI_INFO_NULL

    ndims = 3

    ! Set the size of the global array
    dims(1)=nxm
    dims(2)=nym
    dims(3)=nzm

    ! Create a dataspace for the (global) array in the file
    call h5screate_simple_f(ndims, dims, filespace, hdf_error)

    ! Set the size of the data local to this process
    data_count(1) = nxm
    data_count(2) = sp%xen(2)-sp%xst(2)+1
    data_count(3) = sp%xen(3)-sp%xst(3)+1

    ! Set the offset position of the local data in the global array
    data_offset(1) = 0
    data_offset(2) = sp%xst(2)-1
    data_offset(3) = sp%xst(3)-1

    ! Create a property list and set it up with the MPI comm
    call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf_error)
    call h5pset_fapl_mpio_f(plist_id, comm, info, hdf_error)

    ! Create the HDF5 file
    call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, &
            hdf_error, access_prp=plist_id)
    ! Close the property list used to set the file access property
    call h5pclose_f(plist_id, hdf_error)

    ! Create a dataset with name 'var'
    call h5dcreate_f(file_id, 'var', H5T_NATIVE_DOUBLE, &
            filespace, dset, hdf_error)

    ! Create a dataspace for the array in the (local) memory
    call h5screate_simple_f(ndims, data_count, memspace, hdf_error) 

    ! Get dataspace of file and select slab (subsection) from offsets
    call h5dget_space_f(dset, slabspace, hdf_error)
    call h5sselect_hyperslab_f(slabspace, H5S_SELECT_SET_F, &
            data_offset, data_count, hdf_error)

    ! Create the property defining the collective write
    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
            hdf_error)
    
    ! Write the dataset
    call h5dwrite_f(dset, H5T_NATIVE_DOUBLE, &
        arr, dims, &
        hdf_error, file_space_id = slabspace, mem_space_id = memspace, &
        xfer_prp = plist_id)
    
    ! Close the property list, dataset, dataspace, and file
    call h5pclose_f(plist_id, hdf_error)
    call h5dclose_f(dset, hdf_error)
    call h5sclose_f(filespace, hdf_error)
    call h5sclose_f(memspace, hdf_error)
    call h5sclose_f(slabspace, hdf_error)
    call h5fclose_f(file_id, hdf_error)

end subroutine write_3D_spectrum

end module h5_tools