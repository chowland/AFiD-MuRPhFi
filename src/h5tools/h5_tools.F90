module h5_tools
    use HDF5
    use mpih, only: MPI_INFO_NULL
    use param, only: salinity, phasefield, nx, nxm, nym, nzm, time, tframe, nxmr, nymr, nzmr
    use decomp_2d, only: xstart, xstartr
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
    integer(HID_T), intent(inout) :: file_id

    integer :: i
    character(len=4), dimension(4) :: var_list

    var_list = [character(len=4) :: "vx", "vy", "vz", "temp"]

    do i=1,size(var_list)
        call h5gcreate_f(file_id, trim(var_list(i)), group_id, hdf_error)
        call h5gclose_f(group_id, hdf_error)
    end do

    if (salinity) then
        call h5gcreate_f(file_id, "sal", group_id, hdf_error)
        call h5gclose_f(group_id, hdf_error)
    end if

    if (phasefield) then
        call h5gcreate_f(file_id, "phi", group_id, hdf_error)
        call h5gclose_f(group_id, hdf_error)
    end if
end subroutine

subroutine write_H5_plane(file_id, varname, var, axis)
    integer(HID_T), intent(inout) :: file_id
    character(len=4), intent(in) :: varname
    real, intent(in) :: var(:,:)
    character, intent(in) :: axis

    integer, parameter :: ndims = 2
    integer(HID_T) :: filespace, slabspace, memspace, dset
    integer(HSIZE_T), dimension(2) :: dims, data_count, data_offset
    character(len=5) :: frame
    character(len=10) :: dsetname
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
    if (index('phisal', trim(varname)) /= 0) then
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

end module h5_tools