
subroutine mean_zplane
    use param, only: nx, nzm
    use mpih
    use hdf5
    use decomp_2d, only: xstart,xend,xsize,xstartr,xendr,DECOMP_2D_COMM_CART_X,nrank
    use local_arrays, only: vz,vy,vx,temp
    use mgrd_arrays, only: sal, phi
    use h5_tools
    use means
    implicit none

    character(70) :: filename
    character(4) :: varname
    real, allocatable :: uplane(:,:), vplane(:,:), wplane(:,:), Tplane(:,:)
    integer :: comm, zrank, i, j, bufsize
    integer(HID_T) :: file_id
    logical :: fexist

    filename = trim("outputdir/flowmov/movie_zmean.h5")

    allocate(uplane(1:nx ,xstart(2):xend(2)))
    allocate(vplane(1:nxm,xstart(2):xend(2)))
    allocate(wplane(1:nxm,xstart(2):xend(2)))
    allocate(Tplane(1:nxm,xstart(2):xend(2)))

    call MPI_CART_SUB(DECOMP_2D_COMM_CART_X, (/.false., .true./), comm, ierr)
    call MPI_COMM_RANK(comm, zrank, ierr)

    call zmean(vx(1:nx ,xstart(2):xend(2),xstart(3):xend(3)), uplane, comm, zrank)
    call zmean(vy(1:nxm,xstart(2):xend(2),xstart(3):xend(3)), vplane, comm, zrank)
    call zmean(vz(1:nxm,xstart(2):xend(2),xstart(3):xend(3)), wplane, comm, zrank)
    call zmean(temp(1:nxm,xstart(2):xend(2),xstart(3):xend(3)), Tplane, comm, zrank)

    call MPI_CART_SUB(DECOMP_2D_COMM_CART_X, (/.true., .false./), comm, ierr)

    if (zrank==0) then

        ! Open the movie file (create if it doesn't exist)
        call h5_open_or_create(file_id, filename, comm, fexist)
        ! If file only just created, add groups for data storage
        if (.not. fexist) call h5_add_slice_groups(file_id)

        varname='vx'
        call write_H5_plane(file_id, varname, uplane, 'z')
        varname='vy'
        call write_H5_plane(file_id, varname, vplane, 'z')
        varname='vz'
        call write_H5_plane(file_id, varname, wplane, 'z')
        varname='temp'
        call write_H5_plane(file_id, varname, Tplane, 'z')

        ! Close HDF5 file
        call h5fclose_f(file_id, hdf_error)
    end if

    if (allocated(uplane)) deallocate(uplane)
    if (allocated(vplane)) deallocate(vplane)
    if (allocated(wplane)) deallocate(wplane)
    if (allocated(Tplane)) deallocate(Tplane)

    if (salinity) then
        allocate(Tplane(1:nxmr,xstartr(2):xendr(2)))

        call MPI_CART_SUB(DECOMP_2D_COMM_CART_X, (/.false., .true./), comm, ierr)
        call MPI_COMM_RANK(comm, zrank, ierr)

        call zmeanr(sal(1:nxmr,xstartr(2):xendr(2),xstartr(3):xendr(3)), Tplane, comm, zrank)

        call MPI_CART_SUB(DECOMP_2D_COMM_CART_X, (/.true., .false./), comm, ierr)

        if (zrank==0) then
            ! Open the movie file (create if it doesn't exist)
            call h5_open_or_create(file_id, filename, comm, fexist)
            ! If file only just created, add groups for data storage
            if (.not. fexist) call h5_add_slice_groups(file_id)

            varname='sal'
            call write_H5_plane(file_id, varname, Tplane, 'z')
    
            ! Close HDF5 file
            call h5fclose_f(file_id, hdf_error)
        end if

        if (allocated(Tplane)) deallocate(Tplane)
    end if

    if (phasefield) then
        allocate(Tplane(1:nxmr,xstartr(2):xendr(2)))

        call MPI_CART_SUB(DECOMP_2D_COMM_CART_X, (/.false., .true./), comm, ierr)
        call MPI_COMM_RANK(comm, zrank, ierr)

        call zmeanr(phi(1:nxmr,xstartr(2):xendr(2),xstartr(3):xendr(3)), Tplane, comm, zrank)

        call MPI_CART_SUB(DECOMP_2D_COMM_CART_X, (/.true., .false./), comm, ierr)

        if (zrank==0) then
            ! Open the movie file (create if it doesn't exist)
            call h5_open_or_create(file_id, filename, comm, fexist)
            ! If file only just created, add groups for data storage
            if (.not. fexist) call h5_add_slice_groups(file_id)

            varname='sal'
            call write_H5_plane(file_id, varname, Tplane, 'z')
    
            ! Close HDF5 file
            call h5fclose_f(file_id, hdf_error)
        end if

        if (allocated(Tplane)) deallocate(Tplane)
    end if

end subroutine mean_zplane

subroutine mean_yplane
    use param, only: nx, nym
    use mpih
    use hdf5
    use decomp_2d, only: xstart,xend,xsize,xstartr,xendr,DECOMP_2D_COMM_CART_X,nrank
    use local_arrays, only: vz,vy,vx,temp
    use mgrd_arrays, only: sal, phi
    use h5_tools
    use means
    implicit none

    character(70) :: filename
    character(4) :: varname
    real, allocatable :: uplane(:,:), vplane(:,:), wplane(:,:), Tplane(:,:)
    integer :: comm, yrank, i, j, bufsize
    integer(HID_T) :: file_id
    logical :: fexist

    filename = trim("outputdir/flowmov/movie_ymean.h5")

    allocate(uplane(1:nx ,xstart(3):xend(3)))
    allocate(vplane(1:nxm,xstart(3):xend(3)))
    allocate(wplane(1:nxm,xstart(3):xend(3)))
    allocate(Tplane(1:nxm,xstart(3):xend(3)))

    call MPI_CART_SUB(DECOMP_2D_COMM_CART_X, (/.true., .false./), comm, ierr)
    call MPI_COMM_RANK(comm, yrank, ierr)

    call ymean(vx(1:nx ,xstart(2):xend(2),xstart(3):xend(3)), uplane, comm, yrank)
    call ymean(vy(1:nxm,xstart(2):xend(2),xstart(3):xend(3)), vplane, comm, yrank)
    call ymean(vz(1:nxm,xstart(2):xend(2),xstart(3):xend(3)), wplane, comm, yrank)
    call ymean(temp(1:nxm,xstart(2):xend(2),xstart(3):xend(3)), Tplane, comm, yrank)

    call MPI_CART_SUB(DECOMP_2D_COMM_CART_X, (/.false., .true./), comm, ierr)

    if (yrank==0) then

        ! Open the movie file (create if it doesn't exist)
        call h5_open_or_create(file_id, filename, comm, fexist)
        ! If file only just created, add groups for data storage
        if (.not. fexist) call h5_add_slice_groups(file_id)

        varname='vx'
        call write_H5_plane(file_id, varname, uplane, 'y')
        varname='vy'
        call write_H5_plane(file_id, varname, vplane, 'y')
        varname='vz'
        call write_H5_plane(file_id, varname, wplane, 'y')
        varname='temp'
        call write_H5_plane(file_id, varname, Tplane, 'y')

        ! Close HDF5 file
        call h5fclose_f(file_id, hdf_error)
    end if

    if (allocated(uplane)) deallocate(uplane)
    if (allocated(vplane)) deallocate(vplane)
    if (allocated(wplane)) deallocate(wplane)
    if (allocated(Tplane)) deallocate(Tplane)

    if (salinity) then
        allocate(Tplane(1:nxmr,xstartr(3):xendr(3)))

        call MPI_CART_SUB(DECOMP_2D_COMM_CART_X, (/.true., .false./), comm, ierr)
        call MPI_COMM_RANK(comm, yrank, ierr)

        call ymeanr(sal(1:nxmr,xstartr(2):xendr(2),xstartr(3):xendr(3)), Tplane, comm, yrank)

        call MPI_CART_SUB(DECOMP_2D_COMM_CART_X, (/.false., .true./), comm, ierr)

        if (yrank==0) then
            ! Open the movie file (create if it doesn't exist)
            call h5_open_or_create(file_id, filename, comm, fexist)
            ! If file only just created, add groups for data storage
            if (.not. fexist) call h5_add_slice_groups(file_id)

            varname='sal'
            call write_H5_plane(file_id, varname, Tplane, 'y')
    
            ! Close HDF5 file
            call h5fclose_f(file_id, hdf_error)
        end if

        if (allocated(Tplane)) deallocate(Tplane)
    end if

    if (phasefield) then
        allocate(Tplane(1:nxmr,xstartr(3):xendr(3)))

        call MPI_CART_SUB(DECOMP_2D_COMM_CART_X, (/.true., .false./), comm, ierr)
        call MPI_COMM_RANK(comm, yrank, ierr)

        call ymeanr(phi(1:nxmr,xstartr(2):xendr(2),xstartr(3):xendr(3)), Tplane, comm, yrank)

        call MPI_CART_SUB(DECOMP_2D_COMM_CART_X, (/.false., .true./), comm, ierr)

        if (yrank==0) then
            ! Open the movie file (create if it doesn't exist)
            call h5_open_or_create(file_id, filename, comm, fexist)
            ! If file only just created, add groups for data storage
            if (.not. fexist) call h5_add_slice_groups(file_id)

            varname='sal'
            call write_H5_plane(file_id, varname, Tplane, 'y')
    
            ! Close HDF5 file
            call h5fclose_f(file_id, hdf_error)
        end if

        if (allocated(Tplane)) deallocate(Tplane)
    end if

end subroutine mean_yplane