!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: MakeMovieXCut.F90                              !
!    CONTAINS: subroutine Mkmov_xcut                      !
!                                                         ! 
!    PURPOSE: Write down xcut snapshot for flow vars      !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine Mkmov_xcut
    use mpih
    use hdf5
    use decomp_2d, only: xstart,xend,xstartr,xendr,DECOMP_2D_COMM_CART_X
    use local_arrays, only: vz,vy,vx,temp
    use mgrd_arrays, only: sal,phi
    use h5_tools
    use param, only: IBM
    implicit none
    character(70) :: filename
    character(4) :: varname

    integer(HID_T) :: file_id
    integer :: ic, comm

    logical :: fexist

    ! Select plane - plane next to lower wall
    ic = 1
    if (IBM) ic = nxm/2

    ! Record filename as string
    filename = trim("outputdir/flowmov/movie_xcut.h5")

    ! MPI
    call MPI_CART_SUB(DECOMP_2D_COMM_CART_X, (/.true.,.true./), comm, ierr)
    info = MPI_INFO_NULL

    call h5_open_or_create(file_id, filename, comm, fexist)
    if (.not. fexist) call h5_add_slice_groups(file_id)

    varname = 'vx'
    call write_H5_plane(file_id, varname, vx(max(ic,2), xstart(2):xend(2), xstart(3):xend(3)), 'x')

    varname = 'vy'
    call write_H5_plane(file_id, varname, vy(ic, xstart(2):xend(2), xstart(3):xend(3)), 'x')

    varname = 'vz'
    call write_H5_plane(file_id, varname, vz(ic, xstart(2):xend(2), xstart(3):xend(3)), 'x')

    varname = 'temp'
    call write_H5_plane(file_id, varname, temp(ic, xstart(2):xend(2), xstart(3):xend(3)), 'x')

    call h5fclose_f(file_id, hdf_error)

    if (salinity) then
        !! Repeat on refined grid to save salinity
        ! Select wall plane index for refined grid
        ic = 1
        if (IBM) ic = nxmr/2

        call h5_open_or_create(file_id, filename, comm, fexist)
        if (.not. fexist) call h5_add_slice_groups(file_id)

        varname = 'sal'
        call write_H5_plane(file_id, varname, sal(ic, xstartr(2):xendr(2), xstartr(3):xendr(3)), 'x')
        if (IBM) then
            varname = 'sal2'
            ic = ic + 52
            call write_H5_plane(file_id, varname, sal(ic, xstartr(2):xendr(2), xstartr(3):xendr(3)), 'x')
        end if

        call h5fclose_f(file_id, hdf_error)

    end if

    call MpiBarrier

end subroutine Mkmov_xcut