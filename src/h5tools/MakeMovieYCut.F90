!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: MakeMovieYCut.F90                              !
!    CONTAINS: subroutine Mkmov_ycut                      !
!                                                         ! 
!    PURPOSE: Write down ycut snapshot for flow vars      !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine Mkmov_ycut
    use param, only: nym, nymr
    use mpih
    use hdf5
    use decomp_2d, only: xstart,xend,xstartr,xendr,DECOMP_2D_COMM_CART_X
    use local_arrays, only: vz,vy,vx,temp
    use afid_salinity, only: sal
    use h5_tools
    implicit none
    character(70) :: filename
    character(4) :: varname

    integer(HID_T) :: file_id
    integer :: ic, comm

    logical :: fexist

    ! Select mid-plane
    ic = nym/2 + 1

    ! Record filename as string
    filename = trim("outputdir/flowmov/movie_ycut.h5")

    ! MPI
    call MPI_CART_SUB(DECOMP_2D_COMM_CART_X, (/.false.,.true./), comm, ierr)
    info = MPI_INFO_NULL

    if (ic.le.xend(2) .and. ic.ge.xstart(2)) then

        call h5_open_or_create(file_id, filename, comm, fexist)
        if (.not. fexist) call h5_add_slice_groups(file_id)

        varname = 'vx'
        call write_H5_plane(file_id, varname, vx(1:nx, ic, xstart(3):xend(3)), 'y')

        varname = 'vy'
        call write_H5_plane(file_id, varname, vy(1:nxm, ic, xstart(3):xend(3)), 'y')

        varname = 'vz'
        call write_H5_plane(file_id, varname, vz(1:nxm, ic, xstart(3):xend(3)), 'y')

        varname = 'temp'
        call write_H5_plane(file_id, varname, temp(1:nxm, ic, xstart(3):xend(3)), 'y')

        call h5fclose_f(file_id, hdf_error)

    end if

    if (salinity) then
        !! Repeat on refined grid to save salinity
        ! Select midplane index for refined grid
        ic = nymr/2 + 1
        
        if (ic.le.xendr(2) .and. ic.ge.xstartr(2)) then

            call h5_open_or_create(file_id, filename, comm, fexist)
            if (.not. fexist) call h5_add_slice_groups(file_id)
    
            varname = 'sal'
            call write_H5_plane(file_id, varname, sal(1:nxmr, ic, xstartr(3):xendr(3)), 'y')
    
            call h5fclose_f(file_id, hdf_error)

        end if
    end if

    if (phasefield) then
        !! Repeat on refined grid to save salinity
        ! Select midplane index for refined grid
        ic = nymr/2 + 1
        
        if (ic.le.xendr(2) .and. ic.ge.xstartr(2)) then

            call h5_open_or_create(file_id, filename, comm, fexist)
            if (.not. fexist) call h5_add_slice_groups(file_id)
    
            varname = 'phi'
            call write_H5_plane(file_id, varname, phi(1:nxmr, ic, xstartr(3):xendr(3)), 'y')
    
            call h5fclose_f(file_id, hdf_error)

        end if
    end if

    call MpiBarrier

end subroutine Mkmov_ycut