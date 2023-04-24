!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: MakeMovieZCut.F90                              !
!    CONTAINS: subroutine Mkmov_zcut                      !
!                                                         ! 
!    PURPOSE: Write down zcut snapshot for flow vars      !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine Mkmov_zcut
    use param, only: nzm, nzmr
    use mpih
    use hdf5
    use decomp_2d, only: xstart,xend,xstartr,xendr,DECOMP_2D_COMM_CART_X
    use local_arrays, only: vz,vy,vx,temp, pr
    use mgrd_arrays, only: sal, phi, phic, tempr
    use moisture, only: humid
    use h5_tools
    implicit none
    character(70) :: filename
    character(4) :: varname

    integer(HID_T) :: file_id
    integer :: ic, comm

    logical :: fexist

    ! Select mid-plane
    ic = nzm/2 + 1

    ! Record filename as string
    filename = trim("outputdir/flowmov/movie_zcut.h5")

    ! MPI
    call MPI_CART_SUB(DECOMP_2D_COMM_CART_X, (/.true.,.false./), comm, ierr)
    info = MPI_INFO_NULL

    if(ic.le.xend(3) .and. ic.ge.xstart(3)) then

        call h5_open_or_create(file_id, filename, comm, fexist)
        if (.not. fexist) call h5_add_slice_groups(file_id)

        varname = 'vx'
        call write_H5_plane(file_id, varname, vx(1:nx, xstart(2):xend(2), ic), 'z')

        varname = 'vy'
        call write_H5_plane(file_id, varname, vy(1:nxm, xstart(2):xend(2), ic), 'z')
        
        varname = 'vz'
        call write_H5_plane(file_id, varname, vz(1:nxm, xstart(2):xend(2), ic), 'z')
        
        varname = 'temp'
        call write_H5_plane(file_id, varname, temp(1:nxm, xstart(2):xend(2), ic), 'z')
        
        if (moist) then
            varname = 'qhum'
            call write_H5_plane(file_id, varname, humid(1:nxm, xstart(2):xend(2), ic), 'z')
        end if

        call h5fclose_f(file_id, hdf_error)

    end if

    if (salinity) then
        !! Repeat on refined grid to save salinity
        ! Select midplane index for refined grid
        ic = nzmr/2 + 1

        if(ic.le.xendr(3) .and. ic.ge.xstartr(3)) then

            call h5_open_or_create(file_id, filename, comm, fexist)
            if (.not. fexist) call h5_add_slice_groups(file_id)
    
            varname = 'sal'
            call write_H5_plane(file_id, varname, sal(1:nxmr, xstartr(2):xendr(2), ic), 'z')
            
            call h5fclose_f(file_id, hdf_error)
    
        end if
    end if

    if (phasefield) then
        !! Repeat on refined grid to save phi
        ! Select midplane index for refined grid
        ic = nzmr/2 + 1

        if(ic.le.xendr(3) .and. ic.ge.xstartr(3)) then

            call h5_open_or_create(file_id, filename, comm, fexist)
            if (.not. fexist) call h5_add_slice_groups(file_id)
    
            varname = 'phi'
            call write_H5_plane(file_id, varname, phi(1:nxmr, xstartr(2):xendr(2), ic), 'z')
            
            call h5fclose_f(file_id, hdf_error)
    
        end if
    end if

    call MpiBarrier

end subroutine Mkmov_zcut