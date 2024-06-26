subroutine ReadFlowInterp(prow,pcol)
    use mpih
    use decomp_2d
    use local_arrays
    use param
    use input_grids
    use afid_salinity
    use afid_phasefield
    use afid_moisture, only: humid, InterpInputHum
    use AuxiliaryRoutines
    
    implicit none

    integer, intent(in) :: prow, pcol

    character*70 :: filnam,dsetname

    real :: yleno, zleno

    logical :: fexist

    ! Read old grid information
    filnam = trim("outputdir/continua_master.h5")

    inquire(file=filnam, exist=fexist)
    if (.not.fexist) then
        write(*,*) "Continuation files not found"
        call MpiAbort
    end if

    if (ismaster) then
        dsetname = trim("nx")
        call HdfSerialReadIntScalar(dsetname, filnam, nxo)
        dsetname = trim("ny")
        call HdfSerialReadIntScalar(dsetname, filnam, nyo)
        dsetname = trim("nz")
        call HdfSerialReadIntScalar(dsetname, filnam, nzo)
        dsetname = trim("time")
        call HdfSerialReadRealScalar(dsetname, filnam, time)
        dsetname = trim("istr3")
        call HdfSerialReadIntScalar(dsetname, filnam, istr3o)
        dsetname = trim("str3")
        call HdfSerialReadRealScalar(dsetname, filnam, str3o)
        dsetname = trim("ylen")
        call HdfSerialReadRealScalar(dsetname, filnam, yleno)
        dsetname = trim("zlen")
        call HdfSerialReadRealScalar(dsetname, filnam, zleno)
        if (multires) then
            dsetname = trim("nxr")
            call HdfSerialReadIntScalar(dsetname, filnam, nxro)
            dsetname = trim("nyr")
            call HdfSerialReadIntScalar(dsetname, filnam, nyro)
            dsetname = trim("nzr")
            call HdfSerialReadIntScalar(dsetname, filnam, nzro)
            dsetname = trim("istr3r")
            call HdfSerialReadIntScalar(dsetname, filnam, istr3ro)
        end if
    end if

    call MpiBarrier
    call MpiBcastInt(nxo)
    call MpiBcastInt(nyo)
    call MpiBcastInt(nzo)
    call MpiBcastInt(istr3o)
    call MpiBcastReal(str3o)
    call MpiBcastReal(time)
    call MpiBcastReal(yleno)
    call MpiBcastReal(zleno)
    if (multires) then
        call MpiBcastInt(nxro)
        call MpiBcastInt(nyro)
        call MpiBcastInt(nzro)
        call MpiBcastInt(istr3ro)
    end if

    if ((abs(ylen-yleno)>1e-8) .or. (abs(zlen-zleno)>1e-8)) then
        write(*,*) "Continua domain size does not match bou.in"
        write(*,*) "old Ly, Lz: ",yleno, zleno
        write(*,*) "current Ly, Lz: ", ylen, zlen
        call MpiAbort
    end if

    ! Check whether we are using a new grid
    if ((nx.ne.nxo) .or. (ny.ne.nyo) .or. (nz.ne.nzo) .or. &
        (istr3/=istr3o) .or. (abs(str3-str3o)>1e-8) .or. &
        (multires .and. ( &
            (nxr/=nxro) .or. (nyr/=nyro) .or. &
            (nzr/=nzro) .or. (istr3r/=istr3ro)) &
        )) then

        ! Check that old grid is consistent with the pencil decomposition
        if (ismaster) then
            if ((modulo(nyo-1,prow) /= 0) .or. (modulo(nzo-1,pcol) /= 0)) then
                write(*,*) "********** WARNING **********"
                write(*,*) "Grid size of old input files not a perfect factor of the pencil decomposition"
                write(*,*) "This would cause severe issues with the interpolation to the new grid"
                write(*,*) "Terminating simulation..."
                call MpiAbort
            end if
            if (multires) then
                if ((modulo(nyro-1,prow)/= 0) .or. (modulo(nzro-1,pcol)/= 0)) then
                    write(*,*) "********** WARNING **********"
                    write(*,*) "Grid size of old input files not a perfect factor of the pencil decomposition"
                    write(*,*) "This would cause severe issues with the interpolation to the new grid"
                    write(*,*) "Terminating simulation..."
                    call MpiAbort
                end if
            end if
        end if
        
        call InitInputVars

        call CreateOldGrid

        call CreateInputStencil

        call MpiBarrier

        call partition(nxmo, nymo, nzmo, (/ 1,2,3 /), &
                    xstarto, xendo, xsizeo)

        call InterpInputVel

        if (moist) call InterpInputHum

        if (ismaster) write(*,*) "Velocity and T interpolated"

        if (multires) then
            call CreateSalStencil

            call MpiBarrier

            call partition(nxmro, nymro, nzmro, (/ 1,2,3 /), &
                        xstarto, xendo, xsizeo)

            if (salinity) call InterpInputSal
            if (phasefield) call InterpInputPhi

        end if

        call DeallocateInputVars

    else
        ! No interpolation necessary
        call HdfReadContinua(nz, ny, nx, xstart(2), xend(2), &
                                xstart(3), xend(3), 1, vx)
        call HdfReadContinua(nz, ny, nx, xstart(2), xend(2), &
                                xstart(3), xend(3), 2, vy)
        call HdfReadContinua(nz, ny, nx, xstart(2), xend(2), &
                                xstart(3), xend(3), 3, vz)
        call HdfReadContinua(nz, ny, nx, xstart(2), xend(2), &
                                xstart(3), xend(3), 4, temp)
        if (salinity) then
            call HdfReadContinua(nzr, nyr, nxr, xstartr(2), xendr(2), &
                                xstartr(3), xendr(3), 5, sal)
        end if
        if (phasefield) then
            call HdfReadContinua(nzr, nyr, nxr, xstartr(2), xendr(2), &
                                xstartr(3), xendr(3), 6, phi)
        end if
        filnam = trim("outputdir/continua_pr.h5")
        inquire(file=filnam, exist=fexist)
        if (fexist) then
            call HdfReadContinua(nz, ny, nx, xstart(2), xend(2), &
                                    xstart(3), xend(3), 7, pr)
        end if
        if (moist) then
            call HdfReadContinua(nz, ny, nx, xstart(2), xend(2), &
                                    xstart(3), xend(3), 8, humid)
        end if
    end if

    if (resetlogstime) time = 0.d0
    
    ! ! Increase max sim time by the end time in continua files
    ! tmax = tmax + time

    return

end subroutine ReadFlowInterp