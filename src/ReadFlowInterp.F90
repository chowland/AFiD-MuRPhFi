subroutine ReadFlowInterp
    use mpih
    use decomp_2d
    use local_arrays
    use param
    use input_grids
    use mgrd_arrays
    use AuxiliaryRoutines
    
    implicit none

    character*70 :: filnam,dsetname

    integer :: i, j, k, ic, jc, kc

    real :: yleno, zleno
    real, allocatable, dimension(:,:,:) :: salo, tempo, vxo, vyo, vzo
    real, dimension(4,4,4) :: qv3
    real, dimension(4,4) :: qv2
    real, dimension(4) :: qv1

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
        dsetname = trim("nxr")
        call HdfSerialReadIntScalar(dsetname, filnam, nxro)
        dsetname = trim("ny")
        call HdfSerialReadIntScalar(dsetname, filnam, nyo)
        dsetname = trim("nyr")
        call HdfSerialReadIntScalar(dsetname, filnam, nyro)
        dsetname = trim("nz")
        call HdfSerialReadIntScalar(dsetname, filnam, nzo)
        dsetname = trim("nzr")
        call HdfSerialReadIntScalar(dsetname, filnam, nzro)
        dsetname = trim("time")
        call HdfSerialReadRealScalar(dsetname, filnam, time)
        dsetname = trim("istr3")
        call HdfSerialReadIntScalar(dsetname, filnam, istr3o)
        dsetname = trim("istr3r")
        call HdfSerialReadIntScalar(dsetname, filnam, istr3ro)
        dsetname = trim("str3")
        call HdfSerialReadRealScalar(dsetname, filnam, str3o)
        dsetname = trim("ylen")
        call HdfSerialReadRealScalar(dsetname, filnam, yleno)
        dsetname = trim("zlen")
        call HdfSerialReadRealScalar(dsetname, filnam, zleno)
    end if

    call MpiBarrier
    call MpiBcastInt(nxo)
    call MpiBcastInt(nxro)
    call MpiBcastInt(nyo)
    call MpiBcastInt(nyro)
    call MpiBcastInt(nzo)
    call MpiBcastInt(nzro)
    call MpiBcastInt(istr3o)
    call MpiBcastInt(istr3ro)
    call MpiBcastReal(str3o)
    call MpiBcastReal(time)
    call MpiBcastReal(yleno)
    call MpiBcastReal(zleno)

    if ((abs(ylen-yleno)>1e-8) .or. (abs(zlen-zleno)>1e-8)) then
        write(*,*) "Continua domain size does not match bou.in"
        write(*,*) "old Ly, Lz: ",yleno, zleno
        write(*,*) "current Ly, Lz: ", ylen, zlen
        call MpiAbort
    end if

    if (ismaster) then
        write(*,*) "nx, nxo: ", nx, nxo
        write(*,*) "ny, nyo: ", ny, nyo
        write(*,*) "nz, nzo: ", nz, nzo
        write(*,*) "nxr, nxro: ", nxr, nxro
        write(*,*) "nyr, nyro: ", nyr, nyro
        write(*,*) "nzr, nzro: ", nzr, nzro
    end if

    ! Check whether we are using a new grid
    if ((nx.ne.nxo) .or. (ny.ne.nyo) .or. (nz.ne.nzo) .or. &
        (nxr/=nxro) .or. (nyr/=nyro) .or. (nzr/=nzro) .or. &
        (istr3/=istr3o) .or. (istr3r/=istr3ro) .or. &
        (abs(str3-str3o)>1e-8)) then

        ! Allocate old grids
        call AllocateReal1DArray(xco, 1, nxo)
        call AllocateReal1DArray(xmo, 1, nxo)
        call AllocateReal1DArray(yco, 1, nyo)
        call AllocateReal1DArray(ymo, 1, nyo)
        call AllocateReal1DArray(zco, 1, nzo)
        call AllocateReal1DArray(zmo, 1, nzo)

        call AllocateReal1DArray(xcro, 1, nxro)
        call AllocateReal1DArray(xmro, 0, nxro+1)
        call AllocateReal1DArray(ycro, 1, nyro)
        call AllocateReal1DArray(ymro, 0, nyro+1)
        call AllocateReal1DArray(zcro, 1, nzro)
        call AllocateReal1DArray(zmro, 0, nzro+1)

        call CreateOldGrid

        if (ismaster) then
            call HdfCreateBlankFile("grids_test.h5")
            call HdfSerialWriteReal1D("xco", "grids_test.h5", xco, nxo)
            call HdfSerialWriteReal1D("xcro", "grids_test.h5", xcro, nxro)
            call HdfSerialWriteReal1D("xmo", "grids_test.h5", xmo, nxo)
            call HdfSerialWriteReal1D("xmro", "grids_test.h5", xmro, nxro)
            call HdfSerialWriteReal1D("ycro", "grids_test.h5", ycro, nyro)
            call HdfSerialWriteReal1D("ymro", "grids_test.h5", ymro, nyro)
        end if

        call CreateInputStencil

        if (ismaster) then
            call HdfCreateBlankFile("itp_range.h5")
            call HdfSerialWriteInt1D("irangs", "itp_range.h5", irangs, nx+1)
            call HdfSerialWriteInt1D("jrangs", "itp_range.h5", jrangs, ny+1)
            call HdfSerialWriteInt1D("krangs", "itp_range.h5", krangs, nz+1)
            call HdfSerialWriteInt1D("irangc", "itp_range.h5", irangc, nx+1)
            call HdfSerialWriteInt1D("jrangc", "itp_range.h5", jrangc, ny+1)
            call HdfSerialWriteInt1D("krangc", "itp_range.h5", krangc, nz+1)
        end if

        call MpiBarrier
        xs2o = ceiling(real((xstart(2) - 1)*nymo/nym)) + 1
        xe2o = ceiling(real(xend(2)*nymo/nym))
        xs3o = ceiling(real((xstart(3) - 1)*nzmo/nzm)) + 1
        xe3o = ceiling(real(xend(3)*nzmo/nzm))
        xs2o = max(xs2o, 1)
        xe2o = min(xe2o, nymo)
        xs3o = max(xs3o, 1)
        xe3o = min(xe3o, nzmo)

        write(*,*) "xs2o, xe2o: ",xs2o, xe2o
        write(*,*) "xs2, xe2: ", xstart(2), xend(2)
        write(*,*) "xs3o, xe3o: ",xs3o, xe3o

        call InterpInputVel

        if (ismaster) write(*,*) "Velocity and T interpolated"

        call CreateSalStencil

        if (ismaster) then
            ! call HdfCreateBlankFile("itp_range.h5")
            call HdfSerialWriteInt1D("irangsr", "itp_range.h5", irangs, nx+1)
            call HdfSerialWriteInt1D("jrangsr", "itp_range.h5", jrangs, ny+1)
            call HdfSerialWriteInt1D("krangsr", "itp_range.h5", krangs, nz+1)
            call HdfSerialWriteInt1D("irangcr", "itp_range.h5", irangc, nx+1)
            call HdfSerialWriteInt1D("jrangcr", "itp_range.h5", jrangc, ny+1)
            call HdfSerialWriteInt1D("krangcr", "itp_range.h5", krangc, nz+1)
        end if
 
        call MpiBarrier
        xs2o = ceiling(real((xstartr(2) - 1)*nymro/nymr)) + 1
        xe2o = ceiling(real(xendr(2)*nymro/nymr))
        xs3o = ceiling(real((xstartr(3) - 1)*nzmro/nzmr)) + 1
        xe3o = ceiling(real(xendr(3)*nzmro/nzmr))
        xs2o = max(xs2o, 1)
        xe2o = min(xe2o, nymro)
        xs3o = max(xs3o, 1)
        xe3o = min(xe3o, nzmro)

        write(*,*) "Rxs2o, xe2o: ",xs2o, xe2o
        write(*,*) "Rxs2, xe2: ", xstartr(2), xendr(2)
        write(*,*) "Rxs3o, xe3o: ",xs3o, xe3o

        call InterpInputSal

        ! Deallocate old grids
        call DestroyReal1DArray(xco)
        call DestroyReal1DArray(xmo)
        call DestroyReal1DArray(yco)
        call DestroyReal1DArray(ymo)
        call DestroyReal1DArray(zco)
        call DestroyReal1DArray(zmo)

        call DestroyReal1DArray(xcro)
        call DestroyReal1DArray(xmro)
        call DestroyReal1DArray(ycro)
        call DestroyReal1DArray(ymro)
        call DestroyReal1DArray(zcro)
        call DestroyReal1DArray(zmro)

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
        call HdfReadContinua(nzr, nyr, nxr, xstartr(2), xendr(2), &
                                xstartr(3), xendr(3), 5, sal)
    end if

    if (resetlogstime) time = 0.d0
    
    ! Increase max sim time by the end time in continua files
    tmax = tmax + time

    return

end subroutine ReadFlowInterp