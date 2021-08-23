!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: InitInputVars.F90                              !
!    CONTAINS: subroutine InitInputVars InitInputIndices  !
!                                                         ! 
!    PURPOSE: Initialization routine. Allocates memory    !
!           for arrays used for input interpolation       !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine InitInputVars
    use param
    use input_grids
    use mgrd_arrays
    use decomp_2d
    use AuxiliaryRoutines
    
    implicit none

    ! Allocate old grids
    call AllocateReal1DArray(xco, 1, nxo)
    call AllocateReal1DArray(xmo, 1, nxo)
    call AllocateReal1DArray(yco, 1, nyo)
    call AllocateReal1DArray(ymo, 1, nyo)
    call AllocateReal1DArray(zco, 1, nzo)
    call AllocateReal1DArray(zmo, 1, nzo)

    if (multires) then
        call AllocateReal1DArray(xcro, 1, nxro)
        call AllocateReal1DArray(xmro, 0, nxro+1)
        call AllocateReal1DArray(ycro, 1, nyro)
        call AllocateReal1DArray(ymro, 0, nyro+1)
        call AllocateReal1DArray(zcro, 1, nzro)
        call AllocateReal1DArray(zmro, 0, nzro+1)

        ! Allocate extended index tracking for old refined grid
        call AllocateInt1DArray(irangsr,0,nxro)
        call AllocateInt1DArray(jrangsr,0,nyro)
        call AllocateInt1DArray(krangsr,0,nzro)
    else
        ! If multires is not enabled, allocate all the variables
        ! necessary to build the interpolation stencil
        call AllocateInt1DArray(irangs,0,nxo)
        call AllocateInt1DArray(jrangs,0,nyo)
        call AllocateInt1DArray(krangs,0,nzo)

        call AllocateInt1DArray(irangc,0,nxo)
        call AllocateInt1DArray(jrangc,0,nyo)
        call AllocateInt1DArray(krangc,0,nzo)

        call AllocateReal2DArray(cxvx,1,4,0,nx)
        call AllocateReal2DArray(cxvy,1,4,0,nx)
        call AllocateReal2DArray(cxvz,1,4,0,nx)

        call AllocateReal2DArray(cyvx,1,4,0,ny)
        call AllocateReal2DArray(cyvy,1,4,0,ny)
        call AllocateReal2DArray(cyvz,1,4,0,ny)

        call AllocateReal2DArray(czvx,1,4,0,nz)
        call AllocateReal2DArray(czvy,1,4,0,nz)
        call AllocateReal2DArray(czvz,1,4,0,nz)

        call AllocateReal2DArray(cxrs,1,4,0,nx)
        call AllocateReal2DArray(cyrs,1,4,0,ny)
        call AllocateReal2DArray(czrs,1,4,0,nz)

        call AllocateReal3DArray(tpdv,-1,nx+1,xstart(2)-2,xend(2)+2,xstart(3)-2,xend(3)+2)
    end if

    return

end subroutine InitInputVars