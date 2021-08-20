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
    end if

    return

end subroutine InitInputVars

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine is needed for interpolating an initial condition
! when MULTIRES is not enabled. This subroutine simply allocates
! the index variable needed to perform the interpolation that is
! usually allocated by InitMgrdVariables when MULTIRES is enabled
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine InitInputIndices
    use param
    use mgrd_arrays

    implicit none

    ! Allocate indices for interpolation
    call AllocateInt1dArray(irangs,0,nxo)
    call AllocateInt1dArray(jrangs,0,nyo)
    call AllocateInt1dArray(krangs,0,nzo)

    call AllocateInt1dArray(irangc,0,nxo)
    call AllocateInt1dArray(jrangc,0,nyo)
    call AllocateInt1dArray(krangc,0,nzo)

    return

end subroutine InitInputIndices