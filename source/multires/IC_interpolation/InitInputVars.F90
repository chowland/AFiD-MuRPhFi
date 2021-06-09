!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: InitInputVars.F90                              !
!    CONTAINS: subroutine InitInputVars                   !
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

    return

end subroutine InitInputVars