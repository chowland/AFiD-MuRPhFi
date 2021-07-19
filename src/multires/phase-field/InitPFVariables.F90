!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: InitPFVariables.F90                            !
!    CONTAINS: subroutine InitPFVariables                 !
!                                                         ! 
!    PURPOSE: Initialization routine. Sets to zero all    !
!     variables loaded when using phase-field method      !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine InitPFVariables
    use param
    use mgrd_arrays
    use decomp_2d
    use AuxiliaryRoutines
    implicit none
    
    ! Main array with ghost cells
    call AllocateReal3DArray(phi,1,nxr,xstartr(2)-lvlhalo,xendr(2)+lvlhalo,xstartr(3)-lvlhalo,xendr(3)+lvlhalo)
    ! Refined temperature array
    call AllocateReal3DArray(tempr,1,nxr,xstartr(2)-lvlhalo,xendr(2)+lvlhalo,xstartr(3)-lvlhalo,xendr(3)+lvlhalo)

    ! Arrays without ghost cells
    call AllocateReal3DArray(ruphi,1,nxr,xstartr(2),xendr(2),xstartr(3),xendr(3))
    call AllocateReal3DArray(hphi,1,nxr,xstartr(2),xendr(2),xstartr(3),xendr(3))

    ! Coarse array for phi or d(phi)/dt
    call AllocateReal3DArray(phic,1,nx,xstart(2)-lvlhalo,xend(2)+lvlhalo,xstart(3)-lvlhalo,xend(3)+lvlhalo)

    return
end subroutine InitPFVariables