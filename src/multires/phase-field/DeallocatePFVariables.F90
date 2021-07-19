!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: DeallocatePFVariables.F90                      !
!    CONTAINS: subroutine DeallocatePFVariables           !
!                                                         ! 
!    PURPOSE: Finalization routine. Deallocates all       !
!     variables used for phase-field method               !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine DeallocatePFVariables
    use param
    use mgrd_arrays
    use decomp_2d
    use AuxiliaryRoutines
    implicit none

    ! Main array
    call DestroyReal3DArray(phi)

    ! Array for refined temperature
    call DestroyReal3DArray(tempr)

    ! Arrays without ghost cells
    call DestroyReal3DArray(ruphi)
    call DestroyReal3DArray(hphi)

    ! Coarse array for phi or d(phi)/dt
    call DestroyReal3DArray(phic)

    return
end subroutine DeallocatePFVariables