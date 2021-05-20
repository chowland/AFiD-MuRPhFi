!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: DeallocateSalVariables.F90                     !
!    CONTAINS: subroutine DeallocateSalVariables          !
!                                                         ! 
!    PURPOSE: Finalization routine. Deallocates all       !
!     variables used when including salinity              !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine DeallocateSalVariables
    use param
    use mgrd_arrays
    use decomp_2d
    use AuxiliaryRoutines
    implicit none

    ! Boundary conditions
    call DestroyReal3DArray(salbp)
    call DestroyReal3DArray(saltp)

    ! Main array
    call DestroyReal3DArray(sal)

    ! RHS arrays
    call DestroyReal3DArray(rhsr)
    call DestroyReal3DArray(rusal)
    call DestroyReal3DArray(hsal)

    ! Coarse array
    call DestroyReal3DArray(salc)

    ! Extra T slice for melt condition
    if (melt) then
        call DestroyReal3DArray(tempr)
    end if

    return
end subroutine DeallocateSalVariables