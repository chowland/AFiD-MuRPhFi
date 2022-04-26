!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: DeallocateIBMVars.F90                          !
!    CONTAINS: subroutine DeallocateIBMVariables          !
!                                                         ! 
!    PURPOSE: Finalization routine. Deallocates all       !
!     variables used when including salinity              !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine DeallocateIBMVariables
    use param, only: salinity
    use ibm_param
    use AuxiliaryRoutines
    implicit none

    ! Main array
    call DestroyReal3DArray(forclo)
    if (salinity) then
        call DestroyReal3DArray(forclor)
        if(allocated(solidr)) deallocate(solidr)
    end if

    return
end subroutine DeallocateIBMVariables