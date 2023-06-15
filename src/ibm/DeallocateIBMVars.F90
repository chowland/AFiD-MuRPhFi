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
    ! call DestroyReal3DArray(forclo)
    if (allocated(ibmaskx)) deallocate(ibmaskx)
    if (allocated(ibmasky)) deallocate(ibmasky)
    if (allocated(ibmaskz)) deallocate(ibmaskz)
    if (allocated(ibmaskt)) deallocate(ibmaskt)

    if (salinity) then
        if(allocated(solidr)) deallocate(solidr)
        if (allocated(ibmaskr)) deallocate(ibmaskr)
    end if

    return
end subroutine DeallocateIBMVariables