!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: DeallocateInputVars.F90                        !
!    CONTAINS: subroutine DeallocateInputVars             !
!                                                         ! 
!    PURPOSE: Finalization routine. Deallocates all       !
!     variables used by the input interpolation           !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine DeallocateInputVars
    use param
    use input_grids
    use decomp_2d
    use AuxiliaryRoutines
    implicit none

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

    if (multires) then
        call DestroyInt1DArray(irangsr)
        call DestroyInt1DArray(jrangsr)
        call DestroyInt1DArray(krangsr)
    else
        call DestroyInt1DArray(irangs)
        call DestroyInt1DArray(jrangs)
        call DestroyInt1DArray(krangs)
        call DestroyInt1DArray(irangc)
        call DestroyInt1DArray(jrangc)
        call DestroyInt1DArray(krangc)
    end if

    return

end subroutine DeallocateInputVars

! subroutine DeallocateInputIndices
!     use mgrd_arrays
!     use AuxiliaryRoutines
!     implicit none

!     call DestroyInt1DArray(irangs)
!     call DestroyInt1DArray(jrangs)
!     call DestroyInt1DArray(krangs)
!     call DestroyInt1DArray(irangc)
!     call DestroyInt1DArray(jrangc)
!     call DestroyInt1DArray(krangc)

! end subroutine DeallocateInputIndices