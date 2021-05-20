!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: SetSalBCs.F90                                  !
!    CONTAINS: subroutine SetSalBCs                       !
!                                                         ! 
!    PURPOSE: Initialization routine. Calcuates the       !
!     salinity boundary conditions at the plates          !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine SetSalBCs
    use param
    use decomp_2d
    implicit none
    integer :: ic,jc

    if (rays>=0) then ! unstable S gradient
        do ic=xstartr(3),xendr(3)
            do jc=xstartr(2),xendr(2)
                saltp(1,jc,ic)=0.5d0
                salbp(1,jc,ic)=-0.5d0
            end do
        end do
    else              ! stable S gradient
        do ic=xstartr(3),xendr(3)
            do jc=xstartr(2),xendr(2)
                saltp(1,jc,ic)=-0.5d0
                salbp(1,jc,ic)=0.5d0
            end do
        end do
    end if
    !CJH Add halo for interpolation routine
    call update_halo(saltp,lvlhalo)
    call update_halo(salbp,lvlhalo)

    return
end subroutine SetSalBCs