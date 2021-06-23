!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
!    FILE: CreateICPF.F90                                 !
!    CONTAINS: subroutine CreateICPF                      !
!                                                         !
!    PURPOSE: Initialization routine. Sets initial        !
!     conditions for phase-field                          !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine CreateICPF

    use param
    use mgrd_arrays, only: phi
    use decomp_2d, only: xstartr,xendr
    use mpih

    implicit none

    integer :: i,j,k

    ! do i=xstartr(3),xendr(3)
    !     do j=xstartr(2),xendr(2)
    !         do k=1,nxmr
    !             phi(k,j,i) = 0.5*(1.0 + tanh((xmr(k) - 0.5)/2/pf_eps))
    !         end do
    !     end do
    ! end do

    do i=xstartr(3),xendr(3)
        do j=xstartr(2),xendr(2)
            do k=1,nxmr
                if ((xmr(k) - 0.75)**2 + (ymr(j) - ylen/2)**2 + (zmr(i) - zlen/2)**2 < 0.0225) then
                    phi(k,j,i) = 1.0
                else
                    phi(k,j,i) = 0.0
                end if
            end do
        end do
    end do

    return

end subroutine CreateICPF