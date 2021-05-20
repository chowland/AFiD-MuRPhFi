!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: CreateICSal.F90                                !
!    CONTAINS: subroutine CreateICSal                     !
!                                                         ! 
!    PURPOSE: Initialization routine. Sets initial        !
!     conditions for salinity field                       !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine CreateICSal
    use param
    use mgrd_arrays, only: sal
    use decomp_2d, only: xstartr,xendr
    use mpih
    implicit none
    integer :: i,k,j
    real :: xxx,yyy,eps,varptb,amp

    call random_seed()
    eps=5d-2
    ! Assign linear salinity profile in the nodes k=1 to k=nxmr
    do i=xstartr(3),xendr(3)
        do j=xstartr(2),xendr(2)
            do k=1,nxmr
                call random_number(varptb)
                sal(k,j,i) = salbp(1,j,i) - (salbp(1,j,i) - saltp(1,j,i))*xmr(j)/xcr(nxr)
                if (abs(xmr(k)-0.5) + eps > 0.5) then
                    amp = 0.5 - abs(xmr(k)-0.5) ! CJH Prevent values of |S| exceeding 0.5
                    sal(k,j,i) = sal(k,j,i) + amp*(2.d0*varptb - 1.d0)
                else
                    sal(k,j,i) = sal(k,j,i) + eps*(2.d0*varptb - 1.d0)
                end if
            end do
        end do
    end do

    if (melt) then
        do i=xstartr(3),xendr(3)
            do j=xstartr(2),xendr(2)
                do k=1,nxmr
                    call random_number(varptb)
                    sal(k,j,i) = eps*(2.d0*varptb - 1.d0) * exp(-xmr(k)/0.1)
                end do
            end do
        end do
    end if
    
    return
end subroutine CreateICSal