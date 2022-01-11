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
    integer :: i,k,j, kmid
    real :: xxx,yyy,eps,varptb,amp
    real :: B, gamma, t0, x0, h0

    call random_seed()
    eps=5d-2
    if ((active_S==1) .and. (active_T==1)) then
        ! For DDC initialise with uniform salinity
        do i=xstartr(3),xendr(3)
            do j=xstartr(2),xendr(2)
                do k=1,nxmr
                    call random_number(varptb)
                    sal(k,j,i) = eps*(2.0*varptb - 1.0)
                end do
            end do
        end do
    else if ((RayS < 0) .and. (RayT < 0)) then
        eps = 1e-2
        do i=xstartr(3),xendr(3)
            do j=xstartr(2),xendr(2)
                do k=1,nxmr
                    sal(k,j,i) = -tanh(xmr(k) - 0.5*alx3)
                    call random_number(varptb)
                    sal(k,j,i) = sal(k,j,i) + &
                            cosh(xmr(k) - 0.5*alx3)**(-2)*eps*(2.0*varptb - 1.0)
                end do
            end do
        end do
    else
        ! Assign linear salinity profile in the nodes k=1 to k=nxmr
        do i=xstartr(3),xendr(3)
            do j=xstartr(2),xendr(2)
                do k=1,nxmr
                    call random_number(varptb)
                    sal(k,j,i) = salbp(1,j,i) - (salbp(1,j,i) - saltp(1,j,i))*xmr(k)/xcr(nxr)
                    if (abs(xmr(k)-0.5) + eps > 0.5) then
                        amp = 0.5 - abs(xmr(k)-0.5) ! CJH Prevent values of |S| exceeding 0.5
                        sal(k,j,i) = sal(k,j,i) + amp*(2.d0*varptb - 1.d0)
                    else
                        sal(k,j,i) = sal(k,j,i) + eps*(2.d0*varptb - 1.d0)
                    end if
                end do
            end do
        end do
    end if

    if (melt) then
        B = 0.77512
        do i=xstartr(3),xendr(3)
            do j=xstartr(2),xendr(2)
                do k=1,nxmr
                    ! call random_number(varptb)
                    ! sal(k,j,i) = eps*(2.d0*varptb - 1.d0) * exp(-xmr(k)/0.1)
                    sal(k,j,i) = 1.0 - B*erfc(xmr(k)*sqrt(pecs)/2.0)
                end do
            end do
        end do
    end if

    if (phasefield) then
        if (pf_IC==1) then
            t0 = 1
            x0 = 0.8
            gamma = 0.62428 !0.76023
            h0 = x0 + 2*gamma*sqrt(t0/pecs)
            B = 0.44748 !0.37372
            do i=xstartr(3),xendr(3)
                do j=xstartr(2),xendr(2)
                    do k=1,nxmr
                        ! if (xmr(k) <= h0) then
                            sal(k,j,i) = 1.0 - B*erfc((x0 - xmr(k))*sqrt(pecs/t0)/2.0)
                        ! else
                        !     sal(k,j,i) = 1.0 - B*erfc(-gamma)
                        ! end if
                    end do
                end do
            end do
        else
            kmid = nxmr/2
            do i=xstartr(3),xendr(3)
                do j=xstartr(2),xendr(2)
                    do k=1,kmid
                        call random_number(varptb)
                        sal(k,j,i) = 1.0 + eps*(2.d0*varptb - 1.d0)
                    end do
                    do k=kmid+1,nxmr
                        sal(k,j,i) = 0.0
                    end do
                end do
            end do
        end if
    end if
    
    return
end subroutine CreateICSal