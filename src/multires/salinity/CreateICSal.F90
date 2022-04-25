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
    real :: gamma, t0, x0, h0, A, B, alpha
    real, dimension(11) :: yh, zh

    call random_seed()
    eps=5e-3
    if (IBM) then
        do i=xstartr(3),xendr(3)
            do j=xstartr(2),xendr(2)
                do k=1,nxmr
                    sal(k,j,i) = 0.5*tanh((xmr(k) - 0.5 + 0.01*sin(2.0*pi*ymr(j)/ylen*5.0))/0.01)
                end do
            end do
        end do
    else if ((active_S==1) .and. (active_T==1) .and. (gAxis==1)) then
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
            call read_phase_field_params(A, B, alpha)
            t0 = 1e-3
            x0 = 0.8
            h0 = x0 + 2*alpha*sqrt(t0)
            do i=xstartr(3),xendr(3)
                do j=xstartr(2),xendr(2)
                    do k=1,nxmr
                        sal(k,j,i) = 1.0 - B*erfc((x0 - xmr(k))/sqrt(PraT/PraS*t0)/2.0)
                    end do
                end do
            end do
        else if (pf_IC==2) then
            call read_phase_field_params(A, B, alpha)
            t0 = 1e-3
            h0 = 0.1 - 2*alpha*sqrt(t0)
            do i=xstartr(3),xendr(3)
                do j=xstartr(2),xendr(2)
                    do k=1,nxmr
                        call random_number(varptb)
                        sal(k,j,i) = 1.0 - &
                            B*erfc((ylen/2.0 - h0 - ymr(j))/sqrt(PraT/PraS*t0)/2.0) + &
                            B*erfc((ylen/2.0 + h0 - ymr(j))/sqrt(PraT/PraS*t0)/2.0) &
                            + eps*(2.d0*varptb - 1.d0)
                    end do
                end do
            end do
        else if (pf_IC==3) then
            call read_phase_field_params(A, B, alpha)
            ! Scallop initial condition
            yh = [0.0, ylen/3, 2*ylen/3, ylen, &
                    ylen/6, ylen/2, 5*ylen/6, &
                    0.0, ylen/3, 2*ylen/3, ylen]
            zh(1:4) = 0.0
            zh(5:7) = zlen/2
            zh(8:11) = zlen
            x0 = 0.8
            gamma = 0.9
            do i=xstartr(3),xendr(3)
                do j=xstartr(2),xendr(2)
                    h0 = 0.0
                    do k=1,11
                        h0 = max(h0, x0 - gamma*((ymr(j) - yh(k))**2 + (zmr(i) - zh(k))**2))
                    end do
                    do k=1,nxmr
                        call random_number(varptb)
                        sal(k,j,i) = 1.0 - B*erfc((h0 - xmr(k))/10/pf_eps) &
                                        + eps*(2.d0*varptb - 1.d0)
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