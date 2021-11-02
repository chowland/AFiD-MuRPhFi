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

    integer :: i,j,k,kmid
    real :: r, x0, lambda, h0, t0

    if (pf_IC == 1) then ! 1D freezing validation
        do i=xstartr(3),xendr(3)
            do j=xstartr(2),xendr(2)
                do k=1,nxmr
                    phi(k,j,i) = 0.5*(1.0 - tanh((xmr(k) - 0.1)/2/pf_eps))
                    if (RAYT > 0) phi(k,j,i) = 1.0 - phi(k,j,i)
                end do
            end do
        end do
    else if (pf_IC == 2) then ! Solid disc of radius 0.1 at x=0.5, centre of y-domain
        do i=xstartr(3),xendr(3)
            do j=xstartr(2),xendr(2)
                do k=1,nxmr
                    r = sqrt((xmr(k) - 0.5)**2 + (ymr(j) - ylen/2)**2)
                    phi(k,j,i) = 0.5*(1.0 - tanh((r - 0.1)/2/pf_eps))
                end do
            end do
        end do
    else if (pf_IC == 3) then ! Favier et al (2019) Appendix A3 Validation Case
        do i=xstartr(3),xendr(3)
            do j=xstartr(2),xendr(2)
                do k=1,nxmr
                    phi(k,j,i) = 0.5*(1.0 + tanh((xmr(k) - 0.5)/2/pf_eps))
                end do
            end do
        end do
    else if (pf_IC == 4) then ! Ice block to compare with Neufeld et al. (2010)
        do i=xstartr(3),xendr(3)
            do j=xstartr(2),xendr(2)
                r = sqrt((ymr(j) - ylen/2.0)**2 + (zmr(i) - zlen/2.0)**2)
                do k=1,nxmr
                    phi(k,j,i) = 0.25*(1.0 + sign(1.0,RayT)*tanh((xmr(k) - alx3/2.0)/2.0/pf_eps)) &
                                    *(1.0 - tanh((r - alx3/2.0)/2.0/pf_eps))
                end do
            end do
        end do
    else if (pf_IC == 5) then ! 1D freezing supercooled validation
        do i=xstartr(3),xendr(3)
            do j=xstartr(2),xendr(2)
                do k=1,nxmr
                    phi(k,j,i) = 0.5*(1.0 - tanh((xmr(k) - 0.02)/2/pf_eps))
                end do
            end do
        end do
    end if

    if (salinity) then ! Ice above salty water
        if (pf_IC==1) then
            t0 = 1.0
            x0 = 0.8
            lambda = 0.2080 !0.24041
            h0 = x0 + 2*Lambda*sqrt(t0/pect)
            do i=xstartr(3),xendr(3)
                do j=xstartr(2),xendr(2)
                    do k=1,nxmr
                        phi(k,j,i) = 0.5*(1.0 + tanh((xmr(k) - h0)/2/pf_eps))
                    end do
                end do
            end do
        else
            kmid = nxmr/2
            do i=xstartr(3),xendr(3)
                do j=xstartr(2),xendr(2)
                    do k=1,kmid
                        phi(k,j,i) = 0.0
                    end do
                    do k=kmid+1,nxmr
                        phi(k,j,i) = 1.0
                    end do
                end do
            end do
        end if
    end if

    return

end subroutine CreateICPF