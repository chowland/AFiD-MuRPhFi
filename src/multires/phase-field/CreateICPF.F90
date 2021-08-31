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
    real :: r

    if (pf_IC.eq.1) then ! 1D freezing validation
        do i=xstartr(3),xendr(3)
            do j=xstartr(2),xendr(2)
                do k=1,nxmr
                    ! if (xmr(k) < 0.1) then
                    !     phi(k,j,i) = 1.0
                    ! else
                    !     phi(k,j,i) = 0.0
                    ! end if
                    phi(k,j,i) = 0.5*(1.0 - tanh((xmr(k) - 0.1)/2/pf_eps))
                    if (RAYT > 0) phi(k,j,i) = 1.0 - phi(k,j,i)
                end do
            end do
        end do
    else if (pf_IC.eq.2) then ! Solid disc of radius 0.1 at x=0.5, centre of y-domain
        do i=xstartr(3),xendr(3)
            do j=xstartr(2),xendr(2)
                do k=1,nxmr
                    r = sqrt((xmr(k) - 0.5)**2 + (ymr(j) - ylen/2)**2)
                    phi(k,j,i) = 0.5*(1.0 - tanh((r - 0.1)/2/pf_eps))
                end do
            end do
        end do
    else if (pf_IC.eq.1) then ! Favier et al (2019) Appendix A3 Validation Case
        do i=xstartr(3),xendr(3)
            do j=xstartr(2),xendr(2)
                do k=1,nxmr
                    phi(k,j,i) = 0.5*(1.0 + tanh((xmr(k) - 0.5)/2/pf_eps))
                end do
            end do
        end do
    end if

    if (salinity) then ! Ice above salty water
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

    return

end subroutine CreateICPF