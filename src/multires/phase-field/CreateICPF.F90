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
    real :: r, x0, h0, t0, A, B, alpha
    real, dimension(11) :: yh, zh

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
            call read_phase_field_params(A, B, alpha)
            t0 = 1e-3
            x0 = 0.8
            h0 = x0 + 2*alpha*sqrt(t0)
            do i=xstartr(3),xendr(3)
                do j=xstartr(2),xendr(2)
                    do k=1,nxmr
                        phi(k,j,i) = 0.5*(1.0 + tanh((xmr(k) - h0)/2/pf_eps))
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
                        phi(k,j,i) = 0.5*(tanh((ymr(j) - (ylen/2.0 - h0))/2.0/pf_eps) - &
                                        tanh((ymr(j) - (ylen/2.0 + h0))/2.0/pf_eps))
                    end do
                end do
            end do
        else if (pf_IC==3) then
            ! Scallop initial condition
            yh = [0.0, ylen/3, 2*ylen/3, ylen, &
                    ylen/6, ylen/2, 5*ylen/6, &
                    0.0, ylen/3, 2*ylen/3, ylen]
            zh(1:4) = 0.0
            zh(5:7) = zlen/2
            zh(8:11) = zlen
            x0 = 0.8
            A = 0.9
            do i=xstartr(3),xendr(3)
                do j=xstartr(2),xendr(2)
                    h0 = 0.0
                    do k=1,11
                        h0 = max(h0, x0 - A*((ymr(j) - yh(k))**2 + (zmr(i) - zh(k))**2))
                    end do
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

    if (IBM) then
        do i=xstartr(3),xendr(3)
            do j=xstartr(2),xendr(2)
                h0 = 0.25 + (ymr(j) - 0.5)**2
                do k=1,nxmr
                    phi(k,j,i) = 0.5*(1.0 - tanh((xmr(k) - h0)/2/pf_eps))
                end do
            end do
        end do
    end if

    return

end subroutine CreateICPF

subroutine read_phase_field_params(A, B, alpha)
    implicit none
    real, intent(out) :: A, B, alpha

    integer :: io
    logical :: exists

    inquire(file="pfparam.in", exist=exists)
    if (exists) then
        open(newunit=io, file="pfparam.in", status="old", action="read")
        read(io, *) A, B, alpha
        close(io)
    else
        A = 1.132
        B = 0.3796
        alpha = 3.987e-2
    end if
end subroutine read_phase_field_params