!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: CreateInitialConditions.F90                    !
!    CONTAINS: subroutine CreateInitialConditions         !
!                                                         ! 
!    PURPOSE: Initialization routine. Sets initial        !
!     conditions for velocity and temperature             !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine CreateInitialConditions
    use param
    use local_arrays, only: vy,vx,temp,vz
    use decomp_2d, only: xstart,xend
    use mpih
    implicit none
    integer :: j,k,i,kmid, io
    real :: xxx,yyy,zzz,eps,varptb,amp
    real :: h0,t0,Lambda,r, x0, A, B, alpha
    real, dimension(11) :: yh, zh
    logical :: exists

    call random_seed()

    if ((xminusU /= 0) .or. (xplusU /= 0)) then
        !CJH: Base velocity profile if walls in motion
        do i=xstart(3),xend(3)
            do j=xstart(2),xend(2)
                do k=1,nxm
                    xxx = xm(k)
                    vy(k,j,i) = xminusU + (xplusU - xminusU)*xxx/alx3
                end do
            end do
        end do
    end if

    if (gAxis == 1) then
        if ((RayT < 0) .and. (RayS <0)) then
            !CJH: Stratified shear layer initial condition
            do i=xstart(3),xend(3)
                do j=xstart(2),xend(2)
                    do k=1,nxm
                        vy(k,j,i) = tanh(xm(k) - alx3/2.0)
                    end do
                end do
            end do
        else
            !CJH: RBC initial condition as used in AFiD 1.0
            eps = 0.01
            do i=xstart(3),xend(3)
                do j=xstart(2),xend(2)
                    do k=1,nxm
                        xxx = xc(k)
                        yyy = ym(j)
                        vx(k,j,i) = vx(k,j,i) - eps*xxx**2*(1.0 - xxx)**2*cos(2.0*pi*yyy/ylen)

                        xxx = xm(k)
                        yyy = yc(j)
                        vy(k,j,i) = vy(k,j,i) + 2.0*eps*xxx*(1.0 - xxx)*(1.0 - 2.0*xxx)*sin(2.0*pi*yyy/ylen)/pi
                    end do
                end do
            end do
        end if
    else if (gAxis == 2) then
        if (inslwN==1) then
        !CJH Laminar vertical convection as in Batchelor (1954)
            do i=xstart(3),xend(3)
                do j=xstart(2),xend(2)
                    do k=1,nxm
                        xxx = xm(k)
                        vy(k,j,i) = min(ren,320.0)/12.0*xxx*(2*xxx-1)*(xxx-1)
                    end do
                end do
            end do
        else
            !CJH Ke et al profile: 4t/(1-Pr)[i2erfc(eta) - i2erfc(eta/sqrt(Pr))]
            t0 = 1.4195567
            do i=xstart(3),xend(3)
                do j=xstart(2),xend(2)
                    do k=1,nxm
                        xxx = xm(k)/2*sqrt(pect/t0)
                        vy(k,j,i) = 4*t0/(1 - PraT)*(&
                            0.25*((1 + 2*xxx**2)*erfc(xxx) - 2*xxx/sqrt(pi)*exp(-xxx**2)) - &
                            0.25*((1 + 2*xxx**2/PraT)*erfc(xxx/sqrt(PraT)) &
                                - 2*xxx/sqrt(pi*PraT)*exp(-xxx**2/PraT)) &
                        )
                    end do
                end do
            end do
        end if
    end if

    if (dPdy > 0) then
        !CJH Interpret dPdy as Re_tau
        if (inslwn==1) then
            do i=xstart(3),xend(3)
                do j=xstart(2),xend(2)
                    do k=1,nxm
                        xxx = xm(k)
                        vy(k,j,i) = vy(k,j,i) + 4.0*dPdy**2/ren*xxx*(1.0 - xxx)
                    end do
                end do
            end do
        else
            eps = 0.1
            r = 2.0*pi*4/ylen
            do i=xstart(3),xend(3)
                do j=xstart(2),xend(2)
                    do k=1,nxm
                        xxx = xc(k)
                        yyy = ym(j)
                        vx(k,j,i) = - eps*xxx**2*(1.0 - xxx)**2*cos(r*yyy)

                        xxx = xm(k)
                        yyy = yc(j)
                        vy(k,j,i) = 20.0*xxx*(2.0 - xxx) + 2.0*eps*xxx*(1.0 - xxx)*(1.0 - 2.0*xxx)*sin(r*yyy)/r
                    end do
                end do
            end do
        end if
    end if

    if (dPdz > 0) then
        !CJH Interpret dPdz as Re_b
        do i=xstart(3),xend(3)
            do j=xstart(2),xend(2)
                do k=1,nxm
                    xxx = xm(k)
                    vz(k,j,i) = vz(k,j,i) + 6.0*dPdz/ren*xxx*(1.0 - xxx)
                end do
            end do
        end do
    end if

    ! Set velocity to zero if we are using the phase-field method
    if (melt .or. phasefield .or. IBM) then
        do i=xstart(3),xend(3)
            do j=xstart(2),xend(2)
                do k=1,nxm
                    vx(k,j,i) = 0.d0
                    vy(k,j,i) = 0.d0
                    vz(k,j,i) = 0.d0
                end do
            end do
        end do
    end if

    if ((RayT < 0) .and. (RayS < 0)) then
        !CJH: Stratified shear layer + noise in centre
        eps = 1e-2
        do i=xstart(3),xend(3)
            do j=xstart(2),xend(2)
                do k=1,nxm
                    temp(k,j,i) = tanh(xm(k) - 0.5*alx3)
                    call random_number(varptb)
                    temp(k,j,i) = temp(k,j,i) + &
                            cosh(xm(k) - 0.5*alx3)**(-2)*eps*(2.0*varptb - 1.0)
                end do
            end do
        end do
    else
        ! Assign linear temperature profile in the nodes k=1 to k=nxm
        do i=xstart(3),xend(3)
            do j=xstart(2),xend(2)
                do k=1,nxm
                    xxx = xm(k)
                    temp(k,j,i) = tempbp(1,j,i) + (temptp(1,j,i) - tempbp(1,j,i))*xm(k)/alx3
                end do
            end do
        end do

        ! Add noise in the temperature profile
        eps = 1e-3
        do i=xstart(3),xend(3)
            do j=xstart(2),xend(2)
                do k=1,nxm
                    call random_number(varptb)
                    if (abs(xm(k)-0.5) + eps > 0.5) then
                    amp = 0.5 - abs(xm(k)-0.5) ! CJH Prevent values of |T| exceeding 0.5
                    temp(k,j,i) = temp(k,j,i) + amp*(2.d0*varptb - 1.d0)
                    else
                    temp(k,j,i) = temp(k,j,i) + eps*(2.d0*varptb - 1.d0)
                    end if
                end do
            end do
        end do
    end if

    if (gAxis==3 .and. active_T==0) then
        do i=xstart(3),xend(3)
            do j=xstart(2),xend(2) ! Convergence test
                do k=1,nxm
                    xxx = xm(k) ! Linear profile + sin perturbation
                    temp(k,j,i) = tempbp(1,j,i) + (temptp(1,j,i) - tempbp(1,j,i))*xm(k)/alx3
                    temp(k,j,i) = temp(k,j,i) + sin(2.0*pi*xxx/alx3) - sin(6.0*pi*xxx/alx3)
                end do
            end do
        end do
    end if

    if (gAxis==2 .and. inslwN==0) then  ! Ke et al comparison case
        t0 = 1.4195567
        do i=xstart(3),xend(3)
            do j=xstart(2),xend(2)
                do k=1,nxm
                    amp = 0.0
                    do kmid=0,7
                        amp = amp + sin(2.0**kmid * 2.0*pi*ym(j)/ylen)
                    end do
                    amp = 1.0 + 1e-3*amp + 1e-3*sin(46.0*pi*zm(i)/zlen)
                    temp(k,j,i) = amp*erfc(xm(k)/2*sqrt(pect/t0))
                end do
            end do
        end do
    end if

    if (IBM .and. dPdy/=0) then
        do i=xstart(3),xend(3)
            do j=xstart(2),xend(2)
                do k=1,nxm
                    temp(k,j,i) = 0.0
                end do
            end do
        end do
    end if

    if (phasefield) then

        if (pf_IC == 1) then        ! 1D moving interface example
            h0 = 0.1                ! Freezing if RAYT < 0
            Lambda = 0.620063       ! Melting if RAYT > 0
            t0 = pect * (h0/2/Lambda)**2
            do i=xstart(3),xend(3)
                do j=xstart(2),xend(2)
                    do k=1,nxm
                        xxx = xm(k)
                        if (xxx < h0) then
                            temp(k,j,i) = erf(xxx*sqrt(pect/t0)/2)/erf(Lambda)
                        else
                            temp(k,j,i) = 1.0
                        end if
                        if (RAYT > 0) temp(k,j,i) = 1.0 - temp(k,j,i)
                    end do
                end do
            end do

        else if (pf_IC == 2) then
            do i=xstart(3),xend(3)
                do j=xstart(2),xend(2)
                    do k=1,nxm
                        r = sqrt((xm(k) - 0.5)**2 + (ym(j) - ylen/2)**2)
                        temp(k,j,i) = 0.5*(1.0 + tanh(100.0*(r - 0.1)))
                    end do
                end do
            end do

        else if (pf_IC == 3) then ! Favier et al (2019) Appendix A3 Validation Case
            eps = 0.1
            kmid = nxm/2
            do i=xstart(3),xend(3)
                do j=xstart(2),xend(2)
                    if (nzm > 1) then
                        do k=1,kmid ! If domain 3D, add in z perturbation too
                            xxx = xm(k)
                            yyy = ym(j)
                            zzz = zm(i)

                            temp(k,j,i) = temp(k,j,i) &
                                + eps*sin(4.0*pi*yyy)*cos(4.0*pi*zzz)*sin(2.0*pi*xxx)**2
                        end do
                    else
                        do k=1,kmid
                            xxx = xm(k)
                            yyy = ym(j)

                            temp(k,j,i) = temp(k,j,i) &
                                + eps*sin(4.0*pi*yyy)*sin(2.0*pi*xxx)**2
                        end do
                    end if
                end do
            end do

        elseif (pf_IC == 4) then ! Ice block to compare with Neufeld et al (2010)
            do i=xstart(3),xend(3)      ! Ice at bottom of domain if RayT<0, at the top if RayT>0
                do j=xstart(2),xend(2)
                    r = sqrt((ym(j) - ylen/2.0)**2 + (zm(i) - zlen/2.0)**2)
                    do k=1,nxm
                        temp(k,j,i) = 1.0 - 0.25*(1.0 + sign(1.0,RayT)*tanh((xm(k) - alx3/2.0)/2.0/pf_eps)) &
                                    *(1.0 - tanh((r - alx3/2.0)/2.0/pf_eps))
                    end do
                end do
            end do

        else if (pf_IC == 5) then        ! 1D supercooling example
            h0 = 0.02
            Lambda = 0.060314
            t0 = pect * (h0/2/Lambda)**2
            do i=xstart(3),xend(3)
                do j=xstart(2),xend(2)
                    do k=1,nxm
                        xxx = xm(k)
                        if (xxx > h0) then
                            temp(k,j,i) = erfc(xxx*sqrt(pect/t0)/2.0)/erfc(Lambda)
                        else
                            temp(k,j,i) = 1.0
                        end if
                        ! if (RAYT > 0) temp(k,j,i) = 1.0 - temp(k,j,i)
                    end do
                end do
            end do
        end if

        if (salinity) then
            if (pf_IC==1) then
                call read_phase_field_params(A, B, alpha)
                t0 = 1e-3
                x0 = 0.8
                h0 = x0 + 2*alpha*sqrt(t0)
                do i=xstart(3),xend(3)
                    do j=xstart(2),xend(2)
                        do k=1,nxm
                            if (xm(k) <= h0) then
                                temp(k,j,i) = 1 - A*erfc((x0 - xm(k))/sqrt(t0)/2.0)
                            else
                                temp(k,j,i) = 1 - A*erfc(-alpha)
                            end if
                        end do
                    end do
                end do
            else if (pf_IC==2) then
                call read_phase_field_params(A, B, alpha)
                t0 = 1e-3
                h0 = 0.1 - 2*alpha*sqrt(t0)
                eps = 5e-3
                do i=xstart(3),xend(3)
                    do j=xstart(2),xend(2)
                        do k=1,nxm
                            call random_number(varptb)
                            if (abs(ym(j) - ylen/2.0) <= h0) then
                                temp(k,j,i) = 1.0 - A*erfc(-alpha)
                            else if (ym(j) < ylen/2.0) then
                                temp(k,j,i) = 1.0 - A*erfc((ylen/2.0 - h0 - ym(j))/sqrt(t0)/2.0) &
                                                + eps*(2.d0*varptb - 1.d0)
                            else
                                temp(k,j,i) = 1.0 - A*erfc((ym(j) - ylen/2.0 - h0)/sqrt(t0)/2.0) &
                                + eps*(2.d0*varptb - 1.d0)
                            end if
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
                amp = 0.9
                eps = 5e-3
                do i=xstart(3),xend(3)
                    do j=xstart(2),xend(2)
                        h0 = 0.0
                        do k=1,11
                            h0 = max(h0, x0 - amp*((ym(j) - yh(k))**2 + (zm(i) - zh(k))**2))
                        end do
                        do k=1,nxm
                            call random_number(varptb)
                            if (xm(k) <= h0) then
                                temp(k,j,i) = 1.0 + eps*(2.d0*varptb - 1.d0)
                            else
                                temp(k,j,i) = 1.0 - A*erfc(-alpha)
                            end if
                        end do
                    end do
                end do
            else
                kmid = nxm/2
                do i=xstart(3),xend(3)
                    do j=xstart(2),xend(2)
                        do k=1,kmid
                            temp(k,j,i) = 1.0
                        end do
                        do k=kmid+1,nxm
                            temp(k,j,i) = 0.0
                        end do
                    end do
                end do
            end if
        end if

        if (IBM) then
            do i=xstart(3),xend(3)
                do j=xstart(2),xend(2)
                    h0 = 0.25 + (ym(j) - 0.5)**2
                    do k=1,nxm
                        if (xm(k) < h0) then
                            temp(k,j,i) = 0.0
                        else
                            temp(k,j,i) = 1.0
                        end if
                    end do
                end do
            end do
        end if

    end if

    if (melt) then
        A = 1.08995
        do i=xstart(3),xend(3)
            do j=xstart(2),xend(2)
                do k=1,nxm
                    ! call random_number(varptb)
                    ! temp(k,j,i) = eps*(2.d0*varptb - 1.d0) * exp(-xm(k)/0.1)
                    temp(k,j,i) = 1.0 - A*erfc(xm(k)*sqrt(pect)/2.0)
                end do
            end do
        end do
    end if

    ! FAVIER ET AL. (2019) APPENDIX A1 VALIDATION CASE
    ! if (phasefield) then
    !   ! kmid = nxm/2
    !   eps = 8.041
    !   do i=xstart(3),xend(3)
    !     do j=xstart(2),xend(2)
    !       do k=1,nxm
    !         temp(k,j,i) = (exp(-eps*(xm(k) - 1.0)) - 1.0)/(exp(eps) - 1.0)
    !       end do
    !       ! do k=1,kmid
    !       !   temp(k,j,i) = 1.0 - 2.0*xm(k)
    !       ! end do
    !       ! do k=kmid+1,nxm
    !       !   temp(k,j,i) = 0.1 - 0.2*xm(k)
    !       ! end do
    !     end do
    !   end do
    ! end if

    return
end subroutine CreateInitialConditions