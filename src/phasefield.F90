!> Module adding phase-field model to AFiD
!! to simulate the melting of solid objects
module afid_phasefield
    use param
    use mgrd_arrays
    use decomp_2d, only: xstart, xend, xstartr, xendr, update_halo
    use AuxiliaryRoutines
    use HermiteInterpolations, only: interpolate_xyz_to_coarse, interpolate_xyz_to_coarse_fast, interpolate_xyz_to_refined
    use ibm_param, only: solidr
    implicit none

    real, allocatable, dimension(:,:,:) :: phi      !! Phase-field variable
    real, allocatable, dimension(:,:,:) :: ruphi    !! RK storage array for phase-field (previous substep)
    real, allocatable, dimension(:,:,:) :: hphi     !! RK storage array for phase-field
    real, allocatable, dimension(:,:,:) :: phic     !! Interpolated phase-field on coarse grid (also used to store d(phi)/dt)
    real, allocatable, dimension(:,:,:) :: tempr    !! Interpolated temperature field on refined grid

    real :: pf_A        !! Phase-field Gibbs-Thomson parameter
    real :: pf_eps      !! Phase-field interface thickness
    real :: pf_D        !! Dimensionless diffusivity of the phase-field
    real :: pf_S        !! Stefan number
    real :: pf_Tm       !! Equilibrium melting temperature
    real :: pf_Lambda   !! Dimensionless liquidus slope

    real, allocatable, dimension(:) :: ap3spkr      !! Upper diagonal derivative coefficient for salinity
    real, allocatable, dimension(:) :: ac3spkr      !! Diagonal derivative coefficient for salinity
    real, allocatable, dimension(:) :: am3spkr      !! Lower diagonal derivative coefficient for salinity

    real :: pf_direct_force = 0.9      !! phi contour above which velocity is forced exactly to zero pre-pressure solve

contains

!> Initialise and allocate memory for phase-field variables
subroutine InitPFVariables

    ! Main array with ghost cells
    call AllocateReal3DArray(phi,1,nxr,xstartr(2)-lvlhalo,xendr(2)+lvlhalo,xstartr(3)-lvlhalo,xendr(3)+lvlhalo)
    ! Refined temperature array
    call AllocateReal3DArray(tempr,1,nxr,xstartr(2)-lvlhalo,xendr(2)+lvlhalo,xstartr(3)-lvlhalo,xendr(3)+lvlhalo)

    ! Arrays without ghost cells
    call AllocateReal3DArray(ruphi,1,nxr,xstartr(2),xendr(2),xstartr(3),xendr(3))
    call AllocateReal3DArray(hphi,1,nxr,xstartr(2),xendr(2),xstartr(3),xendr(3))

    ! Coarse array for phi or d(phi)/dt
    call AllocateReal3DArray(phic,1,nx,xstart(2)-lvlhalo,xend(2)+lvlhalo,xstart(3)-lvlhalo,xend(3)+lvlhalo)

    ! Second derivative coefficients
    call AllocateReal1DArray(ap3spkr,1,nxr)
    call AllocateReal1DArray(ac3spkr,1,nxr)
    call AllocateReal1DArray(am3spkr,1,nxr)

end subroutine InitPFVariables

!> Deallocate the variables used for evolving the phase-field
subroutine DeallocatePFVariables

    ! Main array
    call DestroyReal3DArray(phi)

    ! Array for refined temperature
    call DestroyReal3DArray(tempr)

    ! Arrays without ghost cells
    call DestroyReal3DArray(ruphi)
    call DestroyReal3DArray(hphi)

    ! Coarse array for phi or d(phi)/dt
    call DestroyReal3DArray(phic)

    ! Second derivative coefficients
    call DestroyReal1DArray(ap3spkr)
    call DestroyReal1DArray(ac3spkr)
    call DestroyReal1DArray(am3spkr)

end subroutine DeallocatePFVariables

!> Set the initial state of the phase field
subroutine CreateInitialPhase
    ! use afid_salinity, only: sal, PraS
    real :: h0

    if (salinity) then
        !! Ice above salty water (1D_DDMelting example)
        call set_multicomponent_interface(0.8, h0)
        call set_flat_interface(h0, .true.)
    else
        !! 1D freezing/moving example
        !! (RayT > 0: melting; RayT < 0: freezing)
        if (pf_IC==1) then
            h0 = 0.1
            if (RayT > 0) then
                call set_flat_interface(h0, .true.)
            else
                call set_flat_interface(h0, .false.)
            end if
            call set_temperature_interface(h0, .false.)

        !! 2D axisymmetric melting example (disc radius 0.1)
        else if (pf_IC==2) then
            call set_ice_disc(r0=0.1)

        !! Favier (2019) appendix A.3 validation case
        elseif (pf_IC==3) then
            h0 = 0.4
            call set_flat_interface(h0, .true.)
            call add_temperature_mode(amp=1e-2, ymode=10, zmode=2, h0=h0)

        else if (pf_IC==4) then
            call set_ice_sphere(r0=0.1)
        
        !! 1D supercooling example
        elseif (pf_IC==5) then
            h0 = 0.02
            call set_flat_interface(h0, .false.)
            call set_temperature_interface(h0, .true.)
        end if
    end if
    
end subroutine CreateInitialPhase

!> Read a simple input file with three parameters definining
!! the initial state for multicomponent melting
subroutine read_phase_field_params(A, B, alpha)
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

!> Impose a flat interface at height x=h0
!! with ice phase for x>h0 if ice_above=.true. and vice versa
subroutine set_flat_interface(h0, ice_above)
    real, intent(in) :: h0              !! height (in x) of interface
    logical, intent(in) :: ice_above    !! flag to determine whether ice phase is above or below interface

    integer :: i, j, k

    do i=xstartr(3),xendr(3)
        do j=xstartr(2),xendr(2)
            do k=1,nxmr
                if (ice_above) then
                    phi(k,j,i) = 0.5*(1.0 + tanh((xmr(k) - h0)/2/pf_eps))
                else
                    phi(k,j,i) = 0.5*(1.0 - tanh((xmr(k) - h0)/2/pf_eps))
                end if
            end do
        end do
    end do
end subroutine set_flat_interface

!> Impose a diffusive boundary layer profile for temperature
!! on one side of the flat interface at height x=h0
subroutine set_temperature_interface(h0, diffuse_above)
    use local_arrays, only: temp
    real, intent(in) :: h0                  !! Interface position (at sim time t=0)
    logical, intent(in) :: diffuse_above    !! Flag: set diffusion above interface (or not)

    real :: t0, Lambda
    integer :: i, j, k

    ! Set normalising similarity value (see docs for justification)
    if (diffuse_above) then
        Lambda = 0.060314
    else
        Lambda = 0.620063
    end if

    ! Effective time assuming that initial interface diffused from
    ! x=0 at t=0
    t0 = PecT*(h0/2.0/Lambda)**2

    do i=xstart(3),xend(3)
        do j=xstart(2),xend(2)
            do k=1,nxm
                if (diffuse_above) then
                    !! For the 1D supercooling example
                    if (xm(k) > h0) then
                        temp(k,j,i) = erfc(xm(k)*sqrt(pect/t0)/2.0)/erfc(Lambda)
                    else
                        temp(k,j,i) = 1.0
                    end if
                else
                    !! For the 1D freezing example
                    if (xm(k) < h0) then
                        temp(k,j,i) = erf(xm(k)*sqrt(pect/t0)/2)/erf(Lambda)
                    else
                        temp(k,j,i) = 1.0
                    end if
                    !! For the 1D melting case
                    if (RayT > 0) temp(k,j,i) = 1.0 - temp(k,j,i)
                end if
            end do
        end do
    end do
end subroutine set_temperature_interface

!> Set temperature and salinity profiles to a diffusive boundary layer below
!! an ice interface that began as a step profile at height x0
subroutine set_multicomponent_interface(x0, h0)
    use local_arrays, only: temp
    use afid_salinity, only: sal, PraS
    real, intent(in) :: x0  !! initial interface height (at diffusion start)
    real, intent(out) :: h0 !! interface height at simulation time 0
    
    real :: t0, A, B, alpha
    integer :: i, j, k
    
    call read_phase_field_params(A, B, alpha)
    t0 = 1e-3
    h0 = x0 + 2*alpha*sqrt(t0)
    do i=xstartr(3),xendr(3)
        do j=xstartr(2),xendr(2)
            do k=1,nxmr
                sal(k,j,i) = 1.0 - B*erfc((x0 - xmr(k))/sqrt(PraT/PraS*t0)/2.0)
            end do
        end do
    end do
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
end subroutine set_multicomponent_interface

!> Add a modal perturbation to the temperature field in the lower half of the domain
!! Intended for use with phase-field validation case of 2D Melting RBC from Favier et al (2019)
!! Vertical structure is of the form sin^2(pi x/h0)
subroutine add_temperature_mode(amp, ymode, zmode, h0)
    use local_arrays, only: temp
    real, intent(in) :: amp         !! amplitude of perturbation
    integer, intent(in) :: ymode    !! mode number in y
    integer, intent(in) :: zmode    !! mode number in z (only used if nzm>1)
    real, intent(in) :: h0          !! initial interface height
    integer :: i, j, k
    real :: xxx, yyy, zzz

    do i=xstart(3),xend(3)
        zzz = zm(i)/zlen
        do j=xstart(2),xend(2)
            yyy = ym(j)/ylen
            if (nzm > 1) then
                do k=1,nxm ! If domain 3D, add in z perturbation too
                    xxx = xm(k)
                    if (xxx < h0) then
                        temp(k,j,i) = temp(k,j,i) &
                            + amp*sin(2.0*pi*ymode*yyy)*cos(2.0*pi*zmode*zzz)*sin(pi*xxx/h0)**2
                    end if
                end do
            else
                do k=1,nxm
                    xxx = xm(k)
                    if (xxx < h0) then
                        temp(k,j,i) = temp(k,j,i) &
                            + amp*sin(2.0*pi*ymode*yyy)*sin(pi*xxx/h0)**2
                    end if
                end do
            end if
        end do
    end do
end subroutine add_temperature_mode

!> Add an ice disc at the centre of the domain of radius r0
!! with a corresponding temperature profile
subroutine set_ice_disc(r0)
    use local_arrays, only: temp
    real, intent(in) :: r0      !! Initial disc radius

    real :: r
    integer :: i, j, k

    ! Temperature field (0 in disc, 1 out of disc, tanh interface width approx 1e-2)
    do i=xstart(3),xend(3)
        do j=xstart(2),xend(2)
            do k=1,nxm
                r = sqrt((xm(k) - 0.5*alx3)**2 + (ym(j) - 0.5*ylen)**2)
                temp(k,j,i) = 0.5*(1.0 + tanh(100.0*(r - r0)))
            end do
        end do
    end do
    ! Phase-field (0 out of disc, 1 in disc)
    do i=xstartr(3),xendr(3)
        do j=xstartr(2),xendr(2)
            do k=1,nxmr
                r = sqrt((xmr(k) - 0.5*alx3)**2 + (ymr(j) - 0.5*ylen)**2)
                phi(k,j,i) = 0.5*(1.0 - tanh(0.5*(r - r0)/pf_eps))
            end do
        end do
    end do
end subroutine set_ice_disc

!> Add an ice sphere at the centre of the domain of radius r0
!! with a corresponding temperature profile
subroutine set_ice_sphere(r0)
    use local_arrays, only: temp
    real, intent(in) :: r0      !! Initial disc radius

    real :: r
    integer :: i, j, k

    ! Temperature field (0 in disc, 1 out of disc, tanh interface width approx 1e-2)
    do i=xstart(3),xend(3)
        do j=xstart(2),xend(2)
            do k=1,nxm
                r = sqrt((xm(k) - 0.5*alx3)**2 + (ym(j) - 0.5*ylen)**2 + (zm(i) - 0.5*zlen)**2)
                ! temp(k,j,i) = 0.5*(1.0 + tanh(100.0*(r - r0)))
                if (r > r0) then
                    temp(k,j,i) = 1.0
                else
                    temp(k,j,i) = 0.0
                end if
            end do
        end do
    end do
    ! Phase-field (0 out of disc, 1 in disc)
    do i=xstartr(3),xendr(3)
        do j=xstartr(2),xendr(2)
            do k=1,nxmr
                r = sqrt((xmr(k) - 0.5*alx3)**2 + (ymr(j) - 0.5*ylen)**2 + (zmr(i) - 0.5*zlen)**2)
                phi(k,j,i) = 0.5*(1.0 - tanh(0.5*(r - r0)/pf_eps))
            end do
        end do
    end do
end subroutine set_ice_sphere

!> Compute the explicit terms for the phase-field evolution
!! and store the result in `hphi`
subroutine ExplicitPhase
    integer :: ic, jc, kc
    integer :: im, jm
    integer :: ip, jp

    real :: udyrq, udzrq
    real :: pf_B, nlphi
    real :: dyyp, dzzp

    ! Nonlinear term prefactor
    pf_B = pf_D/(pf_eps)**2

    ! Diffusion coefficients
    udzrq = pf_D*dzqr
    udyrq = pf_D*dyqr

    do ic=xstartr(3),xendr(3)
        im = ic - 1
        ip = ic + 1
        do jc=xstartr(2),xendr(2)
            jm = jc - 1
            jp = jc + 1
            do kc=1,nxmr
                ! yy second derivative of phi
                dyyp = (phi(kc,jp,ic) - 2.0*phi(kc,jc,ic) + phi(kc,jm,ic))*udyrq
                ! zz second derivative of phi
                dzzp = (phi(kc,jc,ip) - 2.0*phi(kc,jc,ic) + phi(kc,jc,im))*udzrq
                ! Extra nonlinear terms
                nlphi = pf_B*phi(kc,jc,ic)*(1.0 - phi(kc,jc,ic)) &
                        *(1.0 - 2.0*phi(kc,jc,ic) + pf_A*(tempr(kc,jc,ic) - pf_Tm))

                hphi(kc,jc,ic) = dyyp + dzzp - nlphi
            end do
        end do
    end do

end subroutine ExplicitPhase

!> Compute the implicit terms for the phase-field evolution
subroutine ImplicitPhase
    integer :: jc,kc,ic
    real    :: alpec,dxxp

    alpec=al*pf_D

    do ic=xstartr(3),xendr(3)
        do jc=xstartr(2),xendr(2)
            do kc=1,nxmr

                ! Second xx derivative
                ! Apply lower BC (d/dx(phi)=0)
                if (kc.eq.1) then
                    dxxp = phi(kc+1,jc, ic)*ap3spkr(kc) &
                         + phi(kc  ,jc, ic)*ac3spkr(kc)
                ! Apply upper BC (d/dx(phi)=0)
                elseif(kc.eq.nxmr) then
                    dxxp = phi(kc  ,jc,ic)*ac3spkr(kc) &
                         + phi(kc-1,jc,ic)*am3spkr(kc)
                else
                    dxxp = phi(kc+1,jc,ic)*ap3spkr(kc) &
                         + phi(kc  ,jc,ic)*ac3spkr(kc) &
                         + phi(kc-1,jc,ic)*am3spkr(kc)
                end if

                rhsr(kc,jc,ic) = (ga*hphi(kc,jc,ic) + ro*ruphi(kc,jc,ic) + alpec*dxxp)*dt

                ruphi(kc,jc,ic) = hphi(kc,jc,ic)

            end do
        end do
    end do

!  Solve equation and update salinity
    call SolveImpEqnUpdate_Phi

end subroutine ImplicitPhase

!> Solve the implicit system for the phase-field
!! and update the global variable phi
subroutine SolveImpEqnUpdate_Phi
    integer :: jc,kc,info,ipkv(nxr),ic,nrhs
    real :: betadx,ackl_b
    real :: amkT(nxmr-1),ackT(nxmr),apkT(nxmr-1),appk(nxmr-2)

    betadx=0.5d0*al*dt*pf_D

    ! Construct tridiagonal matrix for LHS
    ! (normalised to prevent floating point errors)
    do kc=1,nxmr
        ackl_b=1.0d0/(1.-ac3spkr(kc)*betadx)
        if (kc > 1) amkT(kc-1) = -am3spkr(kc)*betadx*ackl_b
        ackT(kc)=1.0d0
        if (kc < nxmr) apkT(kc) = -ap3spkr(kc)*betadx*ackl_b
    end do

    ! Factor the tridiagonal matrix
    call dgttrf(nxmr,amkT,ackT,apkT,appk,ipkv,info)
    
    ! Rescale RHS to match rescaling of LHS
    nrhs=(xendr(3)-xstartr(3)+1)*(xendr(2)-xstartr(2)+1)
    do ic=xstartr(3),xendr(3)
        do jc=xstartr(2),xendr(2)
            do kc=1,nxmr
                ackl_b=1.0/(1.0-ac3spkr(kc)*betadx)
                rhsr(kc,jc,ic)=rhsr(kc,jc,ic)*ackl_b
            end do
        end do
    end do

    ! Solve tridiagonal system
    call dgttrs('N',nxmr,nrhs,amkT,ackT,apkT,appk,ipkv,rhsr,nxmr,info)

    ! Update global variable
    do ic=xstartr(3),xendr(3)
        do jc=xstartr(2),xendr(2)
            do kc=1,nxmr
                phi(kc,jc,ic) = phi(kc,jc,ic) + rhsr(kc,jc,ic)
            end do
        end do
    end do

end subroutine SolveImpEqnUpdate_Phi

!> Interpolate the phase-field onto the coarse grid to provide
!! volume penalty forcing for the momentum equation
subroutine InterpPhiMultigrid
    integer  :: ic,jc,kc

    ! Set coarse phase-field array to zero
    phic(:,:,:) = 0.d0

    ! Construct temporary array with extended range for interpolation
    ! (using dphi/dx = 0 BC)
    tpdvr(:,:,:) = 0.d0
    do ic=xstartr(3)-lvlhalo,xendr(3)+lvlhalo
        do jc=xstartr(2)-lvlhalo,xendr(2)+lvlhalo
            tpdvr(0,jc,ic) = phi(1,jc,ic)
            do kc=1,nxmr
                tpdvr(kc,jc,ic) = phi(kc,jc,ic)
            end do
            tpdvr(nxr,jc,ic) = phi(nxmr,jc,ic)
        end do
    end do

    ! Interpolate the refined field to the coarse grid, storing in phic
    if ((xmr(1) < xm(1)) .and. (xmr(nxmr) > xm(nxm))) then
        call interpolate_xyz_to_coarse_fast(tpdvr, phic(1:nxm,:,:), "phi")
    else
        call interpolate_xyz_to_coarse(tpdvr, phic(1:nxm,:,:))
    end if

end subroutine InterpPhiMultigrid

!> Interpolate the temperature field onto the refined grid
!! for calculation of nonlinear terms in phase-field equation
subroutine InterpTempMultigrid
    use local_arrays, only: temp, tempbp, temptp
    integer :: ic,jc,kc

    tempr(:,:,:) = 0.d0
    tpdv(:,:,:) = 0.d0

    ! Fill temporary array with temperature field and BCs
    do ic=xstart(3)-lvlhalo,xend(3)+lvlhalo
        do jc=xstart(2)-lvlhalo,xend(2)+lvlhalo
            tpdv(0,jc,ic) = 2.0*tempbp(1,jc,ic) - temp(1,jc,ic)
            tpdv(nx,jc,ic) = 2.0*temptp(1,jc,ic) - temp(nxm,jc,ic)
            do kc=1,nxm
                tpdv(kc,jc,ic) = temp(kc,jc,ic)
            end do
        end do
    end do

    ! Interpolate to the refined grid, storing in tempr
    call interpolate_xyz_to_refined(tpdv, tempr(1:nxmr,:,:))

end subroutine InterpTempMultigrid

!> Add volume penalty forcing to each component of the momentum
!! equation by modifying the RK explicit terms
subroutine AddVolumePenalty
    use local_arrays, only: qcap, dph, dq, vx, vy, vz
    real :: pf_eta, heta, volpen
    integer :: ic, jc, kc, jmm, imm

    ! Set optimal damping factor following Hester et al (2020)
    pf_eta = ren*(1.51044285*pf_eps)**2
    heta = 0.5d0/pf_eta

    ! x-component
    do ic=xstart(3),xend(3)
        do jc=xstart(2),xend(2)
            do kc=2,nxm
                volpen = (phic(kc,jc,ic) + phic(kc-1,jc,ic))*vx(kc,jc,ic)*heta
                qcap(kc,jc,ic) = qcap(kc,jc,ic) - volpen
            end do
        end do
    end do

    ! y-component
    do ic=xstart(3),xend(3)
        do jc=xstart(2),xend(2)
            jmm = jc - 1
            do kc=1,nxm
                volpen = (phic(kc,jc,ic) + phic(kc,jmm,ic))*vy(kc,jc,ic)*heta
                dph(kc,jc,ic) = dph(kc,jc,ic) - volpen
            end do
        end do
    end do

    ! z-component
    do ic=xstart(3),xend(3)
        imm = ic - 1
        do jc=xstart(2),xend(2)
            do kc=1,nxm
                volpen = (phic(kc,jc,ic) + phic(kc,jc,imm))*vz(kc,jc,ic)*heta
                dq(kc,jc,ic) = dq(kc,jc,ic) - volpen
            end do
        end do
    end do

end subroutine AddVolumePenalty

!> Add the normal interfacial salt flux term to the
!! explicit terms of the salinity evolution when using
!! the phase-field method
!! -1/Pe_S (grad(phi).grad(S)/(1 - phi + delta))
subroutine AddSaltFluxInterface
    use afid_salinity
    integer :: ic, jc, kc
    integer :: im, jm, km
    integer :: ip, jp, kp
    real :: sdx(nxmr), sqdx(nxmr), sqdy, sqdz, kaps
    real :: dpdsx, dpdsy, dpdsz, phfrac
    real :: pf_delta = 1e-6

    do kc=1,nxmr
        sdx(kc) = 0.5*dxr/g3rmr(kc)
        sqdx(kc) = sdx(kc)**2
    end do
    sqdy = (0.5d0*dyr)**2
    sqdz = (0.5d0*dzr)**2
    kaps = 1.d0/pecs

    do ic=xstartr(3),xendr(3)
        im = ic - 1
        ip = ic + 1
        do jc=xstartr(2),xendr(2)
            jm = jc - 1
            jp = jc + 1
            do kc=1,nxmr
                km = kc - 1
                kp = kc + 1

                if (kc==1) then
                    dpdsx = sqdx(kc)*(phi(kp,jc,ic) - phi(kc,jc,ic))&
                                    *(sal(kp,jc,ic) - sal(kc,jc,ic) &
                                + 2.0*SfixS*(sal(kc,jc,ic) - salbp(1,jc,ic)))
                elseif (kc==nxmr) then
                    dpdsx = sqdx(kc)*(phi(kc,jc,ic) - phi(km,jc,ic)) &
                                    *(sal(kc,jc,ic) - sal(km,jc,ic) &
                                + 2.0*SfixN*(saltp(1,jc,ic) - sal(kc,jc,ic)))
                else
                    dpdsx = sqdx(kc)*(phi(kp,jc,ic) - phi(km,jc,ic)) &
                                    *(sal(kp,jc,ic) - sal(km,jc,ic))
                end if
                dpdsy = sqdy*(phi(kc,jp,ic) - phi(kc,jm,ic))&
                            *(sal(kc,jp,ic) - sal(kc,jm,ic))
                dpdsz = sqdz*(phi(kc,jc,ip) - phi(kc,jc,im))&
                            *(sal(kc,jc,ip) - sal(kc,jc,im))

                phfrac = 1.0/(1.0 - phi(kc,jc,ic) + pf_delta)

                hsal(kc,jc,ic) = hsal(kc,jc,ic) - phfrac*kaps*(dpdsx + dpdsy + dpdsz)
            end do
        end do
    end do
end subroutine AddSaltFluxInterface

!> Add forcing terms to rhs of phase-field equation when
!! using salinity to reflect change in local melting temperature
!! with salt concentration
!! -D/eps^2 phi(1 - phi) A\Lambda S
subroutine AdjustMeltPoint
    use afid_salinity, only: sal
    integer :: ic, jc, kc
    real :: bcl

    ! Nonlinear term prefactor
    bcl = pf_D/(pf_eps)**2*pf_A*pf_Lambda

    do ic=xstartr(3),xendr(3)
        do jc=xstartr(2),xendr(2)
            do kc=1,nxmr
                hphi(kc,jc,ic) = hphi(kc,jc,ic) &
                                - bcl*sal(kc,jc,ic)*phi(kc,jc,ic)*(1.0 - phi(kc,jc,ic))
            end do
        end do
    end do
end subroutine AdjustMeltPoint

!> Use the result from the implicit solution of phi to compute a latent heat
!! term, interpolate it onto the coarse grid, and add it to the temperature
!! equation. Interpolation is explicitly written out here, rather than using a
!! function since we are adding to the RK array hro
!! When this routine is called, rhsr should contain phi(l+1)-phi(l) = d/dt(phi)*dt
subroutine AddLatentHeat
    use local_arrays, only: hro
    integer :: ic,jc,kc,icr,jcr,kcr,n

    real, dimension(4,4,4) :: qv3
    real, dimension(4,4) :: qv2
    real, dimension(4) :: qv1

    real :: phi_rhs, aldt

    phi_rhs = 0.d0
    aldt = 1.0/al/dt

    tpdvr(:,:,:) = 0.d0 ! Temporary array with extended range for interpolation

    ! Fill temporary array with dphi/dt
    do ic=xstartr(3),xendr(3)
        do jc=xstartr(2),xendr(2)
            !CJH Note if d/dx(phi)==0 on boundary, then
            ! d/dx(dphi/dt)=0 on boundary
            tpdvr(0,jc,ic) = rhsr(1,jc,ic)
            do kc=1,nxmr
                tpdvr(kc,jc,ic) = rhsr(kc,jc,ic)
            end do
            tpdvr(nxr,jc,ic) = rhsr(nxmr,jc,ic)
        end do
    end do

    ! Fill in halo values
    call update_halo(tpdvr,lvlhalo)
    if (sidewall) then
        if (xstart(2)==1) then
            do n=1,lvlhalo
                tpdvr(:,1-n,:) = tpdvr(:,n,:)
            end do
        end if
        if (xend(2)==nym) then
            do n=1,lvlhalo
                tpdvr(:,nym+n,:) = tpdvr(:,ny-n,:)
            end do
        end if
        if (xstart(3)==1) then
            do n=1,lvlhalo
                tpdvr(:,:,1-n) = tpdvr(:,:,n)
            end do
        end if
        if (xend(3)==nzm) then
            do n=1,lvlhalo
                tpdvr(:,:,nzm+n) = tpdvr(:,:,nz-n)
            end do
        end if
    end if

    ! If the grid is sufficiently fine at the boundary, use a more efficient interpolation
    if ((xmr(1) < xm(1)) .and. (xmr(nxmr) > xm(nxm))) then
        do ic=xstart(3),xend(3)
            icr = krangs(ic)
            do jc=xstart(2),xend(2)
                jcr = jrangs(jc)
                do kc=1,nxm
                    kcr = irangs(kc)

                    qv3 = tpdvr(kcr-2:kcr+1,jcr-2:jcr+1,icr-2:icr+1)
                    qv2(:,:) = qv3(:,:,1)*czsalc(1,ic) + qv3(:,:,2)*czsalc(2,ic) &
                                + qv3(:,:,3)*czsalc(3,ic) + qv3(:,:,4)*czsalc(4,ic)
                    qv1(:) = qv2(:,1)*cysalc(1,jc) + qv2(:,2)*cysalc(2,jc) &
                            + qv2(:,3)*cysalc(3,jc) + qv2(:,4)*cysalc(4,jc)
                        
                    phi_rhs = sum(qv1(1:4)*cxsalc(1:4,kc))
                    hro(kc,jc,ic) = hro(kc,jc,ic) + pf_S*phi_rhs*aldt
                end do
            end do
        end do
    else
        do icr=xstartr(3)-1,xendr(3)
            do jcr=xstartr(2)-1,xendr(2)
                do kcr=0,nxmr
            
                    qv3 = tpdvr(kcr-1:kcr+2,jcr-1:jcr+2,icr-1:icr+2)
                    do ic=max(krangr(icr),xstart(3)),min(krangr(icr+1)-1,xend(3))
                        qv2(:,:) = qv3(:,:,1)*czsalc(1,ic) + qv3(:,:,2)*czsalc(2,ic)&
                                +qv3(:,:,3)*czsalc(3,ic) + qv3(:,:,4)*czsalc(4,ic)
                        do jc=max(jrangr(jcr),xstart(2)),min(jrangr(jcr+1)-1,xend(2))
                            qv1(:) = qv2(:,1)*cysalc(1,jc) + qv2(:,2)*cysalc(2,jc) &
                                    +qv2(:,3)*cysalc(3,jc) + qv2(:,4)*cysalc(4,jc)
                            do kc=max(irangr(kcr),1),min(irangr(kcr+1)-1,nxm)
                                phi_rhs = sum(qv1(1:4)*cxsalc(1:4,kc))

                                hro(kc,jc,ic) = hro(kc,jc,ic) + pf_S*phi_rhs*aldt
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end if
end subroutine AddLatentHeat

!> Add "latent salt" term to the RK forcing array for salinity (hsal),
!! having calculated d/dt(phi) from the implicit solve and stored it in rhsr
subroutine AddLatentSalt
    use afid_salinity, only: sal, hsal
    real :: aldt, phfrac
    real :: pf_delta = 1e-6
    integer :: ic, jc, kc

    aldt = 1.0/al/dt

    do ic=xstartr(3),xendr(3)
        do jc=xstartr(2),xendr(2)
            do kc=1,nxmr
                phfrac = 1.0/(1.0 - phi(kc,jc,ic) + pf_delta)
                hsal(kc,jc,ic) = hsal(kc,jc,ic) &
                                + phfrac*sal(kc,jc,ic)*rhsr(kc,jc,ic)*aldt
            end do
        end do
    end do
end subroutine AddLatentSalt

!> Calculate vertical profiles of phase-field-related statistics and
!! store them in means.h5
subroutine CalcPhiStats

    real, dimension(nxmr) :: phibar    !! Horizontally-averaged phase-field variable
    real, dimension(nxmr) :: phirms    !! Horizontally-averaged rms phase-field

    real :: inyzm   !! 1.0/nymr/nzmr

    character(30) :: dsetname   !! Dataset name for HDF5 file
    character(30) :: filename   !! HDF5 file name for statistic storage
    character( 5) :: nstat      !! Character string of statistic index

    integer :: i, j, k

    inyzm = 1.0/nymr/nzmr

    filename = trim("outputdir/means.h5")

    phibar(:) = 0.0;    phirms(:) = 0.0;

    do i=xstartr(3),xendr(3)
        do j=xstartr(2),xendr(2)
            do k=1,nxmr
                phibar(k) = phibar(k) + phi(k,j,i)
                phirms(k) = phirms(k) + phi(k,j,i)**2
            end do
        end do
    end do

    call MpiSumReal1D(phibar,nxmr)
    call MpiSumReal1D(phirms,nxmr)

    do k=1,nxmr
        phibar(k) = phibar(k)*inyzm
        phirms(k) = sqrt(phirms(k)*inyzm)
    end do

    ! Store index as character string
    write(nstat,"(i5.5)")nint(time/tout)

    if (ismaster) then
        dsetname = trim("phibar/"//nstat)
        call HdfSerialWriteReal1D(dsetname, filename, phibar, nxmr)
        dsetname = trim("phirms/"//nstat)
        call HdfSerialWriteReal1D(dsetname, filename, phirms, nxmr)
    end if

    call MpiBarrier

end subroutine CalcPhiStats

!> Create the groups in the means.h5 file to store the
!! statistics related to the phase-field
subroutine CreatePhaseH5Groups(filename)
    use HDF5

    character(30), intent(in) :: filename
    integer(HID_T) :: file_id, group_id
    integer :: hdf_error

    call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, hdf_error)

    call h5gcreate_f(file_id, "phibar", group_id, hdf_error)
    call h5gclose_f(group_id, hdf_error)
    call h5gcreate_f(file_id, "phirms", group_id, hdf_error)
    call h5gclose_f(group_id, hdf_error)

    call h5fclose_f(file_id, hdf_error)

end subroutine CreatePhaseH5Groups

!> Implicit solve for vx when using direct forcing constraint
subroutine SolveImpEqnUpdate_X_pf
    use local_arrays, only: vx, rhs
    integer :: ic, jc, kc
    real :: fkl(1:nx)       ! RHS vector for implicit solve
    real :: philoc          ! local value of phi on vx-grid
    real :: amkl(1:nxm)     ! Lower diagonal of LHS matrix
    real :: ackl(1:nx)      ! Diagonal of LHS matrix
    real :: apkl(1:nxm)     ! Upper diagonal of LHS matrix
    real :: appk(1:nxm-1)
    integer :: ipkv(1:nx), info
    real :: ackl_b, betadx

    betadx = beta*al

    do ic=xstart(3),xend(3)
        do jc=xstart(2),xend(2)
            
            amkl(:) = 0.0
            ackl(:) = 1.0
            apkl(:) = 0.0

            fkl(1) = -vx(1,jc,ic)
            
            do kc=2,nxm
                philoc = 0.5*(phic(kc,jc,ic) + phic(kc+1,jc,ic))
                if (philoc > pf_direct_force) then ! Solid phase
                    amkl(kc-1) = 0.d0
                    ackl(kc) = 1.d0
                    apkl(kc) = 0.d0
                    fkl(kc) = -vx(kc,jc,ic)
                else ! Liquid phase
                    ackl_b=1.0d0/(1.0d0-ac3ck(kc)*betadx)
                    amkl(kc-1)=-am3ck(kc)*betadx*ackl_b
                    ackl(kc)=1.d0
                    apkl(kc)=-ap3ck(kc)*betadx*ackl_b
                    fkl(kc) = rhs(kc,jc,ic)*ackl_b
                end if
            end do

            fkl(nx) = -vx(nx,jc,ic)

            call dgttrf(nx, amkl, ackl, apkl, appk, ipkv, info)
            call dgttrs('N',nx,1,amkl,ackl,apkl,appk,ipkv,fkl,nx,info)

            do kc=2,nxm
                vx(kc,jc,ic) = vx(kc,jc,ic) + fkl(kc)
            end do

        end do
    end do

end subroutine SolveImpEqnUpdate_X_pf

!> Implicit solve for vy, vz, setting zero in solid
subroutine SolveImpEqnUpdate_YZ_pf(q, rhs, axis)
    real, intent(inout) :: q(1:nx,xstart(2)-lvlhalo:xend(2)+lvlhalo, &
                            & xstart(3)-lvlhalo:xend(3)+lvlhalo)    !> variable to update
    real, intent(inout) :: rhs(1:nx,xstart(2):xend(2),xstart(3):xend(3))    !> rhs to update with
    character, intent(in) :: axis   !> component axis of velocity (to identify grid stagger)
    integer :: ic, jc, kc
    real :: fkl(nxm)       ! RHS vector for implicit solve
    real :: philoc          ! local value of phi on vx-grid
    real :: amkl(nxm-1)   ! Lower diagonal of LHS matrix
    real :: ackl(nxm)     ! Diagonal of LHS matrix
    real :: apkl(nxm-1)   ! Upper diagonal of LHS matrix
    integer :: ipkv(1:nxm), info
    real :: ackl_b, betadx
    real :: appk(nxm-2)

    betadx = beta*al

    do ic=xstart(3),xend(3)
        do jc=xstart(2),xend(2)
            
            amkl(:) = 0.0
            ackl(:) = 1.0
            apkl(:) = 0.0

            do kc=1,nxm
                if (axis=="y") then
                    philoc = 0.5*(phic(kc,jc,ic) + phic(kc,jc+1,ic))
                elseif (axis=="z") then
                    philoc = 0.5*(phic(kc,jc,ic) + phic(kc,jc,ic+1))
                else
                    philoc = phic(kc,jc,ic)
                end if
                if (philoc > pf_direct_force) then ! Solid phase
                    if (kc > 1) amkl(kc-1) = 0.d0
                    ackl(kc) = 1.d0
                    if (kc < nxm) apkl(kc) = 0.d0
                    fkl(kc) = -q(kc,jc,ic)
                else ! Liquid phase
                    ackl_b=1.0d0/(1.0d0-ac3sk(kc)*betadx)
                    if (kc > 1) amkl(kc-1)=-am3sk(kc)*betadx*ackl_b
                    ackl(kc)=1.d0
                    if (kc < nxm) apkl(kc)=-ap3sk(kc)*betadx*ackl_b
                    fkl(kc) = rhs(kc,jc,ic)*ackl_b
                end if
            end do

            call dgttrf(nxm, amkl, ackl, apkl, appk, ipkv, info)
            call dgttrs('N',nxm,1,amkl,ackl,apkl,appk,ipkv,fkl,nxm,info)

            do kc=1,nxm
                q(kc,jc,ic) = q(kc,jc,ic) + fkl(kc)
            end do
        end do
    end do

end subroutine SolveImpEqnUpdate_YZ_pf

end module afid_phasefield