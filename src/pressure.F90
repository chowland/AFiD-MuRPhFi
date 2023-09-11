!> Module containing all subroutines and arrays needed to compute the
!! pressure correction needed to ensure a divergence-free velocity field
module afid_pressure
    use param
    use fftw_params
    use local_arrays, only: pr, dph, dphhalo
    use decomp_2d
    use decomp_2d_fft
    use AuxiliaryRoutines
    use mpih
    implicit none

    logical :: sidewall     !! Flag to determine whether to impose sidewalls in y (using a DCT)

    real, allocatable, dimension(:) :: ak1      !! Modified wavenumber in z
    real, allocatable, dimension(:) :: ak2      !! Modified wavenumber in y

    real, allocatable, dimension(:) :: amphk    !! Lower diagonal Poisson coefficient for pressure solve
    real, allocatable, dimension(:) :: acphk    !! Centre diagonal Poisson coefficient for pressure solve
    real, allocatable, dimension(:) :: apphk    !! Upper diagonal Poisson coefficient for pressure solve

contains

!> Allocate the arrays used for the solution of the pressure field
subroutine InitPressureVars

    !! Modified wavenumbers
    call AllocateReal1DArray(ak1,1,nz)
    call AllocateReal1DArray(ak2,1,ny)
    
    !! Poisson solver coefficients
    call AllocateReal1DArray(amphk,1,nx)
    call AllocateReal1DArray(acphk,1,nx)
    call AllocateReal1DArray(apphk,1,nx)
    
end subroutine InitPressureVars

!> Free up memory from arrays used in pressure correction
subroutine DeallocatePressureVars

    !! Modified wavenumbers
    call DestroyReal1DArray(ak1)
    call DestroyReal1DArray(ak2)

    !! Poisson solver coefficients
    call DestroyReal1DArray(amphk)
    call DestroyReal1DArray(acphk)
    call DestroyReal1DArray(apphk)

    !! Temporary mid-FFT arrays
    if(allocated(dphc)) deallocate(dphc)
    if(allocated(rz1)) deallocate(rz1)
    if(allocated(cz1)) deallocate(cz1)
    if(allocated(ry1)) deallocate(ry1)
    if(allocated(cy1)) deallocate(cy1)
    
end subroutine DeallocatePressureVars

!> Compute the metric terms and modified wavenumbers used in the
!! solution of the pressure correction
subroutine InitPressureSolver
    integer :: nymh, nymp, i, j, nzmh, nzmp
    integer :: kc, km, kp
    real :: ugmmm, a33icc, a33icp
    real :: ao(1:nz), ap(1:ny)

    !! Index counters for wavenumbers
    nymh=nym/2+1
    nymp=nymh+1
    nzmh=nzm/2+1
    nzmp=nzmh+1

    ! call AllocateReal1DArray(ao,1,nz)
    ! call AllocateReal1DArray(ap,1,ny)
    
    !! Construct modified wavenumbers
    !! z
    do i=1,nzmh
        ao(i)=(i-1)*2.d0*pi
    end do
    do i=nzmp,nzm
        ao(i)=-(nzm-i+1)*2.d0*pi
    end do
    do i=1,nzm
        ak1(i)=2.d0*(1.d0-dcos(ao(i)/nzm))*(float(nzm)/zlen)**2
    end do
    
    !! y
    do j=1,nymh
        ap(j)=(j-1)*2.d0*pi
    end do
    do j=nymp,nym
        ap(j)=-(nym-j+1)*2.d0*pi
    end do
    do j=1,nym
        ak2(j)=2.d0*(1.d0-dcos(ap(j)/nym))*(float(nym)/ylen)**2
    end do

    ! call DestroyReal1DArray(ao)
    ! call DestroyReal1DArray(ap)

    !! Initialise tridiagonal matrices for Poisson solver
    do kc=1,nxm
        km = kmv(kc)
        kp = kpv(kc)
        
        a33icc = kmc(kc)*dxq/d3xm(km)
        a33icp = kpc(kc)*dxq/d3xm(kc)
        ugmmm = 1.0d0/d3xc(kc)
        
        amphk(kc) = a33icc*ugmmm
        apphk(kc) = a33icp*ugmmm
        acphk(kc) = -(amphk(kc)+apphk(kc))
    end do

    !! Initialise pencil transposes for the pressure solver
    call decomp_2d_fft_init

    !! Allocate temporary mid-FFT arrays
    allocate(ry1(ph%yst(1):ph%yen(1), &
                 ph%yst(2):ph%yen(2), &
                 ph%yst(3):ph%yen(3)))
    allocate(rz1(ph%zst(1):ph%zen(1), &
                 ph%zst(2):ph%zen(2), &
                 ph%zst(3):ph%zen(3)))
    allocate(cy1(sp%yst(1):sp%yen(1), &
                 sp%yst(2):sp%yen(2), &
                 sp%yst(3):sp%yen(3)))
    allocate(cz1(sp%zst(1):sp%zen(1), &
                 sp%zst(2):sp%zen(2), &
                 sp%zst(3):sp%zen(3)))
    allocate(dphc(sp%xst(1):sp%xen(1), &
                  sp%xst(2):sp%xen(2), &
                  sp%xst(3):sp%xen(3)))

end subroutine InitPressureSolver

!> Compute the local divergence of the intermediate velocity field at
!! every point for the pressure correction step
!! div(u)/al/dt is stored in the array dph after this step
subroutine CalcLocalDivergence
    use local_arrays, only: vx, vy, vz

    integer :: ic, ip, jc, jp, kc, kp
    real :: usdtal, dqcap   

    usdtal = 1.d0/(dt*al)

    do ic=xstart(3),xend(3)
        ip = ic+1
        do jc=xstart(2),xend(2)
            jp = jc+1
            do kc=1,nxm
                kp = kc+1
                dqcap = (vz(kc,jc,ip)-vz(kc,jc,ic))*dz &
                       +(vy(kc,jp,ic)-vy(kc,jc,ic))*dy &
                       +(vx(kp,jc,ic)-vx(kc,jc,ic))*udx3m(kc)
                dph(kc,jc,ic)=dqcap*usdtal
            end do
        end do
    end do
end subroutine CalcLocalDivergence

!> Compute the pressure correction by solving a Poisson equation
subroutine SolvePressureCorrection
    integer :: nymh

    nymh=nym/2+1
    
    !! Compute a 2D Fourier transform in y and z
    call transpose_x_to_y(dph,ry1,ph)
    ! If calling for the first time, plan the Fourier transform
    if (.not. planned) call PlanFourierTransform
    call dfftw_execute_dft_r2c(fwd_guruplan_y,ry1,cy1)
    call transpose_y_to_z(cy1,cz1,sp)
    call dfftw_execute_dft(fwd_guruplan_z,cz1,cz1)
    
    ! Normalise the transformed array. FFTW does not do this automatically
    cz1 = cz1 / (nzm*nym)
    
    ! Transpose back to an x-pencil before the tridiagonal solve
    call transpose_z_to_x(cz1,dphc,sp)
    
    ! Solve the tridiagonal system with complex coefficients
    call SolveTridiagonalPressure

    call transpose_x_to_z(dphc,cz1,sp)
    
    call dfftw_execute_dft(bwd_guruplan_z,cz1,cz1)
    
    call transpose_z_to_y(cz1,cy1,sp)
    
    call dfftw_execute_dft_c2r(bwd_guruplan_y,cy1,ry1)
    
    call transpose_y_to_x(ry1,dph,ph)

end subroutine SolvePressureCorrection

!> Solve the tridiagonal system for the x-derivatives of the Poisson equation
!! for the Fourier transformed pressure
subroutine SolveTridiagonalPressure
    integer :: i, j, k, info
    complex :: acphT_b
    complex :: appph(nxm-2)
    complex :: acphT(nxm), apph(nxm-1), amph(nxm-1)
    integer :: phpiv(nxm)

    ! Construct tridiagonal system for each Fourier mode
    do i=sp%xst(3),sp%xen(3)
        do j=sp%xst(2),sp%xen(2)
            do k=1,nxm
                ! Normalise RHS & LHS to avoid floating point errors
                acphT_b = 1.0/(acphk(k) - ak2(j) - ak1(i))
                dphc(k,j,i) = dphc(k,j,i)*acphT_b
                if (k < nxm) apph(k  ) = apphk(k)*acphT_b
                if (k > 1)   amph(k-1) = amphk(k)*acphT_b
                ! Small perturbation needed to prevent singular matrix
                ! when using uniform grid:
                acphT(k) = 1.0d0 + 1.0d-15 
            end do
            
            ! Factor the tridiagonal matrix
            call zgttrf(nxm, amph, acphT, apph, appph, phpiv, info)
            
            if (info.gt.0) then
                print*,'Singular value found in LAPACK routine zgttrf: info=',info
                print*,'Please try to adjust either NX or STR3 in bou.in'
                call MPI_Abort(MPI_COMM_WORLD,1,ierr)
            endif
            
            ! Solve the tridiagonal system
            call zgttrs('N',nxm,1,amph,acphT,apph,appph,phpiv,dphc(1:nxm,j,i),nxm,info)
            
        end do
    end do
end subroutine SolveTridiagonalPressure

!> Remove the divergent component of velocity from the velocity field
!! using the solved pressure correction
subroutine RemoveDivergence
    use local_arrays, only: vx, vy, vz
    real :: usukm, udy, udz, locdph
    integer :: ic, im, jc, jm, kc, km

    udy = al*dt*dy
    udz = al*dt*dz

    do ic=xstart(3),xend(3)
        im = ic-1
        do jc=xstart(2),xend(2)
            jm = jc-1
            do kc=1,nxm
                km = kmv(kc)
                
                usukm = al*dt*udx3c(kc)
                locdph = dphhalo(kc,jc,ic)
                
                vx(kc,jc,ic) = vx(kc,jc,ic) &
                            - (locdph - dphhalo(km,jc,ic))*usukm
                vy(kc,jc,ic) = vy(kc,jc,ic) &
                            - (locdph - dphhalo(kc,jm,ic))*udy
                vz(kc,jc,ic) = vz(kc,jc,ic) &
                            - (locdph - dphhalo(kc,jc,im))*udz
            end do
       end do
    end do
end subroutine RemoveDivergence

!> Update the total pressure using the solved pressure correction
subroutine CorrectPressure
    integer :: ic, jc, kc
    integer :: im, jm, km
    integer :: ip, jp, kp
    real :: be, amm, acc, app

    be = al*beta
    do ic=xstart(3),xend(3)
        im = ic - 1
        ip = ic + 1
        do jc=xstart(2),xend(2)
            jm = jc - 1
            jp = jc + 1
            do kc=1,nxm
                kp = kpv(kc)
                km = kmv(kc)

                amm = amphk(kc)
                acc = acphk(kc)
                app = apphk(kc)

                pr(kc,jc,ic) = pr(kc,jc,ic) + dphhalo(kc,jc,ic) - be*( &
                    ! dzz(p')
                    (dphhalo(kc,jc,ip) - 2.0*dphhalo(kc,jc,ic) + dphhalo(kc,jc,im))*dzq + &
                    ! dyy(p')
                    (dphhalo(kc,jp,ic) - 2.0*dphhalo(kc,jc,ic) + dphhalo(kc,jm,ic))*dyq + &
                    ! dxx(p')
                    (dphhalo(kp,jc,ic)*app + dphhalo(kc,jc,ic)*acc + dphhalo(km,jc,ic)*amm))
            end do
        end do
    end do
end subroutine CorrectPressure

!> Create FFTW plan
subroutine PlanFourierTransform

    type(fftw_iodim),dimension(1) :: iodim
    type(fftw_iodim),dimension(2) :: iodim_howmany
    integer :: info
    
    iodim(1) % n = nzm
    iodim(1) % is = (sp%zen(1)-sp%zst(1)+1)*(sp%zen(2)-sp%zst(2)+1)
    iodim(1) % os = (sp%zen(1)-sp%zst(1)+1)*(sp%zen(2)-sp%zst(2)+1)

    iodim_howmany(1) % n = (sp%zen(1)-sp%zst(1)+1)
    iodim_howmany(1) % is = 1
    iodim_howmany(1) % os = 1
    iodim_howmany(2) % n = (sp%zen(2)-sp%zst(2)+1)
    iodim_howmany(2) % is = (sp%zen(1)-sp%zst(1)+1)
    iodim_howmany(2) % os = (sp%zen(1)-sp%zst(1)+1)

    ! Construct forward plan for z transform
    fwd_guruplan_z = fftw_plan_guru_dft(1, iodim, &
            2, iodim_howmany, cz1, cz1, &
            FFTW_FORWARD, FFTW_ESTIMATE)
    
    iodim(1) % n = nzm
    ! Construct backward plan for z transform
    bwd_guruplan_z = fftw_plan_guru_dft(1, iodim, &
            2,iodim_howmany,cz1,cz1, &
            FFTW_BACKWARD,FFTW_ESTIMATE)
    
    if (.not.c_associated(bwd_guruplan_z)) then
        if (ismaster) print*,'Failed to create guru plan. You should'
        if (ismaster) print*,'link with FFTW3 before MKL'
        if (ismaster) print*,'Please check linking order.'
        call MPI_Abort(MPI_COMM_WORLD,1,info)
    end if
    
    iodim(1) % n = nym
    iodim(1) % is = ph%yen(1)-ph%yst(1)+1
    iodim(1) % os = sp%yen(1)-sp%yst(1)+1
    
    iodim_howmany(1) % n = (ph%yen(1)-ph%yst(1)+1)
    iodim_howmany(1) % is = 1
    iodim_howmany(1) % os = 1
    iodim_howmany(2) % n = (ph%yen(3)-ph%yst(3)+1)
    iodim_howmany(2) % is = (ph%yen(1)-ph%yst(1)+1) &
                          * (ph%yen(2)-ph%yst(2)+1)
    iodim_howmany(2) % os = (sp%yen(1)-sp%yst(1)+1) &
                          * (sp%yen(2)-sp%yst(2)+1)

    ! Construct forward plan for y transform
    fwd_guruplan_y = fftw_plan_guru_dft_r2c(1, iodim, &
            2, iodim_howmany, ry1, cy1, &
            FFTW_ESTIMATE)
    
    iodim(1) % n = nym
    iodim(1) % is = sp%yen(1)-sp%yst(1)+1
    iodim(1) % os = ph%yen(1)-ph%yst(1)+1

    iodim_howmany(1) % n=(sp%yen(1)-sp%yst(1)+1)
    iodim_howmany(1) % is=1
    iodim_howmany(1) % os=1
    iodim_howmany(2) % n=(sp%yen(3)-sp%yst(3)+1)
    iodim_howmany(2) % is=(sp%yen(1)-sp%yst(1)+1) &
                        * (sp%yen(2)-sp%yst(2)+1)
    iodim_howmany(2) % os=(ph%yen(1)-ph%yst(1)+1) &
                        * (ph%yen(2)-ph%yst(2)+1)
    
    ! Construct backward plan for y transform
    bwd_guruplan_y = fftw_plan_guru_dft_c2r(1,iodim, &
            2,iodim_howmany,cy1,ry1, &
            FFTW_ESTIMATE)

    ! Save that we have planned the FFT
    planned=.true.

end subroutine PlanFourierTransform

end module afid_pressure