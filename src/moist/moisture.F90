module moisture
    use param
    use decomp_2d, only: xstart, xend, nrank, update_halo
    use AuxiliaryRoutines
    use local_arrays, only: temp, vx, vy, vz, rhs, hro
    implicit none

    real, allocatable, dimension(:,:,:) :: humid    !! Specific humidity field
    real, allocatable, dimension(:,:,:) :: qsat     !! Saturation vapour field
    real, allocatable, dimension(:,:,:) :: rkhumid  !! RK storage array for humidity
    real, allocatable, dimension(:,:,:) :: hkhumid  !! RK storage array for humidity (previous substep)

    real :: gamma_q     !! Coefficient for latent heat effect of condensing vapour
    real :: tau_q       !! Time-scale for condensation
    real :: kappa_q     !! Diffusivity of humidity field
    real :: alpha_q     !! Saturation coefficient
    real :: beta_q      !! Lapse rate coefficient
    real :: pecq        !! Peclet number for humidity (inverse dimensionless diffusivity)

    integer :: qfixS    !! Flag for whether humidity takes fixed value at lower plate
    integer :: qfixN    !! Flag for whether humidity takes fixed value at upper plate

    real, allocatable, dimension(:,:,:) :: humbp    !! Humidity boundary value (bottom plate)
    real, allocatable, dimension(:,:,:) :: humtp    !! Humidity boundary value (top plate)

contains

!> Subroutine to allocate memory for all moisture-related variables
subroutine InitMoistVariables
    ! Allocate array for 3D specific humidity field
    call AllocateReal3DArray(humid,1,nxm,xstart(2)-lvlhalo,xend(2)+lvlhalo,xstart(3)-lvlhalo,xend(3)+lvlhalo)
    ! Allocate array for 3D saturation values
    call AllocateReal3DArray(qsat,1,nxm,xstart(2),xend(2),xstart(3),xend(3))

    ! Runge-Kutta storage arrays
    call AllocateReal3DArray(rkhumid,1,nxm,xstart(2),xend(2),xstart(3),xend(3))
    call AllocateReal3DArray(hkhumid,1,nxm,xstart(2),xend(2),xstart(3),xend(3))

    call AllocateReal3DArray(humbp,1,1,xstart(2),xend(2),xstart(3),xend(3))
    call AllocateReal3DArray(humtp,1,1,xstart(2),xend(2),xstart(3),xend(3))

    ! Parameter values to match Vallis et al (2019) examples

    beta_q = 1.2
    alpha_q = 3.0
    pecq = ren          ! Set Peclet equal to Reynolds
    kappa_q = 1.0/pecq
    gamma_q = 0.19
    tau_q = 5e-5*sqrt(rayt*prat)    ! Vallis et al use diffusive scaling for nondimensionalisation

    qfixN = 1
    qfixS = 1

end subroutine InitMoistVariables
    
!> Deallocate the variables used for Rainy-Benard model
subroutine DeallocateMoistVariables

    call DestroyReal3DArray(humid)
    call DestroyReal3DArray(qsat)
    call DestroyReal3DArray(rkhumid)
    call DestroyReal3DArray(hkhumid)

end subroutine DeallocateMoistVariables

subroutine SetHumidityBCs
    integer :: ic, jc
    
    ! NOTE: tempbp = 0  => humbp = 1
    do ic=xstart(3),xend(3)
        do jc=xstart(2),xend(2)
            humbp(1,jc,ic) = exp(alpha_q*tempbp(1,jc,ic))
            humtp(1,jc,ic) = exp(alpha_q*(temptp(1,jc,ic) - beta_q))
        end do
    end do
end subroutine SetHumidityBCs

subroutine CreateInitialHumidity
    integer :: ic, jc, kc
    real :: rnum

    call random_seed()

    do ic=xstart(3),xend(3)
        do jc=xstart(2),xend(2)
            do kc=1,nxm
                call random_number(rnum)
                humid(kc,jc,ic) = humbp(1,jc,ic) + (humtp(1,jc,ic) - humbp(1,jc,ic))*xm(kc)
                ! humid(kc,jc,ic) = humid(kc,jc,ic) + 1e-3*rnum
                ! humid(kc,jc,ic) = temp(kc,jc,ic)
            end do
        end do
    end do

    call update_halo(humid,lvlhalo)

    call UpdateSaturation

end subroutine CreateInitialHumidity

!>
!! Compute the explicit terms for the specific humidity equation
!! and store the result in rkhumid
subroutine ExplicitHumidity
    integer :: ip, jp
    integer :: ic, jc, kc
    integer :: im, jm
    real :: udy, udz, itau
    real :: dyyq, dzzq, condensation, hqx, hqy, hqz

    udy = 0.5*dy
    udz = 0.5*dz

    itau = 1.0/tau_q

    do ic=xstart(3),xend(3)
        im = ic-1
        ip = ic+1
        do jc=xstart(2),xend(2)
            jm = jc-1
            jp = jc+1
            do kc=1,nxm

                ! x-advection d/dx(vx * q)
                if (kc==1) then
                    hqx = ( &
                          vx(kc+1,jc,ic)*(humid(kc+1,jc,ic) + humid(kc,jc,ic)) &
                        - vx(kc  ,jc,ic)*2.0*humbp(1,jc,ic) &
                    )*0.5*udx3m(kc)
                elseif (kc==nxm) then
                    hqx = ( &
                            vx(kc+1,jc,ic)*2.0*humtp(1,jc,ic) &
                          - vx(kc  ,jc,ic)*(humid(kc,jc,ic) + humid(kc-1,jc,ic)) &
                    )*0.5*udx3m(kc)
                else
                    hqx = ( &
                        vx(kc+1,jc,ic)*(humid(kc+1,jc,ic) + humid(kc  ,jc,ic)) &
                        - vx(kc  ,jc,ic)*(humid(kc  ,jc,ic) + humid(kc-1,jc,ic)) &
                    )*0.5*udx3m(kc)
                end if

                ! y-advection d/dy(vy * q)
                hqy = ( &
                      vy(kc,jp,ic)*(humid(kc,jp,ic) + humid(kc,jc,ic)) &
                    - vy(kc,jc,ic)*(humid(kc,jc,ic) + humid(kc,jm,ic)) &
                )*udy

                ! z-advection d/dz(vz * q)
                hqz = ( &
                      vz(kc,jc,ip)*(humid(kc,jc,ip) + humid(kc,jc,ic)) &
                    - vz(kc,jc,ic)*(humid(kc,jc,ic) + humid(kc,jc,im)) &
                )*udz

                ! yy second derivative of humidity
                dyyq = (humid(kc,jp,ic) - 2.0*humid(kc,jc,ic) + humid(kc,jm,ic))*dyq

                ! zz second derivative of humidity
                dzzq = (humid(kc,jc,ip) - 2.0*humid(kc,jc,ic) + humid(kc,jc,im))*dzq

                rkhumid(kc,jc,ic) = kappa_q*(dyyq + dzzq)

                ! add -(q-q_s)/tau H(q-q_s) to rhs
                if (humid(kc,jc,ic) > qsat(kc,jc,ic)) then
                    condensation = (qsat(kc,jc,ic) - humid(kc,jc,ic))*itau
                    rkhumid(kc,jc,ic) = rkhumid(kc,jc,ic) + condensation
                end if

            end do
        end do
    end do

end subroutine ExplicitHumidity

!> Add the effect of condensation to the buoyancy equation
!! gamma/tau*(q - q_s) H(q-q_s)
subroutine AddCondensation
    integer :: ic,jc,kc
    real :: condensation, gtau

    gtau = gamma_q/tau_q

    do ic=xstart(3),xend(3)
        do jc=xstart(2),xend(2)
            do kc=1,nxm
                if (humid(kc,jc,ic) > qsat(kc,jc,ic)) then
                    condensation = gtau*(humid(kc,jc,ic) - qsat(kc,jc,ic))
                    hro(kc,jc,ic) = hro(kc,jc,ic) + condensation
                end if
            end do
        end do
    end do
    
end subroutine AddCondensation

!> Update the saturation variable qsat using the buoyancy field temp
subroutine UpdateSaturation
    integer :: ic, jc, kc

    ! Recall for this model we use the variable temp for the *buoyancy* field, not the temperature

    do ic=xstart(3),xend(3)
        do jc=xstart(2),xend(2)
            do kc=1,nxm
                qsat(kc,jc,ic) = exp(alpha_q*(temp(kc,jc,ic) - beta_q*xm(kc)))
            end do
        end do
    end do
end subroutine UpdateSaturation

!> Compute the implicit terms for the humidity evolution
subroutine ImplicitHumidity
    integer :: ic, jc, kc
    real :: dxxq, alpec

    alpec = al/pecq

    do ic=xstart(3),xend(3)
        do jc=xstart(2),xend(2)
            do kc=1,nxm
                ! Second xx derivative
                ! Apply lower BC
                if (kc==1) then
                    dxxq = humid(kc+1,jc,ic)*ap3ssk(kc) &
                         + humid(kc  ,jc,ic)*ac3ssk(kc) &
                         - (ap3ssk(kc) + ac3ssk(kc))*humbp(1,jc,ic)*qfixS
                ! Apply upper BC
                elseif (kc==nxm) then
                    dxxq = -(am3ssk(kc) + ac3ssk(kc))*humtp(1,jc,ic)*qfixN &
                         + humid(kc  ,jc,ic)*ac3ssk(kc) &
                         + humid(kc-1,jc,ic)*am3ssk(kc)
                ! Compute dxxq in the interior
                else
                    dxxq = humid(kc+1,jc,ic)*ap3ssk(kc) &
                         + humid(kc  ,jc,ic)*ac3ssk(kc) &
                         + humid(kc-1,jc,ic)*am3ssk(kc)
                end if

                rhs(kc,jc,ic) = (ga*rkhumid(kc,jc,ic) + ro*hkhumid(kc,jc,ic) + alpec*dxxq)*dt
                
                hkhumid(kc,jc,ic) = rkhumid(kc,jc,ic)

            end do
        end do
    end do

    call SolveImpEqnUpdate_Humidity

end subroutine ImplicitHumidity

!> Solve the implicit system for the humidity
!! and update the global variable
subroutine SolveImpEqnUpdate_Humidity
    integer :: ic, jc, kc, info, ipkv(nxm), nrhs
    real :: amkl(nxm), ackl(nxm), apkl(nxm), ackl_b
    real :: amkq(nxm-1), ackq(nxm), apkq(nxm-1), appk(nxm-2), betadx

    betadx = 0.5d0*al*dt/pecq

    ! Construct tridiagonal matrix for LHS
    do kc=1,nxm
        ackl_b = 1.d0/(1.d0 - ac3ssk(kc)*betadx)
        amkl(kc) = -am3ssk(kc)*betadx*ackl_b
        ackl(kc) = 1.d0
        apkl(kc) = -ap3ssk(kc)*betadx*ackl_b
    end do
    
    amkq = amkl(2:nxm)
    ackq = ackl(1:nxm)
    apkq = apkl(1:(nxm-1))

    call dgttrf(nxm, amkq, ackq, apkq, appk, ipkv, info)

    nrhs = (xend(3)-xstart(3)+1)*(xend(2)-xstart(2)+1)
    do ic=xstart(3),xend(3)
        do jc=xstart(2),xend(2)
            do kc=1,nxm
                ackl_b = 1.0/(1.0 - ac3ssk(kc)*betadx)
                rhs(kc,jc,ic) = rhs(kc,jc,ic)*ackl_b
            end do
        end do
    end do

    call dgttrs('N', nxm, nrhs, amkq, ackq, apkq, appk, ipkv, rhs(1:nxm,:,:), nxm, info)

    do ic=xstart(3),xend(3)
        do jc=xstart(2),xend(2)
            do kc=1,nxm
                humid(kc,jc,ic) = humid(kc,jc,ic) + rhs(kc,jc,ic)
            end do
        end do
    end do
    
end subroutine SolveImpEqnUpdate_Humidity

!> Calculate vertical profiles of moisture-related statistics and
!! store them in means.h5
subroutine CalcMoistStats

    real, dimension(nxm) :: qbar    !! Horizontally-averaged specific humidity
    real, dimension(nxm) :: qrms    !! Horizontally-averaged rms humidity
    real, dimension(nxm) :: qsbar   !! Horizontally-averaged saturation humidity
    real, dimension(nxm) :: qrel    !! Horizontally-averaged relative humidity (q/qs)

    real, dimension(nxm) :: vxq     !! Advective flux of moisture (x)
    real, dimension(nxm) :: vyq     !! Advective flux of moisture (y)
    real, dimension(nxm) :: vzq     !! Advective flux of moisture (z)

    real :: inyzm   !! 1.0/nym/nzm

    character(30) :: dsetname   !! Dataset name for HDF5 file
    character(30) :: filename   !! HDF5 file name for statistic storage
    character( 5) :: nstat      !! Character string of statistic index

    integer :: i, j, k, ip, jp

    inyzm = 1.0/nym/nzm

    filename = trim("outputdir/means.h5")

    do i=xstart(3),xend(3)
        ip = i + 1
        do j=xstart(2),xend(2)
            jp = j + 1
            do k=1,nxm
                qbar(k) = qbar(k) + humid(k,j,i)
                qsbar(k) = qsbar(k) + qsat(k,j,i)
                qrms(k) = qrms(k) + humid(k,j,i)**2
                qrel(k) = qrel(k) + humid(k,j,i)/qsat(k,j,i)
                vxq(k) = vxq(k) + 0.5*(vx(k,j,i) + vx(k+1,j,i))*humid(k,j,i)
                vyq(k) = vyq(k) + 0.5*(vy(k,j,i) + vy(k,j,i))*humid(k,j,i)
                vzq(k) = vzq(k) + 0.5*(vz(k,j,i) + vz(k,j,i))*humid(k,j,i)
            end do
        end do
    end do

    call MpiSumReal1D(qbar,nxm)
    call MpiSumReal1D(qsbar,nxm)
    call MpiSumReal1D(qrms,nxm)
    call MpiSumReal1D(qrel,nxm)
    call MpiSumReal1D(vxq,nxm)
    call MpiSumReal1D(vyq,nxm)
    call MpiSumReal1D(vzq,nxm)

    do k=1,nxm
        qbar(k) = qbar(k)*inyzm
        qsbar(k) = qsbar(k)*inyzm
        qrms(k) = sqrt(qrms(k)*inyzm)
        qrel(k) = qrel(k)*inyzm
        vxq(k) = vxq(k)*inyzm
        vyq(k) = vyq(k)*inyzm
        vzq(k) = vzq(k)*inyzm
    end do

    ! Store index as character string
    write(nstat,"(i5.5)")nint(time/tout)

    
    if (ismaster) then
        dsetname = trim("qbar/"//nstat)
        call HdfSerialWriteReal1D(dsetname, filename, qbar, nxm)
        dsetname = trim("qsbar/"//nstat)
        call HdfSerialWriteReal1D(dsetname, filename, qsbar, nxm)
        dsetname = trim("qrms/"//nstat)
        call HdfSerialWriteReal1D(dsetname, filename, qrms, nxm)
        dsetname = trim("qrel/"//nstat)
        call HdfSerialWriteReal1D(dsetname, filename, qrel, nxm)
        dsetname = trim("vxq/"//nstat)
        call HdfSerialWriteReal1D(dsetname, filename, vxq, nxm)
        dsetname = trim("vyq/"//nstat)
        call HdfSerialWriteReal1D(dsetname, filename, vyq, nxm)
        dsetname = trim("vzq/"//nstat)
        call HdfSerialWriteReal1D(dsetname, filename, vzq, nxm)
    end if

    call MpiBarrier

end subroutine CalcMoistStats

!> Create the groups in the means.h5 file to store the
!! moist-related statistics
subroutine CreateMoistH5Groups(filename)
    use HDF5
    
    character(30), intent(in) :: filename
    integer(HID_T) :: file_id, group_id
    integer :: hdf_error

    call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, hdf_error)

    call h5gcreate_f(file_id, "qbar", group_id, hdf_error)
    call h5gclose_f(group_id, hdf_error)
    call h5gcreate_f(file_id, "qsbar", group_id, hdf_error)
    call h5gclose_f(group_id, hdf_error)
    call h5gcreate_f(file_id, "qrms", group_id, hdf_error)
    call h5gclose_f(group_id, hdf_error)
    call h5gcreate_f(file_id, "qrel", group_id, hdf_error)
    call h5gclose_f(group_id, hdf_error)
    call h5gcreate_f(file_id, "vxq", group_id, hdf_error)
    call h5gclose_f(group_id, hdf_error)
    call h5gcreate_f(file_id, "vyq", group_id, hdf_error)
    call h5gclose_f(group_id, hdf_error)
    call h5gcreate_f(file_id, "vzq", group_id, hdf_error)
    call h5gclose_f(group_id, hdf_error)

    call h5fclose_f(file_id, hdf_error)

end subroutine CreateMoistH5Groups

end module moisture