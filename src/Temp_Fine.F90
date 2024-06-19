!> Module adding evolution of Temperature on the rifine grid to AFiD
!!  is evolved on a refined grid, and therefore
!! this module depends on the multiple-resolution and
!! interpolation modules
module afid_Termperature_Fine
    use param
    use local_arrays
    use mgrd_arrays
    use decomp_2d, only: xstart, xend, xstartr, xendr, update_halo
    use AuxiliaryRoutines
    use HermiteInterpolations, only: interpolate_xyz_to_coarse, interpolate_xyz_to_coarse_fast
    use ibm_param, only: solidr
    implicit none

    real, allocatable, dimension(:,:,:) :: temp_fine      !! Temperature field
    real, allocatable, dimension(:,:,:) :: temp_fine_corse     !! Interpolated Temperature field on coarse grid
    real, allocatable, dimension(:,:,:) :: rutemp_fine     !! RK storage array for salinity (previous substep)
    real, allocatable, dimension(:,:,:) :: htemp_fine      !! RK storage array for salinity

    real :: byctemp_fine   !! Buoyancy prefactor for temp_fine 

  
    real, allocatable, dimension(:,:,:) :: temp_fine_bp    !! temp_fine  boundary value (lower plate)
    real, allocatable, dimension(:,:,:) :: temp_fine_tp    !! temp_fine  boundary value (upper plate)
   
    real, allocatable, dimension(:,:) ::temp_ap3sskr      !! Upper diagonal derivative coefficient for salinity
    real, allocatable, dimension(:,:) :: temp_ac3sskr      !! Diagonal derivative coefficient for salinity
    real, allocatable, dimension(:,:) :: temp_am3sskr 
contains


!> Subroutine to allocate memory for temp_fine -related variables
subroutine Init_Termperature_Fine_Variables
    use  local_arrays

    ! Boundary planes
    call AllocateReal3DArray(temp_fine_bp,1,1,xstartr(2)-lvlhalo,xendr(2)+lvlhalo,xstartr(3)-lvlhalo,xendr(3)+lvlhalo)
    call AllocateReal3DArray(temp_fine_tp,1,1,xstartr(2)-lvlhalo,xendr(2)+lvlhalo,xstartr(3)-lvlhalo,xendr(3)+lvlhalo)

    ! Main arrays with ghost cells
    call AllocateReal3DArray(temp_fine,1,nxr,xstartr(2)-lvlhalo,xendr(2)+lvlhalo,xstartr(3)-lvlhalo,xendr(3)+lvlhalo)
  

    ! Coarse array
    call AllocateReal3DArray(temp_fine_corse,1,nx,xstart(2)-lvlhalo,xend(2)+lvlhalo,xstart(3)-lvlhalo,xend(3)+lvlhalo)

    call AllocateReal3DArray(htemp_fine,1,nxr,xstartr(2),xendr(2),xstartr(3),xendr(3))
    call AllocateReal3DArray(rutemp_fine,1,nxr,xstartr(2),xendr(2),xstartr(3),xendr(3))

      ! Second derivative coefficients
    call AllocateReal2DArray(temp_ap3sskr,1,nxr,1,2)
    call AllocateReal2DArray(temp_ac3sskr,1,nxr,1,2)
    call AllocateReal2DArray(temp_am3sskr,1,nxr,1,2)

    if(.not.salinity)then
    call AllocateReal3DArray(vxr,1,nxr,xstartr(2)-lvlhalo,xendr(2)+lvlhalo,xstartr(3)-lvlhalo,xendr(3)+lvlhalo)
    call AllocateReal3DArray(vyr,1,nxr,xstartr(2)-lvlhalo,xendr(2)+lvlhalo,xstartr(3)-lvlhalo,xendr(3)+lvlhalo)
    call AllocateReal3DArray(vzr,1,nxr,xstartr(2)-lvlhalo,xendr(2)+lvlhalo,xstartr(3)-lvlhalo,xendr(3)+lvlhalo)
end if
end subroutine Init_Termperature_Fine_Variables

!> Deallocate the variables used for evolving temp_fine 
subroutine Deallocate_Termperature_Fine_Variables
    use  local_arrays

    ! Boundary planes
    call DestroyReal3DArray(temp_fine_bp)
    call DestroyReal3DArray(temp_fine_tp)

    ! Main array
    call DestroyReal3DArray(temp_fine)


    ! Coarse array
    call DestroyReal3DArray(temp_fine_corse)

    call DestroyReal3DArray(rutemp_fine)
    call DestroyReal3DArray(htemp_fine)

! Second derivative coefficients
    call DestroyReal2DArray(temp_ap3sskr)
    call DestroyReal2DArray(temp_ac3sskr)
    call DestroyReal2DArray(temp_am3sskr)
    if(.not.salinity)then
        call DestroyReal3DArray(vxr)
        call DestroyReal3DArray(vyr)
        call DestroyReal3DArray(vzr)
    end if

end subroutine Deallocate_Termperature_Fine_Variables

subroutine Set_Termperature_Fine_BCs
    use GridModule
    integer :: i, j,ii
    real, dimension(nymr):: m_BC_Hot,m_BC_Cold
    real :: pontzero_hot,pontzero_cold

    pontzero_hot =  0.01 * FixValueBCRegion_Length * YLEN
    pontzero_cold = YLEN - 0.01 * FixValueBCRegion_Length * YLEN
  if(FixValueBCRegion_Length==0)then
    do i=xstartr(3),xendr(3)
        do j=xstartr(2),xendr(2)
            temp_fine_tp(1,j,i) = -0.5
            temp_fine_bp(1,j,i) = 0.5d0
        end do 
    end do
  end if 

if (FixValueBCRegion_Length == 0) then
    do i=xstartr(3),xendr(3)
        do j=xstartr(2),xendr(2)
            temp_fine_tp(1, j, i) =0.0
            temp_fine_bp(1, j, i) = ycr(j)/YLEN !Linear
        end do
    end do
end if 
  if (str_BC /= 0) then
    call Smooth_non_uniform_BC(m_BC_Hot, nymr,pontzero_hot)
    call Smooth_non_uniform_BC(m_BC_Cold, nymr,pontzero_cold)
    do j=xstartr(2),xendr(2)
            if (ymr(j).ge.pontzero_hot)then
                m_BC_Hot(j) = 0
            end if 
            if (ymr(j) .le. pontzero_cold)then
                m_BC_Cold(j) = 0
            end if 
        end do 
    end if 

    if  (FixValueBCRegion_Length/=0) then  
        do i=xstartr(3),xendr(3)
            do j=xstartr(2),xendr(2)
                if (str_BC == 0) then
                if (j<nymr/2) then
                    if ( FixValueBCRegion_Nord_or_Sud==0) then
                        temp_fine_tp(1,j,i) = 0.0
                        temp_fine_bp(1,j,i) = 0.5d0
                    else  
                        temp_fine_tp(1,j,i) = 0.5d0
                        temp_fine_bp(1,j,i) = 0.0
                    end if 
                else 
                    if ( FixValueBCRegion_Nord_or_Sud==0) then
                        temp_fine_tp(1,j,i) = 0.0
                        temp_fine_bp(1,j,i) = - 0.5d0
                   else  
                         temp_fine_tp(1,j,i) =  - 0.5d0
                        temp_fine_bp(1,j,i) = 0.0
                   end if 
                end if 
            else 

                if (j<nymr/2) then
                    if ( FixValueBCRegion_Nord_or_Sud==0) then
                        temp_fine_tp(1,j,i) = 0.0
                        temp_fine_bp(1,j,i) =  m_BC_Hot(j)
                    else  
                        temp_fine_tp(1,j,i) = m_BC_Hot(j)
                        temp_fine_bp(1,j,i) = 0.0
                    end if 
                else 
                    if ( FixValueBCRegion_Nord_or_Sud==0) then
                        temp_fine_tp(1,j,i) = 0.0
                        temp_fine_bp(1,j,i) = m_BC_Cold(j)
                   else  
                        temp_fine_tp(1,j,i) =  m_BC_Cold(j)
                        temp_fine_bp(1,j,i) = 0.0
                   end if 
                end if 


            end if 
            end do
        end do
    end if
   ! Update halo for interpolation routine
   call update_halo(temp_fine_tp,lvlhalo)
   call update_halo(temp_fine_bp,lvlhalo)

   ! Extend to sidewall halos (could be modified for other BCs...)
   if (sidewall) then
       if (xstartr(2)==1) then
           do j=1,lvlhalo
            temp_fine_tp(1,1-j,:) = temp_fine_tp(1,j,:)
            temp_fine_bp(1,1-j,:) = temp_fine_bp(1,j,:)
           end do
       end if
       if (xendr(2)==nymr) then
           do j=1,lvlhalo
            temp_fine_tp(1,nymr+j,:) = temp_fine_tp(1,nymr+1-j,:)
            temp_fine_bp(1,nymr+j,:) = temp_fine_bp(1,nymr+1-j,:)
           end do
       end if

       if (xstartr(3)==1) then
           do i=1,lvlhalo
            temp_fine_tp(1,:,1-i) = temp_fine_tp(1,:,i)
            temp_fine_bp(1,:,1-i) = temp_fine_bp(1,:,i)
           end do
       end if
       if (xendr(3)==nzmr) then
           do i=1,lvlhalo
            temp_fine_tp(1,:,nzmr+i) = temp_fine_tp(1,:,nzmr+1-i)
            temp_fine_bp(1,:,nzmr+i) = temp_fine_bp(1,:,nzmr+1-i)
           end do
       end if
   end if
end subroutine Set_Termperature_Fine_BCs
subroutine CreateInitia_Termperature_Fine

    integer :: i, j, k
    real :: varptb
    

    do i=xstartr(3),xendr(3)
        do j=xstartr(2),xendr(2)
            do k=1,nxmr
                call random_number(varptb)
                varptb = (varptb - 0.5) * 0.01
                temp_fine(k,j,i) = 0 !+varptb
            end do
        end do
    end do

end subroutine CreateInitia_Termperature_Fine

subroutine Explicit_Termperature_Fine
    use  local_arrays
    !use afid_salinity,only: vxr,vyr,vzr
    integer :: ic, jc, kc
    integer :: im, jm, km
    integer :: ip, jp, kp
    
    real :: udyr, udzr, udyrq, udzrq
    real :: aldt
    real, dimension(1:nxmr) :: sdx
    real :: htemp_finex, htemp_finey, htemp_finez
    real :: dyytemp_fine, dzztemp_fine
    ! Advection coefficients
    udyr = 0.5d0*dyr
    udzr = 0.5d0*dzr
    ! Diffusion coefficients
    udyrq = dyqr/pect
    udzrq = dzqr/pect

    ! x-advection coefficients
    do kc=1,nxmr
        sdx(kc) = 0.5*dxr/g3rmr(kc)
    end do
    ! Time advancing pre-factor
    aldt = 1.0/al/dt
    do ic=xstartr(3),xendr(3)
        im = ic - 1
        ip = ic + 1
        do jc=xstartr(2),xendr(2)
            jm = jc - 1
            jp = jc + 1
            do kc=1,nxmr
                km = kc - 1
                kp = kc + 1

                ! x-advection d/dx (vx * S)
                if (kc==1) then
                    htemp_finex = ( &
                          vxr(kp,jc,ic)*(temp_fine(kp,jc,ic) + temp_fine(kc,jc,ic)) &
                        - vxr(kc,jc,ic)*2.d0*temp_fine_bp(1,jc,ic) &
                    )*udx3mr(kc)*0.5d0
                elseif (kc==nxmr) then
                    htemp_finex = ( &
                          vxr(kp,jc,ic)*2.d0*temp_fine_tp(1,jc,ic) &
                        - vxr(kc,jc,ic)*(temp_fine(kc,jc,ic) + temp_fine(km,jc,ic)) &
                    )*udx3mr(kc)*0.5d0
                else
                    htemp_finex = ( &
                          vxr(kp,jc,ic)*(temp_fine(kp,jc,ic) + temp_fine(kc,jc,ic)) &
                        - vxr(kc,jc,ic)*(temp_fine(kc,jc,ic) + temp_fine(km,jc,ic)) &
                    )*udx3mr(kc)*0.5d0
                end if

                ! y-advection d/dy(vy * S)
                htemp_finey = ( &
                      vyr(kc,jp,ic)*(temp_fine(kc,jp,ic) + temp_fine(kc,jc,ic)) &
                    - vyr(kc,jc,ic)*(temp_fine(kc,jc,ic) + temp_fine(kc,jm,ic)) &
                )*udyr

                ! z-advection d/dz(vz * S)
                htemp_finez = ( &
                      vzr(kc,jc,ip)*(temp_fine(kc,jc,ip) + temp_fine(kc,jc,ic)) &
                    - vzr(kc,jc,ic)*(temp_fine(kc,jc,ic) + temp_fine(kc,jc,im)) &
                )*udzr


                    ! yy second derivative of temp_fine
                    dyytemp_fine = (temp_fine(kc,jp,ic) - 2.0*temp_fine(kc,jc,ic) + temp_fine(kc,jm,ic))*udyrq
                    ! zz second derivative of temp_fine
                    dzztemp_fine = (temp_fine(kc,jc,ip) - 2.0*temp_fine(kc,jc,ic) + temp_fine(kc,jc,im))*udzrq

                ! Sum explicit terms
                htemp_fine(kc,jc,ic) = -(htemp_finex + htemp_finey + htemp_finez) + dyytemp_fine + dzztemp_fine
            end do
        end do
    end do


end subroutine Explicit_Termperature_Fine



subroutine Implicit_Termperature_Fine
    use param,only: pect,TfixN,TfixS,Robin,FixValueBCRegion_Nord_or_Sud,FixValueBCRegion_Length,al,ga
    integer :: ic, jc, kc,ii
    real :: dxx_Termperature_Fine, alpet,  FlagBC_Nord, FlagBC_Sud
    alpet = al/pect
    if (FixValueBCRegion_Length==0) then
        FlagBC_Sud = TfixS
        FlagBC_Nord = TfixN
    else if  (FixValueBCRegion_Length/=0 .and. FixValueBCRegion_Nord_or_Sud==0) then
        FlagBC_Nord = TfixN
    else if  (FixValueBCRegion_Length/=0 .and. FixValueBCRegion_Nord_or_Sud==1) then
        FlagBC_Sud = TfixS
    end if
    
    do ic=xstartr(3),xendr(3)
        do jc=xstartr(2),xendr(2)
            do kc=1,nxmr
                if (FixValueBCRegion_Length/=0) then
                    if (FixValueBCRegion_Nord_or_Sud==0) then
                        if (ymr(jc) <= 0.01 * FixValueBCRegion_Length * YLEN .or. &
                            ymr(jc) >= YLEN - 0.01 * FixValueBCRegion_Length * YLEN) then
                            ii = 1
                            FlagBC_Sud = 1
                         else 
                            ii = 2
                            FlagBC_Sud = 0
                         end if
                    else if (FixValueBCRegion_Nord_or_Sud == 1) then
                        if (ymr(jc) <= 0.01 * FixValueBCRegion_Length * YLEN .or. &
                            ymr(jc) >= YLEN - 0.01 * FixValueBCRegion_Length * YLEN) then
                            ii = 1
                            FlagBC_Nord = 1
                         else 
                            ii = 2
                            FlagBC_Nord = 0
                         end if
                    end if
                else 
                ii = 1
                end if
                ! Second xx derivative
                ! Apply lower BC
                if (kc==1) then
                    dxx_Termperature_Fine= temp_fine(kc+1,jc,ic)*temp_ap3sskr(kc,ii) &
                        + temp_fine(kc  ,jc,ic)*temp_ac3sskr(kc,ii) &
                        - (temp_ap3sskr(kc,ii) + temp_ac3sskr(kc,ii))*temp_fine_bp(1,jc,ic)*FlagBC_Sud
           
                ! Apply upper BC
                elseif (kc==nxmr) then
 
                        dxx_Termperature_Fine= temp_fine(kc  ,jc,ic)*temp_ac3sskr(kc,ii) &
                        + temp_fine(kc-1,jc,ic)*temp_am3sskr(kc,ii) &
                        - (temp_am3sskr(kc,ii) + temp_ac3sskr(kc,ii))*temp_fine_tp(1,jc,ic)*FlagBC_Nord

                else
                    dxx_Termperature_Fine= temp_fine(kc+1,jc,ic)*temp_ap3sskr(kc,ii) &
                        + temp_fine(kc  ,jc,ic)*temp_ac3sskr(kc,ii) &
                        + temp_fine(kc-1,jc,ic)*temp_am3sskr(kc,ii)
                end if
          
                rhsr(kc,jc,ic) = (ga*htemp_fine(kc,jc,ic) + ro*rutemp_fine(kc,jc,ic) + alpet*dxx_Termperature_Fine)*dt
                rutemp_fine(kc,jc,ic) = htemp_fine(kc,jc,ic)
     
            end do
        end do
    end do
    call SolveImpEqnUpdate_Termperature_Fine
end subroutine Implicit_Termperature_Fine


subroutine SolveImpEqnUpdate_Termperature_Fine
    use param,only: pect,TfixN,TfixS
    real :: betadx, ackl_b
    integer :: ic, jc, kc, nrhs, ipkv(nxr), info,ii
    real :: amkT(nxmr-1), ackT(nxmr), apkT(nxmr-1), appk(nxmr-2)

    betadx = 0.5d0*al*dt/pect

    ! Construct tridiagonal matrix for LHS
    do jc=xstart(2),xend(2)
        if (FixValueBCRegion_Length==0) then
            ii = 1
        else 
            if (ymr(jc) <= 0.01 * FixValueBCRegion_Length * YLEN .or. &
                 ymr(jc) >= YLEN - 0.01 * FixValueBCRegion_Length * YLEN) then
                
                    ii = 1
            else 
                    ii = 2
            end if        

        end if
    do kc=1,nxmr
  
        ackl_b = 1.0d0/(1. - temp_ac3sskr(kc,ii)*betadx)
        if (kc > 1) amkT(kc-1) = -temp_am3sskr(kc,ii)*betadx*ackl_b
        ackT(kc) = 1.d0
        if (kc < nxmr) apkT(kc) = -temp_ap3sskr(kc,ii)*betadx*ackl_b
    end do
    end do
    ! Factor the tridiagonal matrix
    call dgttrf(nxmr,amkT,ackT,apkT,appk,ipkv,info)

    ! Rescale RHS to match rescaling of LHS
    nrhs=(xendr(3)-xstartr(3)+1)*(xendr(2)-xstartr(2)+1)
    do ic=xstartr(3),xendr(3)
        do jc=xstartr(2),xendr(2)
            if (FixValueBCRegion_Length==0) then
                ii = 1
             
            else 
                if (ymr(jc) <= 0.01 * FixValueBCRegion_Length * YLEN .or. &
                     ymr(jc) >= YLEN - 0.01 * FixValueBCRegion_Length * YLEN) then
                    
                        ii = 1
                else 
                        ii = 2
                end if        
    
            end if

            do kc=1,nxmr

                ackl_b = 1.0/(1.0 - temp_ac3sskr(kc,ii)*betadx)
                rhsr(kc,jc,ic) = rhsr(kc,jc,ic)*ackl_b
            end do
        end do
    end do

    ! Solve tridiagonal system
    call dgttrs('N',nxmr,nrhs,amkT,ackT,apkT,appk,ipkv,rhsr,nxmr,info)
    ! Update global variable
    do ic=xstartr(3),xendr(3)
        do jc=xstartr(2),xendr(2)
           do kc=1,nxmr
            temp_fine(kc,jc,ic) = temp_fine(kc,jc,ic) + rhsr(kc,jc,ic)
            end do
         end do
     end do

end subroutine SolveImpEqnUpdate_Termperature_Fine



subroutine Interpl_temp_fineMultigrid
    use param
    integer :: icr, jcr, kcr
    real ::  FlagBC_Nord, FlagBC_Sud,dxr_Top
    ! Set coarse salinity array to zero
    temp_fine_corse(:,:,:) = 0.d0
    dxr_Top = alx3 - xmr(nxmr)

    ! Extend refined array in wall-normal direction to give sufficient points
    ! for cubic interpolation



    



    if (FixValueBCRegion_Length==0) then
        FlagBC_Sud = TfixS
        FlagBC_Nord = TfixN
    else if  (FixValueBCRegion_Length/=0 .and. FixValueBCRegion_Nord_or_Sud==0) then
        FlagBC_Nord = TfixN
    else if  (FixValueBCRegion_Length/=0 .and. FixValueBCRegion_Nord_or_Sud==1) then
        FlagBC_Sud = TfixS
    end if
    do icr=xstartr(3)-lvlhalo,xendr(3)+lvlhalo
        do jcr=xstartr(2)-lvlhalo,xendr(2)+lvlhalo
            do kcr=1,nxmr
                temp_pdvr(kcr,jcr,icr) = temp_fine(kcr,jcr,icr)
            end do
            



            if (FixValueBCRegion_Length/=0) then
                if (FixValueBCRegion_Nord_or_Sud==0) then
                    if (ymr(jcr) <= 0.01 * FixValueBCRegion_Length * YLEN .or. &
                        ymr(jcr) >= YLEN - 0.01 * FixValueBCRegion_Length * YLEN) then
                        FlagBC_Sud = 1
                     else 
                        FlagBC_Sud = 0
                     end if
                else if (FixValueBCRegion_Nord_or_Sud == 1) then
                    if (ymr(jcr) <= 0.01 * FixValueBCRegion_Length * YLEN .or. &
                        ymr(jcr) >= YLEN - 0.01 * FixValueBCRegion_Length * YLEN) then
                        FlagBC_Nord = 1
                     else 
                        FlagBC_Nord = 0
                     end if
                end if
            end if
            if (FlagBC_Sud==1) then
                temp_pdvr(0,jcr,icr) = 2.0*temp_fine_bp(1,jcr,icr) - temp_fine(1,jcr,icr)
            else
                temp_pdvr(0,jcr,icr) = temp_fine(1,jcr,icr)
             
            end if
            if (FlagBC_Nord==1) then
                temp_pdvr(nxr,jcr,icr) = 2.0*temp_fine_tp(1,jcr,icr) - temp_fine(nxmr,jcr,icr)
            else
                temp_pdvr(nxr,jcr,icr) = temp_fine(nxmr,jcr,icr)
            end if
            
        end do
    end do
    ! Interpolate the refined field to the coarse grid, storing in salc
    !call interpolate_xyz_to_coarse_fast(temp_pdvr, temp_fine_corse(1:nxm,:,:), "sal")

        call interpolate_xyz_to_coarse(temp_pdvr, temp_fine_corse(1:nxm,:,:))
  
 
end subroutine Interpl_temp_fineMultigrid

!> Add buoyancy contribution from the salinity to one of the
!! momentum forcing arrays
subroutine AddBuoyancy_temp_fine(rkv)
    use param, only:byct
    real, dimension(:,xstart(2):,xstart(3):), intent(inout) :: rkv
    integer :: ic, jc, kc

    do ic=xstart(3),xend(3)
        do jc=xstart(2),xend(2)
            do kc=1,nxm
                rkv(kc,jc,ic) = rkv(kc,jc,ic) + byct*temp_fine_corse(kc,jc,ic)
            end do
        end do
    end do
 
end subroutine AddBuoyancy_temp_fine

end module afid_Termperature_Fine