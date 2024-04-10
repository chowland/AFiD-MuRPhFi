!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: SetTempBCs.F90                                 !
!    CONTAINS: subroutine SetTempBCs                      !
!                                                         ! 
!    PURPOSE: Initialization routine. Calcuates the       !
!     temperature boundary conditions at the plates       !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


function smoothstep(edge0, edge1, x)
    real :: edge0, edge1, x
    real :: smoothstep

    x = clamp((x - edge0) / (edge1 - edge0), 0.0, 1.0)

    smoothstep = x * x * (3.0 - 2.0 * x)
end function smoothstep

function clamp(x, lowerlimit, upperlimit)
    real :: x, lowerlimit, upperlimit
    real :: clamp

    clamp = min(max(x, lowerlimit), upperlimit)
end function clamp

subroutine SetTempBCs
    use param
    use decomp_2d
    use afid_moisture, only: beta_q
    use afid_salinity, only: RayS     !2DHorizontalConvection 

    implicit none
    integer :: ic,jc, ii
    real, dimension(nym)::  values
    real :: ProfileTemp,smoothstep

    
    if (rayt>=0) then ! unstable T gradient
        if (inslwN==0) then !Single heated wall case
            do ic=xstart(3),xend(3)
                do jc=xstart(2),xend(2)
                    temptp(1,jc,ic)=0.0
                    tempbp(1,jc,ic)=1.0
                end do
            end do
        else
            do ic=xstart(3),xend(3)
                do jc=xstart(2),xend(2)
                    temptp(1,jc,ic)=-0.5d0
                    tempbp(1,jc,ic)=0.5d0
                end do
            end do
        end if
    else              ! stable T gradient
        do ic=xstart(3),xend(3)
            do jc=xstart(2),xend(2)
                temptp(1,jc,ic)=0.5d0
                tempbp(1,jc,ic)=-0.5d0
            end do
        end do
    end if
    
    if (phasefield) then
        if (rayt>=0) then
            do ic=xstart(3),xend(3)
                do jc=xstart(2),xend(2)
                    temptp(1,jc,ic) = 0.d0
                    tempbp(1,jc,ic) = 1.d0
                end do
            end do
        else
            do ic=xstart(3),xend(3)
                do jc=xstart(2),xend(2)
                    temptp(1,jc,ic) = 1.d0
                    tempbp(1,jc,ic) = 0.d0
                end do
            end do
        end if
    end if

    if (moist) then
        do ic=xstart(3),xend(3)
            do jc=xstart(2),xend(2)
                tempbp(1,jc,ic) = 0.0
                temptp(1,jc,ic) = beta_q - 1.0
            end do
        end do
    end if

   
    !2DHorizontalConvection 
    do ii = 1, nym
        ProfileTemp = -1.0 + real(ii - 1) * 2.0/real(nym - 1) 
        values(ii) = (1.0 + tanh(5*ProfileTemp)) / 2.0  ! JFM Limiting regimes of turbulent horizontal convection. Part I: Intermediate and low Prandtl numbers
        
        values(ii) = xm(ii)/YLEN !Linear
    end do
    if (RayS<0) then 
        do ic=xstart(3),xend(3)
            do jc=xstart(2),xend(2)
                tempbp(1,jc,ic)= values(jc)
               ! if (jc<nym/2) then
                !temptp(1,jc,ic)=0.0d0
                !else 
                !temptp(1,jc,ic)=abs(1)
                !end if 
                temptp(1,jc,ic)=0.0
            end do
        end do
    end if


    if ((active_S==1) .and. (active_T==1) .and. (gAxis==1)) then
        do ic=xstart(3),xend(3)
            do jc=xstart(2),xend(2)
                tempbp(1,jc,ic)= 0
                temptp(1,jc,ic) = smoothstep(YLEN/8, YLEN -YLEN/8, yc(jc))
            end do
        end do
    end if 
    call update_halo(temptp,lvlhalo)
    call update_halo(tempbp,lvlhalo)
    
    return
end subroutine SetTempBCs
