!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: SetTempBCs.F90                                 !
!    CONTAINS: subroutine SetTempBCs                      !
!                                                         ! 
!    PURPOSE: Initialization routine. Calcuates the       !
!     temperature boundary conditions at the plates       !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine SetTempBCs
    use param
    use decomp_2d
    use afid_moisture, only: beta_q
    use afid_salinity, only: RayS     !2DHorizontalConvection 
    use GridModule


    implicit none
    integer :: ic,jc, ii
    real, dimension(nym):: m_BC_Hot,m_BC_Cold
    real :: pontzero_hot,pontzero_cold



    pontzero_hot =  0.01 * FixValueBCRegion_Length * YLEN
    pontzero_cold = YLEN - 0.01 * FixValueBCRegion_Length * YLEN

    
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
        !ProfileTemp = -1.0 + real(ii - 1) * 2.0/real(nym - 1) 
        !values(ii) = (1.0 + tanh(5*ProfileTemp)) / 2.0  ! JFM Limiting regimes of turbulent horizontal convection. Part I: Intermediate and low Prandtl numbers
        
        !values(ii) = yc(ii)/YLEN !Linear
    end do
    if (FixValueBCRegion_Length == 0) then
        do ic = xstart(3), xend(3)
            do jc = xstart(2), xend(2)
                temptp(1, jc, ic) =0.0
                tempbp(1, jc, ic) = yc(jc)/YLEN !Linear
            end do
        end do
    end if 
    if (str_BC /= 0) then
        call Smooth_non_uniform_BC(m_BC_Hot, nym,pontzero_hot)
        call Smooth_non_uniform_BC(m_BC_Cold, nym,pontzero_cold)

        do jc = xstart(2), xend(2)
            if (ym(jc).ge.pontzero_hot)then
                m_BC_Hot(jc) = 0
            endif 
            if (ym(jc).le.pontzero_cold)then
                m_BC_Cold(jc) = 0
            end if 
        end do 
    end if 
    if (FixValueBCRegion_Length /= 0) then
        do ic = xstart(3), xend(3)
            do jc = xstart(2), xend(2)
                if (str_BC == 0) then
                    if (jc < nym / 2) then
                        if (FixValueBCRegion_Nord_or_Sud == 0) then
                            tempbp(1, jc, ic) = 0.5d0
                            temptp(1, jc, ic) = 0.0d0
                        else
                            temptp(1, jc, ic) = 0.5d0
                            tempbp(1, jc, ic) = 0.0d0
                        end if
                    else
                        if (FixValueBCRegion_Nord_or_Sud == 0) then
                            tempbp(1, jc, ic) = -0.5d0
                            temptp(1, jc, ic) = 0.0d0
                        else
                            temptp(1, jc, ic) = -0.5d0
                            tempbp(1, jc, ic) = 0.0d0
                        end if
                    end if
                else
                    if (jc < nym / 2) then
                        if (FixValueBCRegion_Nord_or_Sud == 0) then
                            tempbp(1, jc, ic) = m_BC_Hot(jc)
                            temptp(1, jc, ic) = 0.0d0
                        else
                            temptp(1, jc, ic) = m_BC_Hot(jc)
                            tempbp(1, jc, ic) = 0.0d0
                        end if
                    else
                        if (FixValueBCRegion_Nord_or_Sud == 0) then
                            tempbp(1, jc, ic) = m_BC_Cold(jc)
                            temptp(1, jc, ic) = 0.0d0
                        else
                            temptp(1, jc, ic) = m_BC_Cold(jc)
                            tempbp(1, jc, ic) = 0.0d0
                        end if
                    end if
                end if
            end do
        end do
    end if
    
    call update_halo(temptp,lvlhalo)
    call update_halo(tempbp,lvlhalo)
    
    return
end subroutine SetTempBCs