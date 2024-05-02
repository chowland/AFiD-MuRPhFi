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
    implicit none
    integer :: ic,jc,n
    
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
    
    call update_halo(temptp,lvlhalo)
    call update_halo(tempbp,lvlhalo)

    !! CJH: May want updating to think about how to properly deal with corners
    if (sidewall) then
        if (xstart(2)==1) then
            do n=1,lvlhalo
                tempbp(1,1-n,:) = tempbp(1,n,:)
                temptp(1,1-n,:) = temptp(1,n,:)
            end do
        end if
        if (xend(2)==nym) then
            do n=1,lvlhalo
                tempbp(1,nym+n,:) = tempbp(1,ny-n,:)
                temptp(1,nym+n,:) = temptp(1,ny-n,:)
            end do
        end if
        if (xstart(3)==1) then
            do n=1,lvlhalo
                tempbp(1,:,1-n) = tempbp(1,:,n)
                temptp(1,:,1-n) = temptp(1,:,n)
            end do
        end if
        if (xend(3)==nzm) then
            do n=1,lvlhalo
                tempbp(1,:,nzm+n) = tempbp(1,:,nz-n)
                temptp(1,:,nzm+n) = temptp(1,:,nz-n)
            end do
        end if
    end if
    
    return
end subroutine SetTempBCs