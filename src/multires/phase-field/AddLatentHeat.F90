!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: AddLatentHeat.F90                              !
!    CONTAINS: subroutine AddLatentHeat                   !
!                                                         ! 
!    PURPOSE: Take the RHS of phase-field equation,       !
!             interpolate it onto the coarse grid, and    !
!             add it to the temperature equation          !
!             for the latent heat term                    !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine AddLatentHeat

    use param
    use local_arrays, only: hro
    use mgrd_arrays
    use mpih
    use decomp_2d
    use AuxiliaryRoutines
    implicit none

    integer :: ic,jc,kc,icr,jcr,kcr

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

    return

end subroutine AddLatentHeat