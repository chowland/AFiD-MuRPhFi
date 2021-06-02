!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: InterpPhiMgrd.F90                              !
!    CONTAINS: subroutine InterpPhiMgrd                   !
!                                                         ! 
!    PURPOSE: Interpolates phase-field variable onto      !
!               the coarse base grid                      !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine InterpPhiMgrd

    use param
    use mgrd_arrays, only: phi,phic,cxsalc,cysalc,czsalc,irangs,jrangs,krangs,tpdvr
    use mpih
    use decomp_2d
    use AuxiliaryRoutines
    implicit none

    integer  :: ic,jc,kc, icr,jcr,kcr

    real,dimension(4,4,4) :: qv3 
    real,dimension(4,4) :: qv2
    real,dimension(4) :: qv1

    ! Interpolate to coarse grid here. A better option is to apply
    ! a box filter.
    phic(:,:,:) = 0.d0

    tpdvr(:,:,:) = 0.d0 ! Temporary array with extended range for interpolation
    do ic=xstartr(3)-lvlhalo,xendr(3)+lvlhalo
        do jc=xstartr(2)-lvlhalo,xendr(2)+lvlhalo
            tpdvr(0,jc,ic) = 2.d0*0.d0 & !phibp(1,jc,ic) &
                                - phi(1,jc,ic)
            tpdvr(-1,jc,ic) = tpdvr(0,jc,ic) !CJH This needs checking in the stencil
            do kc=1,nxmr
                tpdvr(kc,jc,ic) = phi(kc,jc,ic)
            end do
            tpdvr(nxr,jc,ic) = 2.d0*1.d0 &!phitp(1,jc,ic) &
                                - phi(nxmr,jc,ic)
            tpdvr(nxr+1,jc,ic) = tpdvr(nxr,jc,ic)
        end do
    end do


    do ic=xstart(3),xend(3)
        icr = krangs(ic)-1
        do jc=xstart(2),xend(2)
            jcr = jrangs(jc)-1
            do kc=1,nxm
                kcr = irangs(kc)-1

                qv3 = tpdvr(kcr-1:kcr+2,jcr-1:jcr+2,icr-1:icr+2)

                qv2(:,:) = qv3(:,:,1)*czsalc(1,ic)+qv3(:,:,2)*czsalc(2,ic)&
                            +qv3(:,:,3)*czsalc(3,ic)+qv3(:,:,4)*czsalc(4,ic)
                qv1(:)   = qv2(:,1)*cysalc(1,jc)+qv2(:,2)*cysalc(2,jc) &
                            +qv2(:,3)*cysalc(3,jc)+qv2(:,4)*cysalc(4,jc)
                phic(kc,jc,ic) = sum(qv1(1:4)*cxsalc(1:4,kc))

            enddo
        enddo
    enddo

    return

end subroutine InterpPhiMgrd