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
    use mgrd_arrays, only: rhsr, cxsalc, cysalc, czsalc, irangs, jrangs, krangs, tpdvr
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
            !CJH Note if phi fixed value on boundary, then
            ! dphi/dt=0 on boundary
            tpdvr(0,jc,ic) = 0.d0 - rhsr(1,jc,ic)
            tpdvr(-1,jc,ic) = tpdvr(0,jc,ic)
            do kc=1,nxmr
                tpdvr(kc,jc,ic) = rhsr(kc,jc,ic)
            end do
            tpdvr(nxr,jc,ic) = 0.d0 - rhsr(nxmr,jc,ic)
            tpdvr(nxr+1,jc,ic) = tpdvr(nxr,jc,ic)
        end do
    end do

    ! Fill in halo values
    call update_halo(tpdvr,lvlhalo)

    do ic=xstart(3),xend(3)
        icr = krangs(ic) - 1
        do jc=xstart(2),xend(2)
            jcr = jrangs(jc) - 1
            do kc=1,nxm
                kcr = irangs(kc) - 1

                qv3 = tpdvr(kcr-1:kcr+2,jcr-1:jcr+2,icr-1:icr+2)

                qv2(:,:) = qv3(:,:,1)*czsalc(1,ic)+qv3(:,:,2)*czsalc(2,ic)&
                            +qv3(:,:,3)*czsalc(3,ic)+qv3(:,:,4)*czsalc(4,ic)
                qv1(:)   = qv2(:,1)*cysalc(1,jc)+qv2(:,2)*cysalc(2,jc) &
                            +qv2(:,3)*cysalc(3,jc)+qv2(:,4)*cysalc(4,jc)
                phi_rhs = sum(qv1(1:4)*cxsalc(1:4,kc))

                hro(kc,jc,ic) = hro(kc,jc,ic) + pf_S*phi_rhs*aldt
            end do
        end do
    end do

    return

end subroutine AddLatentHeat