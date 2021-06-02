!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: InterpTempMgrd.F90                             !
!    CONTAINS: subroutine InterpTempMgrd                  !
!                                                         ! 
!    PURPOSE: Interpolates temperature onto refined grid  !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine InterpTempMgrd

    use param
    use local_arrays, only: temp
    use mgrd_arrays, only: tempr,cxrs,cyrs,czrs,irangs,jrangs,krangs,tpdv
    use mpih
    use decomp_2d
    use AuxiliaryRoutines
    implicit none

    integer :: ic,jc,kc,icr,jcr,kcr

    real,dimension(4,4,4) :: qv3 
    real,dimension(4,4) :: qv2
    real,dimension(4) :: qv1

    tempr(:,:,:) = 0.d0

    tpdv(:,:,:) = 0.d0 ! Temporary array with exended range in x for interpolation

    ! Fill temporary array with temperature field and BCs
    do ic=xstart(3)-lvlhalo,xend(3)+lvlhalo
        do jc=xstart(2)-lvlhalo,xend(2)+lvlhalo
            tpdv(0,jc,ic) = tempbp(1,jc,ic)
            tpdv(nx,jc,ic) = temptp(1,jc,ic)
            do kc=1,nxm
                tpdv(kc,jc,ic) = temp(kc,jc,ic)
            end do
        end do
    end do

    ! call update_halo(tpdv,lvlhalo)

    do ic=xstart(3)-1,xend(3)
        do jc=xstart(2)-1,xend(2)
            do kc=0,nxm
                qv3 = tpdv(kc-1:kc+2,jc-1:jc+2,ic-1:ic+2)

                do icr=max(krangs(ic),1),min(krangs(ic+1)-1,nzmr)
                    qv2(:,:) = qv3(:,:,1)*czrs(1,icr)+qv3(:,:,2)*czrs(2,icr) &
                              +qv3(:,:,3)*czrs(3,icr)+qv3(:,:,4)*czrs(4,icr)
                    do jcr=max(jrangs(jc),1),min(jrangs(jc+1)-1,nymr)
                        qv1(:) = qv2(:,1)*cyrs(1,jcr)+qv2(:,2)*cyrs(2,jcr) &
                                +qv2(:,3)*cyrs(3,jcr)+qv2(:,4)*cyrs(4,jcr)
                        do kcr=max(irangs(kc),1),min(irangs(kc+1)-1,nxmr)
                            tempr(kcr,jcr,icr) = sum(qv1(1:4)*cxrs(1:4,kcr))
                        end do
                    end do
                end do
            end do
        end do
    end do

    return
end subroutine InterpTempMgrd