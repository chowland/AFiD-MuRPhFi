subroutine InterpSalMgrd

    use param
    use mgrd_arrays, only: sal,salc,cxsalc,cysalc,czsalc,irangr,jrangr,krangr,tpdvr
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
    salc(:,:,:) = 0.d0

    do icr=xstartr(3)-lvlhalo,xendr(3)+lvlhalo
        do jcr=xstartr(2)-lvlhalo,xendr(2)+lvlhalo
            do kcr=1,nxmr
                tpdvr(kcr,jcr,icr) = sal(kcr,jcr,icr)
            end do
            if (SfixS==1) then
                tpdvr(0,jcr,icr) = 2.0*salbp(1,jcr,icr) - sal(1,jcr,icr)
            else
                tpdvr(0,jcr,icr) = sal(1,jcr,icr)
            end if
            if (SfixN==1) then
                tpdvr(nxr,jcr,icr) = 2.0*saltp(1,jcr,icr) - sal(nxmr,jcr,icr)
            else
                tpdvr(nxr,jcr,icr) = sal(nxmr,jcr,icr)
            end if
        end do
    end do

    do icr=xstartr(3)-1,xendr(3)
        do jcr=xstartr(2)-1,xendr(2)
            do kcr=0,nxmr
                
                qv3 = tpdvr(kcr-1:kcr+2,jcr-1:jcr+2,icr-1:icr+2)
                do ic=max(krangr(icr),1),min(krangr(icr+1)-1,nzm)
                    qv2(:,:) = qv3(:,:,1)*czsalc(1,ic) + qv3(:,:,2)*czsalc(2,ic) &
                             + qv3(:,:,3)*czsalc(3,ic) + qv3(:,:,4)*czsalc(4,ic)
                    do jc=max(jrangr(jcr),1),min(jrangr(jcr+1)-1,nym)
                        qv1(:) = qv2(:,1)*cysalc(1,jc) + qv2(:,2)*cysalc(2,jc) &
                               + qv2(:,3)*cysalc(3,jc) + qv2(:,4)*cysalc(4,jc)
                        do kc=max(irangr(kcr),1),min(irangr(kcr+1)-1,nxm)
                            salc(kc,jc,ic) = sum(qv1(1:4)*cxsalc(1:4,kc))
                        end do
                    end do
                end do
                
            enddo
        enddo
    enddo

    !CJH Not appropriate for staggered grid
    ! do ic=xstart(3),xend(3)
    !  do jc=xstart(2),xend(2)
    !     salc(1,jc,ic) = 0.d0
    !  enddo
    ! enddo

    return
end subroutine InterpSalMgrd