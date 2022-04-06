subroutine InterpSalMgrd

    use param
    use mgrd_arrays, only: sal,salc,cxsalc,cysalc,czsalc,irangr,jrangr,krangr,tpdvr
    use mpih
    use decomp_2d
    use AuxiliaryRoutines
    use HermiteInterpolations, only: interpolate_xyz_to_coarse, interpolate_xyz_to_coarse_fast
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

    if ((xmr(1) < xm(1)) .and. (xmr(nxmr) > xm(nxm))) then
        call interpolate_xyz_to_coarse_fast(tpdvr, salc(1:nxm,:,:))
    else
        call interpolate_xyz_to_coarse(tpdvr, salc(1:nxm,:,:))
    end if

    !CJH Not appropriate for staggered grid
    ! do ic=xstart(3),xend(3)
    !  do jc=xstart(2),xend(2)
    !     salc(1,jc,ic) = 0.d0
    !  enddo
    ! enddo

    return
end subroutine InterpSalMgrd