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
    use mgrd_arrays, only: phi,phic,tpdvr
    use mpih
    use decomp_2d
    use AuxiliaryRoutines
    use HermiteInterpolations, only: interpolate_xyz_to_coarse, interpolate_xyz_to_coarse_fast
    implicit none

    integer  :: ic,jc,kc

    ! Interpolate to coarse grid here. A better option is to apply
    ! a box filter.
    phic(:,:,:) = 0.d0

    tpdvr(:,:,:) = 0.d0 ! Temporary array with extended range for interpolation
    do ic=xstartr(3)-lvlhalo,xendr(3)+lvlhalo
        do jc=xstartr(2)-lvlhalo,xendr(2)+lvlhalo
            tpdvr(0,jc,ic) = phi(1,jc,ic)
            do kc=1,nxmr
                tpdvr(kc,jc,ic) = phi(kc,jc,ic)
            end do
            tpdvr(nxr,jc,ic) = phi(nxmr,jc,ic)
        end do
    end do

    if ((xmr(1) < xm(1)) .and. (xmr(nxmr) > xm(nxm))) then
        call interpolate_xyz_to_coarse_fast(tpdvr, phic(1:nxm,:,:), "phi")
    else
        call interpolate_xyz_to_coarse(tpdvr, phic(1:nxm,:,:))
    end if
 
    return

end subroutine InterpPhiMgrd