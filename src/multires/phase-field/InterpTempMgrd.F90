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
    use HermiteInterpolations, only: interpolate_xyz_to_refined
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
            tpdv(0,jc,ic) = 2.0*tempbp(1,jc,ic) - temp(1,jc,ic)
            tpdv(nx,jc,ic) = 2.0*temptp(1,jc,ic) - temp(nxm,jc,ic)
            do kc=1,nxm
                tpdv(kc,jc,ic) = temp(kc,jc,ic)
            end do
        end do
    end do

    call interpolate_xyz_to_refined(tpdv, tempr(1:nxmr,:,:))

    return
end subroutine InterpTempMgrd