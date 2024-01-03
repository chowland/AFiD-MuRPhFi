subroutine InterpVelMgrd

    use param
    use local_arrays, only: vx,vy,vz
    use mgrd_arrays
    use afid_salinity, only: vxr, vyr, vzr
    use mpih
    use decomp_2d
    use AuxiliaryRoutines
    use HermiteInterpolations
    implicit none
    
    integer  :: ic,jc,kc, ip,jp,kp, icr,jcr,kcr
    integer  :: jc0,jcr0, ic0,icr0
    integer  :: comm_col,comm_row,ierror

    real,dimension(4,4) :: qv2
    real,dimension(4) :: qv1

    real,allocatable,dimension(:,:) :: vyxzc,vyxzr, vzxyc,vzxyr

    real  :: lxa, lya, lza

    !-- Allocate temporary arrays for velocities and gradients
    real,allocatable,dimension(:,:) :: dvyloc,dvybot, dvzloc,dvzbot

    tpdv(:,:,:) = 0.d0   ! Temporary gradient array - coarse
    tpdvr(:,:,:) = 0.d0  ! Temporary gradient array - refined
    vxr(:,:,:) = 0.d0
    vyr(:,:,:) = 0.d0
    vzr(:,:,:) = 0.d0

!=========================================================
!     Interpolation of vx

    !-- Interpolate dvx/dx to cell center after transpose
    do ic=xstart(3),xend(3)
        ip=ic+1
        do jc=xstart(2),xend(2)
            jp=jc+1
            !-- Interior points
            do kc=1,nxm
                kp=kc+1
                tpdv(kc,jc,ic) = (vx(kp,jc,ic) - vx(kc,jc,ic))*udx3m(kc)
            end do
        !-- Boundary points, enforce continuity
        !CJH Note indices 0 & nx are now beyond boundaries in CreateMgrdStencil
        tpdv( 0,jc,ic) = -tpdv(1,jc,ic)
        tpdv(nx,jc,ic) = -tpdv(nxm,jc,ic)
        end do
    end do

    call update_halo(tpdv,2)

    call interpolate_xyz_to_refined(tpdv,tpdvr(1:nxmr,:,:))

    !-- Integrate vxr using interpolated gradients
    vxr(1,:,:)=0.d0
    do icr=xstartr(3),xendr(3)
        do jcr=xstartr(2),xendr(2)
            do kcr=1,nxmr
                lxa = xcr(kcr+1) - xcr(kcr)
                vxr(kcr+1,jcr,icr) = vxr(kcr,jcr,icr) + tpdvr(kcr,jcr,icr)*lxa
            end do
        end do
    end do

!=========================================================
!     Interpolation of vy

    tpdv(:,:,:) = 0.d0
    tpdvr(:,:,:) = 0.d0

    !-- Interpolate dvy/dy to cell center after transpose
    do ic=xstart(3),xend(3)
        do jc=xstart(2),xend(2)
            jp=jc+1
            !-- Interior points
            do kc=1,nxm
                tpdv(kc,jc,ic) = (vy(kc,jp,ic) - vy(kc,jc,ic))*dy
            end do
            !-- Boundary points, enforce zero velocity / gradient
            tpdv( 0,jc,ic) = (1-2*inslwS)*tpdv(1,jc,ic)
            tpdv(nx,jc,ic) = (1-2*inslwN)*tpdv(nxm,jc,ic)
        end do
    end do

    call update_halo(tpdv,2)

    call interpolate_xyz_to_refined(tpdv, tpdvr(1:nxmr,:,:))

    !-- Interpolate vyr at an arbitrary x-z plane | tpvr | 1st guess
    call AllocateReal2DArray(vyxzc,-1,nxm+2,xstart(3)-2,xend(3)+2)
    call AllocateReal2DArray(vyxzr,1,nxmr,xstartr(3),xendr(3))
    jc0 = 1 !global index
    jcr0 = 1
    vyxzc(:,:)=0.d0
    if (jc0.ge.xstart(2).and.jc0.le.xend(2)) then
        do ic=xstart(3)-2,xend(3)+2
            do kc=1,nxm
                vyxzc(kc,ic) = vy(kc,jc0,ic)
            end do
            ! x boundaries
            vyxzc(0,ic) = (1-2*inslwS)*vyxzc(1,ic)
            vyxzc(nx,ic) = (1-2*inslwN)*vyxzc(nxm,ic)
        end do

        do ic=xstart(3)-1,xend(3)
            do kc=0,nxm

                qv2=vyxzc(kc-1:kc+2,ic-1:ic+2)

                do icr=max(krangs(ic),xstartr(3)),min(krangs(ic+1)-1,xendr(3))
                    qv1(:) = qv2(:,1)*czvy(1,icr) + qv2(:,2)*czvy(2,icr) &
                           + qv2(:,3)*czvy(3,icr) + qv2(:,4)*czvy(4,icr)
                    do kcr=max(irangs(kc),1),min(irangs(kc+1)-1,nxmr)
                        vyr(kcr,jcr0,icr) = sum(qv1*cxvy(:,kcr))
                    end do
                end do

            end do
        end do
    end if
    
    call MPI_CART_SUB(DECOMP_2D_COMM_CART_X,(/.true.,.false./),comm_col,ierror)

    lya = 1./dyr
    !-- Integrate along each y-line 
    do icr=xstartr(3),xendr(3)
        do jcr=xstartr(2),xendr(2)
            do kcr=1,nxmr
                vyr(kcr,jcr+1,icr) = vyr(kcr,jcr,icr) + tpdvr(kcr,jcr,icr)*lya
            enddo
        enddo
    enddo

    call AllocateReal2DArray(dvyloc,1,nxmr,xstartr(3),xendr(3)) !Local
    call AllocateReal2DArray(dvybot,1,nxmr,xstartr(3),xendr(3)) !Scan

    dvyloc=vyr(1:nxmr,xendr(2)+1,xstartr(3):xendr(3))
    dvybot=0.d0

    !-- Scan velocity in y direction
    call MPI_SCAN(dvyloc,dvybot,nxmr*xsizer(3),MDP,MPI_SUM,comm_col,ierror)
    do icr=xstartr(3),xendr(3)
        do jcr=xstartr(2),xendr(2)+1
            do kcr=1,nxmr
                vyr(kcr,jcr,icr) = vyr(kcr,jcr,icr) + dvybot(kcr,icr) - dvyloc(kcr,icr)
            end do
        end do
    end do
    call MPI_Comm_free(comm_col,ierror) !CS Probably can be avoided by creating comm once

!=========================================================
!     Interpolation of vz

    tpdv(:,:,:) = 0.d0
    tpdvr(:,:,:) = 0.d0

    !-- Interpolate dvz/dz to cell center after transpose
    do ic=xstart(3),xend(3)
        ip = ic + 1
        do jc=xstart(2),xend(2)
            !-- Interior points
            do kc=1,nxm
                tpdv(kc,jc,ic) = (vz(kc,jc,ip) - vz(kc,jc,ic))*dz
            enddo
            !-- Boundary points, enforce zero velocity / gradient
            tpdv( 0,jc,ic) = (1-2*inslwS)*tpdv(1,jc,ic)
            tpdv(nx,jc,ic) = (1-2*inslwN)*tpdv(nxm,jc,ic)
        end do
    end do

    call update_halo(tpdv,2)

    call interpolate_xyz_to_refined(tpdv, tpdvr(1:nxmr,:,:))

    !-- Interpolate vzr at an arbitrary x-y plane | tpvr | 1st guess
    call AllocateReal2DArray(vzxyc,-1,nxm+2,xstart(2)-2,xend(2)+2)
    call AllocateReal2DArray(vzxyr,1,nxmr,xstartr(2),xendr(2))
    ic0 = 1 !global index
    icr0 = 1
    vzxyc(:,:) = 0.d0
    if (ic0>=xstart(3) .and. ic0<xend(3)) then
        do jc=xstart(2)-2,xend(2)+2
            do kc=1,nxm
                vzxyc(kc,jc) = vz(kc,jc,ic0)
            end do
            ! Boundary points
            vzxyc( 0,jc) = (1-2*inslwS)*vzxyc(1,jc)
            vzxyc(nx,jc) = (1-2*inslwN)*vzxyc(nxm,jc)
        end do

        do jc=xstart(2)-1,xend(2)
            do kc=0,nxm

                qv2=vzxyc(kc-1:kc+2,jc-1:jc+2)

                do jcr=max(jrangs(jc),xstartr(2)),min(jrangs(jc+1)-1,xendr(2))
                    qv1(:) = qv2(:,1)*cyvz(1,jcr) + qv2(:,2)*cyvz(2,jcr) &
                           + qv2(:,3)*cyvz(3,jcr) + qv2(:,4)*cyvz(4,jcr)
                    do kcr=max(irangs(kc),1),min(irangs(kc+1)-1,nxmr)
                        vzr(kcr,jcr,icr0) = sum(qv1*cxvz(:,kcr))
                    end do
                end do

            end do
        end do
    end if

    call MPI_CART_SUB(DECOMP_2D_COMM_CART_X,(/.false.,.true./),comm_row,ierror)

    lza = 1./dzr
    !-- Integrate along each y-line 
    do icr=xstartr(3),xendr(3)
        do jcr=xstartr(2),xendr(2)
            do kcr=1,nxmr
            vzr(kcr,jcr,icr+1)=vzr(kcr,jcr,icr)+tpdvr(kcr,jcr,icr)*lza
            end do
        end do
    end do

    call AllocateReal2DArray(dvzloc,1,nxmr,xstartr(2),xendr(2)) !Local
    call AllocateReal2DArray(dvzbot,1,nxmr,xstartr(2),xendr(2)) !Scan

    dvzloc=vzr(1:nxmr,xstartr(2):xendr(2),xendr(3)+1)
    dvzbot=0.d0

    !-- Scan velocity in z direction
    call MPI_SCAN(dvzloc,dvzbot,nxmr*xsizer(2),MDP,MPI_SUM,comm_row,ierror)
    do icr=xstartr(3),xendr(3)+1
        do jcr=xstartr(2),xendr(2)
            do kcr=1,nxmr
            vzr(kcr,jcr,icr)=vzr(kcr,jcr,icr)+dvzbot(kcr,jcr)-dvzloc(kcr,jcr)
            end do
        end do
    end do
    call MPI_Comm_free(comm_row,ierror) !CS Probably can be avoided by creating comm once

!=========================================================

    return
end subroutine InterpVelMgrd