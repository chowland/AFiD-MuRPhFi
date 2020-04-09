      subroutine InterpVelMgrd

      use param
      use local_arrays, only: vx,vy,vz
      use mgrd_arrays
      use mpih
      use decomp_2d
      use AuxiliaryRoutines
      implicit none
       
      integer  :: ic,jc,kc, ip,jp,kp, icr,jcr,kcr
      integer  :: jc0,jcr0 !ic0,icr0, , ipr,jpr

      real,dimension(4,4,4) :: qv3 
      real,dimension(4,4) :: qv2
      real,dimension(4) :: qv1

      !real,dimension(1:m1mr,1:m2mr) :: dwloc,dwbot,dwall
      !real,allocatable,dimension(:,:,:) :: tpdv, tpdvr!, tpvr
      real,allocatable,dimension(:,:) :: vyxzc ! q1yzc,

      !real udx1r, udx2r
      real  :: lxa ! dl1q, dl2q, dlf1, dlf2, lza

      !-- Allocate temporary arrays for velocities and gradients
      real,allocatable,dimension(:,:,:) :: vyrT
      real,allocatable,dimension(:,:,:) :: tpdvryT !,tpdvrzT
      !call AllocateReal3DArray(tpvr,1,nxmr,xstartr(2),xendr(2),xstartr(3),xendr(3))

      tpdv(:,:,:) = 0.d0
      tpdvr(:,:,:) = 0.d0
      !tpvr = 0.d0

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
         tpdv(kc,jc,ic)=(vx(kp,jc,ic)-vx(kc,jc,ic))*udx3m(kc)
        enddo

        !-- Boundary points, enforce continuity
        tpdv( 0,jc,ic)=-(vy( 0,jp,ic)-vy( 0,jc,ic))*dy &
                       -(vz( 0,jc,ip)-vz( 0,jc,ic))*dz 
        tpdv(nx,jc,ic)=-(vy(nx,jp,ic)-vy(nx,jc,ic))*dy &
                       -(vz(nx,jc,ip)-vz(nx,jc,ic))*dz 

       enddo
      enddo

      call update_halo(tpdv,2)   !CS Are the corners updated? Might need to check.

      !-- Now interpolate gradients to refined grid
      do ic=xstart(3)-1,xend(3)
       do jc=xstart(2)-1,xend(2)
        do kc=0,nxm

         qv3=tpdv(kc-1:kc+2,jc-1:jc+2,ic-1:ic+2)

         do icr=max(krangs(ic),1),min(krangs(ic+1)-1,nzmr)
          qv2(:,:) = qv3(:,:,1)*czrs(1,icr)+qv3(:,:,2)*czrs(2,icr) &
                    +qv3(:,:,3)*czrs(3,icr)+qv3(:,:,4)*czrs(4,icr)
          do jcr=max(jrangs(jc),1),min(jrangs(jc+1)-1,nymr)
           qv1(:) = qv2(:,1)*cyrs(1,jcr)+qv2(:,2)*cyrs(2,jcr) &
                   +qv2(:,3)*cyrs(3,jcr)+qv2(:,4)*cyrs(4,jcr)
           do kcr=max(irangs(kc),1),min(irangs(kc+1)-1,nxmr)
            tpdvr(kcr,jcr,icr) = sum(qv1*cxrs(:,kcr))
           enddo
          enddo
         enddo

        enddo
       enddo
      enddo

      !-- Integrate vxr using interpolated gradients
      vxr(1,:,:)=0.d0
      do icr=xstartr(3),xendr(3)
       do jcr=xstartr(2),xendr(2)
        do kcr=1,nxmr
         lxa = xcr(kcr+1)-xcr(kcr)
         vxr(kcr+1,jcr,icr) = vxr(kcr,jcr,icr)+tpdvr(kcr,jcr,icr)*lxa
        enddo
       enddo
      enddo
      
      !-- Enforce zero net flux in x !TODO

!!!!       !call DestroyReal3DArray(tpdv)

!=========================================================
!     Interpolation of vy

      !-- Interpolate dvy/dy to cell center after transpose
      do ic=xstart(3),xend(3)
       do jc=xstart(2),xend(2)
        jp=jc+1

        !-- Interior points
        do kc=1,nxm
         kp=kc+1
         tpdv(kc,jc,ic)=(vy(kc,jp,ic)-vy(kc,jc,ic))*dy
        enddo

       enddo
      enddo

      call update_halo(tpdv,2)   !CS Are the corners updated? Might need to check.

      !-- Now interpolate gradients to refined grid
      do ic=xstart(3)-1,xend(3)
       do jc=xstart(2)-1,xend(2)
        do kc=0,nxm

         qv3=tpdv(kc-1:kc+2,jc-1:jc+2,ic-1:ic+2)

         do icr=max(krangs(ic),1),min(krangs(ic+1)-1,nzmr)
          qv2(:,:) = qv3(:,:,1)*czrs(1,icr)+qv3(:,:,2)*czrs(2,icr) &
                    +qv3(:,:,3)*czrs(3,icr)+qv3(:,:,4)*czrs(4,icr)
          do jcr=max(jrangs(jc),1),min(jrangs(jc+1)-1,nymr)
           qv1(:) = qv2(:,1)*cyrs(1,jcr)+qv2(:,2)*cyrs(2,jcr) &
                   +qv2(:,3)*cyrs(3,jcr)+qv2(:,4)*cyrs(4,jcr)
           do kcr=max(irangs(kc),1),min(irangs(kc+1)-1,nxmr)
            tpdvr(kcr,jcr,icr) = sum(qv1*cxrs(:,kcr))
           enddo
          enddo
         enddo

        enddo
       enddo
      enddo

      allocate(vyrT(ystartr(1):yendr(1),ystartr(2):yendr(2),ystartr(3):yendr(3)))
      allocate(tpdvryT(ystartr(1):yendr(1),ystartr(2):yendr(2),ystartr(3):yendr(3)))
      !write(*,*)'ss',nrank,ystartr(1:3)
      !write(*,*)'ee',nrank,yendr(1:3)
      !write(*,*) 'x', size(tpdvr,1), size(tpdvr,2), size(tpdvr,3)
      !write(*,*) 'y', size(tpdvryT,1), size(tpdvryT,2), size(tpdvryT,3)
      !-- Transpose x to y and integrate
      call transpose_x_to_y(vyr,vyrT)
      call transpose_x_to_y(tpdvr,tpdvryT)

      !-- Interpolate vyr at an arbitrary x-z plane | tpvr | 1st guess
      jc0 = 1
      jcr0 = 1
      do ic=ystartr(3),yendr(3)
       do kc=ystartr(1)-1,yendr(1)+1
         vyxzc(kc,ic)=vy(kc,jc0,ic)
       enddo
      enddo
      !TODO halo_updates in other pencil directions XXX
      !-- Integrate along each y-line 

      !-- Construct 2nd guess of vyr (vyr) at an arbitrary x-z plane
      !-- Integrate along each y-line 

      !-- Average of two integrations
      !-- Transpose y to x

!=========================================================
!     Interpolation of vz


!=========================================================
      call update_halo(vxr,lvlhalo)
      !call update_halo(vyr)
      !call update_halo(vzr)

      return
      end subroutine InterpVelMgrd
