      subroutine InterpVelMgrd

      use param
      use local_arrays, only: vx,vy,vz
      use mgrd_arrays
      use mpih
      use decomp_2d
      use AuxiliaryRoutines
      implicit none
       
      integer  :: ic,jc,kc, ip,jp,kp, icr,jcr,kcr
      integer  :: icc,jcc,kcc
      integer  :: jc0,jcr0, ic0,icr0
      integer  :: comm_col,comm_row,comm,ierror,chk

      real,dimension(4,4,4) :: qv3 
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
         tpdv(kc,jc,ic)=(vx(kp,jc,ic)-vx(kc,jc,ic))*udx3m(kc)
        enddo

        !-- Boundary points, enforce continuity
        !tpdv( 0,jc,ic)=-(vy( 0,jp,ic)-vy( 0,jc,ic))*dy &
        !               -(vz( 0,jc,ip)-vz( 0,jc,ic))*dz 
        tpdv( 0,jc,ic)=-(-vy( 1,jp,ic)+vy( 1,jc,ic))*dy &
                       -(-vz( 1,jc,ip)+vz( 1,jc,ic))*dz 
        tpdv(nx,jc,ic)=-(vy(nx,jp,ic)-vy(nx,jc,ic))*dy &
                       -(vz(nx,jc,ip)-vz(nx,jc,ic))*dz 

       enddo
      enddo

      call update_halo(tpdv,2)

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
            tpdvr(kcr,jcr,icr) = sum(qv1(1:4)*cxrs(1:4,kcr))
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
      !call DestroyReal3DArray(tpdv)

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
         tpdv(kc,jc,ic)=(vy(kc,jp,ic)-vy(kc,jc,ic))*dy
        enddo

        !-- Boundary points, enforce zero velocity
        tpdv(0,jc,ic)=(-vy(1,jp,ic)+vy(1,jc,ic))*dy
        tpdv(nx,jc,ic)=(-vy(nxm,jp,ic)+vy(nxm,jc,ic))*dy

       enddo
      enddo

      call update_halo(tpdv,2)

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

      !-- Interpolate vyr at an arbitrary x-z plane | tpvr | 1st guess
      call AllocateReal2DArray(vyxzc,-1,nxm+2,xstart(3)-2,xend(3)+2)
      call AllocateReal2DArray(vyxzr,1,nxmr,xstartr(3),xendr(3))
      jc0 = 1 !global index
      jcr0 = 1
      vyxzc(:,:)=0.d0
      if (jc0.ge.xstart(2).and.jc0.le.xend(2)) then
       do ic=xstart(3)-2,xend(3)+2
        do kc=1,nxm+1 !0,nxm+1
         vyxzc(kc,ic)=vy(kc,jc0,ic) !CS Halo updates can be optimised. Otherwise lvlhalo=2 required
        enddo
       enddo

       do ic=xstart(3)-1,xend(3)
        do kc=0,nxm

         qv2=vyxzc(kc-1:kc+2,ic-1:ic+2)

         do icr=max(krangs(ic),1),min(krangs(ic+1)-1,nzmr)
          qv1(:) = qv2(:,1)*czvy(1,icr)+qv2(:,2)*czvy(2,icr) &
                  +qv2(:,3)*czvy(3,icr)+qv2(:,4)*czvy(4,icr)
          do kcr=max(irangs(kc),1),min(irangs(kc+1)-1,nxmr)
           vyr(kcr,jcr0,icr) = sum(qv1*cxvy(:,kcr))
          enddo
         enddo

        enddo
       enddo
      endif
      
      call MPI_CART_SUB(DECOMP_2D_COMM_CART_X,(/.true.,.false./),comm_col,ierror)
      !call MPI_Comm_rank(comm, row_id, ierror)

      lya = 1./dyr
      !-- Integrate along each y-line 
      do icr=xstartr(3),xendr(3)
       do jcr=xstartr(2),xendr(2)
        do kcr=1,nxmr
         vyr(kcr,jcr+1,icr)=vyr(kcr,jcr,icr)+tpdvr(kcr,jcr,icr)*lya
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
         vyr(kcr,jcr,icr)=vyr(kcr,jcr,icr)+dvybot(kcr,icr)-dvyloc(kcr,icr)
        enddo
       enddo
      enddo
      call MPI_Comm_free(comm_col,ierror) !CS Probably can be avoided by creating comm once
      ! call MPI_BARRIER(comm,ierr)

      !CS Is the following necessary? Maybe not.
      !-- Construct 2nd guess of vyr (vyr) at an arbitrary x-z plane
      !-- Integrate along each y-line 
      !-- Average of two integrations

!=========================================================
!     Interpolation of vz

      tpdv(:,:,:) = 0.d0
      tpdvr(:,:,:) = 0.d0

      !-- Interpolate dvz/dz to cell center after transpose
      do ic=xstart(3),xend(3)
       ip=ic+1
       do jc=xstart(2),xend(2)

        !-- Interior points
        do kc=1,nxm
         tpdv(kc,jc,ic)=(vz(kc,jc,ip)-vz(kc,jc,ic))*dz
        enddo

        !-- Boundary points, enforce zero velocity
        tpdv(0,jc,ic)=(-vz(1,jc,ip)+vz(1,jc,ic))*dz
        tpdv(nx,jc,ic)=(-vz(nxm,jc,ip)+vz(nxm,jc,ic))*dz

       enddo
      enddo

      call update_halo(tpdv,2)   !CS Are the corners updated? Might need to check.

      !-- Now interpolate gradients to refined grid
      do ic=xstart(3)-1,xend(3)
       do jc=xstart(2)-1,xend(2)
        do kc=0,nxm

         qv3=tpdv(kc-1:kc+2,jc-1:jc+2,ic-1:ic+2)

         do icr=max(krangs(ic),1),min(krangs(ic+1)-1,nzmr)  !CS Is this correct?
          qv2(:,:) = qv3(:,:,1)*czrs(1,icr)+qv3(:,:,2)*czrs(2,icr) &
                    +qv3(:,:,3)*czrs(3,icr)+qv3(:,:,4)*czrs(4,icr)
          do jcr=max(jrangs(jc),1),min(jrangs(jc+1)-1,nymr) !CS Is this correct?
           qv1(:) = qv2(:,1)*cyrs(1,jcr)+qv2(:,2)*cyrs(2,jcr) &
                   +qv2(:,3)*cyrs(3,jcr)+qv2(:,4)*cyrs(4,jcr)
           do kcr=max(irangs(kc),1),min(irangs(kc+1)-1,nxmr) !CS Is this correct?
            tpdvr(kcr,jcr,icr) = sum(qv1*cxrs(:,kcr))
           enddo
          enddo
         enddo

        enddo
       enddo
      enddo
      
      !-- Interpolate vzr at an arbitrary x-y plane | tpvr | 1st guess
      call AllocateReal2DArray(vzxyc,-1,nxm+2,xstart(2)-2,xend(2)+2)
      call AllocateReal2DArray(vzxyr,1,nxmr,xstartr(2),xendr(2))
      ic0 = 1 !global index
      icr0 = 1
      vzxyc(:,:)=0.d0
      if (ic0.ge.xstart(3).and.ic0.le.xend(3)) then
       do jc=xstart(2)-2,xend(2)+2
        do kc=1,nxm+1 !0,nxm+1
         vzxyc(kc,jc)=vz(kc,jc,ic0) !CS Halo updates can be optimised. Otherwise lvlhalo=2 required
        enddo
       enddo

       do jc=xstart(2)-1,xend(2)
        do kc=0,nxm

         qv2=vzxyc(kc-1:kc+2,jc-1:jc+2)

         do jcr=max(jrangs(jc),1),min(jrangs(jc+1)-1,nymr)
          qv1(:) = qv2(:,1)*cyvz(1,jcr)+qv2(:,2)*cyvz(2,jcr) &
                  +qv2(:,3)*cyvz(3,jcr)+qv2(:,4)*cyvz(4,jcr)
          do kcr=max(irangs(kc),1),min(irangs(kc+1)-1,nxmr)
           vzr(kcr,jcr,icr0) = sum(qv1*cxvz(:,kcr))
          enddo
         enddo

        enddo
       enddo
      endif

      call MPI_CART_SUB(DECOMP_2D_COMM_CART_X,(/.false.,.true./),comm_row,ierror)
      !call MPI_Comm_rank(comm, row_id, ierror)

      lza = 1./dzr
      !-- Integrate along each y-line 
      do icr=xstartr(3),xendr(3)
       do jcr=xstartr(2),xendr(2)
        do kcr=1,nxmr
         vzr(kcr,jcr,icr+1)=vzr(kcr,jcr,icr)+tpdvr(kcr,jcr,icr)*lza
        enddo
       enddo
      enddo

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
        enddo
       enddo
      enddo
      call MPI_Comm_free(comm_row,ierror) !CS Probably can be avoided by creating comm once
      !call MPI_BARRIER(comm,ierr)

      !CS Is the following necessary? Maybe not.
      !-- Construct 2nd guess of vyr (vyr) at an arbitrary x-z plane
      !-- Integrate along each y-line 
      !-- Average of two integrations

!=========================================================
      !call update_halo(vxr,lvlhalo)
      !call update_halo(vyr,lvlhalo)
      !call update_halo(vzr,lvlhalo)

      return
      end subroutine InterpVelMgrd
