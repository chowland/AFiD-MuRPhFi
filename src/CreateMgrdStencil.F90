!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: MgrdAuxiliaryRoutines.F90                      !
!    CONTAINS: subroutines CreateMgrdStencil,             !
!               InterpVelMgrd, Interp*Mgrd                !
!                                                         ! 
!    PURPOSE: Auxiliary routines for multigrid interp     !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine CreateMgrdStencil

      use param
      use mgrd_arrays
      use mpih
      implicit none
      integer :: jc,kc,ic,i,j,k
      integer :: icr, jcr, kcr

      real xxl(0:nxr), yyl(-1:nyr), zzl(-1:nzr)
      real xxs(-1:nx+1), yys(-1:ny+1), zzs(-1:nz+1)
      real h00, h01, h10, h11
      real lxm,lxp, lym,lyp, lzm,lzp, lxa,lya,lza
      real dlc, dlm, dlp

      !===================================================
      !  First create index look-up arrays for fine grid
      !===================================================
      !-- irangs, jrangs, krangs are at staggered location
      !-- irangc, jrangc, krangc are at cell-face location
      IF(nxm.eq.nxmr) then

      irangs(0) = 1
      irangc(0) = 1
      do ic=1,nxm  
        irangs(ic) = ic
        irangc(ic) = ic !CHECK
      enddo
      irangs(nx) = nx
      irangc(nx) = nx

      ELSE

      irangs(0) = 1
      irangc(0) = 1
      do ic=1,nxm
        do i=irangs(ic-1),nxmr
          if(xmr(i).lt.xm(ic) .and. xmr(i+1).ge.xm(ic))then
            irangs(ic) = i+1
          endif
        enddo
        do i=irangc(ic-1),nxmr
          if(xcr(i).lt.xc(ic) .and. xcr(i+1).ge.xc(ic))then
            irangc(ic) = i+1
          endif
        enddo
      enddo
      irangs(nx) = nxr
      irangc(1) = 1
      irangc(nx) = nxr

      ENDIF 

      IF(nym.eq.nymr) then
         
      jrangs(0) = 1
      jrangc(0) = 1
      do jc=1,nym
       jrangs(jc) = jc
       jrangc(jc) = jc !CHECK
      enddo
      jrangs(ny) = nyr
      jrangc(ny) = nyr
      
      ELSE

      jrangs(0) = 1
      jrangc(0) = 1
      do jc=1,nym
        do j=jrangs(jc-1),nymr
          if(ymr(j).lt.ym(jc) .and. ymr(j+1).ge.ym(jc))then
            jrangs(jc) = j+1
          endif
        enddo
        do j=jrangc(jc-1),nymr
          if(ycr(j).lt.yc(jc) .and. ycr(j+1).ge.yc(jc))then
            jrangc(jc) = j+1
          endif
        enddo
      enddo
      jrangs(ny) = nyr
      jrangc(1) = 1
      jrangc(ny) = nyr

      ENDIF 

      IF(nzm.eq.nzmr) then

      krangs(0) = 1
      krangc(0) = 1
      do kc=1,nzm
       krangs(kc) = kc
       krangc(kc) = kc
      enddo
      krangs(nz) = nzr
      krangc(nz) = nzr

      ELSE

      krangs(0) = 1
      krangc(0) = 1
      do kc=1,nzm
        do k=krangs(kc-1),nzmr
          if(zmr(k).lt.zm(kc) .and. zmr(k+1).ge.zm(kc))then
            krangs(kc) = k+1
          endif
        enddo
        do k=krangc(kc-1),nzmr
          if(zcr(k).lt.zc(kc) .and. zcr(k+1).ge.zc(kc))then
            krangc(kc) = k+1
          endif
        enddo
      enddo
      krangs(nz) = nzr
      krangc(1) = 1
      krangc(nz) = nzr

      ENDIF

! #ifdef DEBUG
      if(ismaster)then
        open(202,file='outputdir/itp_rangs.txt',status='unknown')
        do ic=0,nx
          write(202,*)ic, irangs(ic)
        enddo
        write(202,*)' '
        write(202,*)' '
        do jc=0,ny
          write(202,*)jc, jrangs(jc)
        enddo
        write(202,*)' '
        write(202,*)' '
        do kc=0,nz
          write(202,*)kc, krangs(kc)
        enddo
        close(202)
        open(203,file='outputdir/itp_rangc.txt',status='unknown')
        do ic=0,nx
          write(203,*)ic, irangc(ic)
        enddo
        write(203,*)' '
        write(203,*)' '
        do jc=0,ny
          write(203,*)jc, jrangc(jc)
        enddo
        write(203,*)' '
        write(203,*)' '
        do kc=0,nz
          write(203,*)kc, krangc(kc)
        enddo
        close(203)
      endif
! #endif

      !=======================================
      !--- Interpolate Coefficients for Vx ---
      !=======================================
      !-- Set-up large (l) and small (s) arrays
      zzl(1:nzmr) = zmr(1:nzmr)
      zzl(0) = 2.d0*zzl(1) - zzl(2)
      zzl(nzr) = 2.d0*zzl(nzmr) - zzl(nzmr-1)
      zzs(1:nzm) = zm(1:nzm)
      zzs(0) = 2.d0*zzs(1) - zzs(2)
      zzs(-1) = 2.d0*zzs(0) - zzs(1)
      zzs(nz) = 2.d0*zzs(nzm) - zzs(nzm-1)
      zzs(nz+1) = 2.d0*zzs(nz) - zzs(nzm)

      yyl(1:nymr) = ymr(1:nymr)
      yyl(0) = 2.d0*yyl(1) - yyl(2)
      yyl(nyr) = 2.d0*yyl(nymr) - yyl(nymr-1)
      yys(1:nym) = ym(1:nym)
      yys(0) = 2.d0*yys(1) - yys(2)
      yys(-1) = 2.d0*yys(0) - yys(1)
      yys(ny) = 2.d0*yys(nym) - yys(nym-1)
      yys(ny+1) = 2.d0*yys(ny) - yys(nym)

      xxl(1:nxr) = xcr(1:nxr)
      xxs(1:nx) = xc(1:nx)

      !-- Get weights for gradient calculations
      !-- Now construct Hermite basis function in x | cxvx
      do ic=1,nxm
         if(ic.eq.1)then
            dlc = xxs(ic+1)-xxs(ic)
            dlp = xxs(ic+2)-xxs(ic+1)
            lxa = 1.d0/dlc
            do icr=max(irangc(ic),1),min(irangc(ic+1)-1,nxmr)
               lxm = (xxl(icr) - xxs(ic))*lxa
               lxp = 1.d0 - lxm
               h00=(1.d0+2.d0*lxm)*lxp*lxp
               h10=lxm*lxp*lxp
               h01=(1.d0+2.d0*lxp)*lxm*lxm
               h11=-lxp*lxm*lxm
               cxvx(1,icr)=0.d0
               cxvx(2,icr)=h00-h10-h11*dlp/(dlp+dlc)
               cxvx(3,icr)=h10+h01+h11*(dlp-dlc)/dlp
               cxvx(4,icr)=h11*dlc*dlc/dlp/(dlp+dlc)
            enddo
         elseif(ic.eq.nxm)then
            dlc = xxs(ic+1)-xxs(ic)
            dlm = xxs(ic)-xxs(ic-1)
            lxa = 1.d0/dlc
            do icr=max(irangc(ic),1),min(irangc(ic+1)-1,nxmr)
               lxm = (xxl(icr) - xxs(ic))*lxa
               lxp = 1.d0 - lxm
               h00=(1.d0+2.d0*lxm)*lxp*lxp
               h10=lxm*lxp*lxp
               h01=(1.d0+2.d0*lxp)*lxm*lxm
               h11=-lxp*lxm*lxm
               cxvx(1,icr)=-h10*dlc*dlc/dlm/(dlc+dlm)
               cxvx(2,icr)=h00-h11+h10*(dlc-dlm)/dlm
               cxvx(3,icr)=h01+h11+h10*dlm/(dlm+dlc)
               cxvx(4,icr)=0.d0
            enddo
         else
            dlc = xxs(ic+1)-xxs(ic)
            dlm = xxs(ic)-xxs(ic-1)
            dlp = xxs(ic+2)-xxs(ic+1)
            lxa = 1.d0/dlc
            !Original do icr=max((ic-1)*3+1,1),min(ic*3,nxmr)
            do icr=max(irangc(ic),1),min(irangc(ic+1)-1,nxmr)
               lxm = (xxl(icr) - xxs(ic))*lxa
               lxp = 1.d0 - lxm
               h00=(1.d0+2.d0*lxm)*lxp*lxp
               h10=lxm*lxp*lxp
               h01=(1.d0+2.d0*lxp)*lxm*lxm
               h11=-lxp*lxm*lxm
               cxvx(1,icr)=-h10*dlc*dlc/dlm/(dlc+dlm)
               cxvx(2,icr)=h00-h11*dlp/(dlp+dlc)+h10*(dlc-dlm)/dlm
               cxvx(3,icr)=h01+h10*dlm/(dlm+dlc)+h11*(dlp-dlc)/dlp
               cxvx(4,icr)=h11*dlc*dlc/dlp/(dlp+dlc)
            enddo
         endif
      enddo
      !-- Now construct Hermite basis function in y | cyvx
      do jc=0,nym
         dlc = yys(jc+1)-yys(jc)
         dlm = yys(jc)-yys(jc-1)
         dlp = yys(jc+2)-yys(jc+1)
         lya=1.d0/dlc
         do jcr=max(jrangs(jc),1),min(jrangs(jc+1)-1,nymr)
            lym = (yyl(jcr) - yys(jc))*lya
            lyp = 1.d0 - lym
            h00=(1.d0+2.d0*lym)*lyp*lyp
            h10=lym*lyp*lyp
            h01=(1.d0+2.d0*lyp)*lym*lym
            h11=-lyp*lym*lym
            cyvx(1,jcr)=-h10*dlc*dlc/dlm/(dlc+dlm)
            cyvx(2,jcr)=h00-h11*dlp/(dlp+dlc)+h10*(dlc-dlm)/dlm
            cyvx(3,jcr)=h01+h10*dlm/(dlm+dlc)+h11*(dlp-dlc)/dlp
            cyvx(4,jcr)=h11*dlc*dlc/dlp/(dlp+dlc)
         enddo
      enddo
      !-- Now construct Hermite basis function in z | czvx
      do kc=0,nzm
         dlc = zzs(kc+1)-zzs(kc)
         dlm = zzs(kc)-zzs(kc-1)
         dlp = zzs(kc+2)-zzs(kc+1)
         lza=1.d0/dlc
         do kcr=max(krangs(kc),1),min(krangs(kc+1)-1,nzmr)
            lzm = (zzl(kcr) - zzs(kc))*lza
            lzp = 1.d0 - lzm
            h00=(1.d0+2.d0*lzm)*lzp*lzp
            h10=lzm*lzp*lzp
            h01=(1.d0+2.d0*lzp)*lzm*lzm
            h11=-lzp*lzm*lzm
            czvx(1,kcr)=-h10*dlc*dlc/dlm/(dlc+dlm)
            czvx(2,kcr)=h00-h11*dlp/(dlp+dlc)+h10*(dlc-dlm)/dlm
            czvx(3,kcr)=h01+h10*dlm/(dlm+dlc)+h11*(dlp-dlc)/dlp
            czvx(4,kcr)=h11*dlc*dlc/dlp/(dlp+dlc)
         enddo
      enddo

      !=======================================
      !--- Interpolate Coefficients for Vy ---
      !=======================================
      !-- Set-up large (l) and small (s) arrays
      zzl(1:nzmr) = zmr(1:nzmr)
      zzl(0) = 2.d0*zzl(1) - zzl(2)
      zzl(nzr) = 2.d0*zzl(nzmr) - zzl(nzmr-1)
      zzs(1:nzm) = zm(1:nzm)
      zzs(0) = 2.d0*zzs(1) - zzs(2)
      zzs(-1) = 2.d0*zzs(0) - zzs(1)
      zzs(nz) = 2.d0*zzs(nzm) - zzs(nzm-1)
      zzs(nz+1) = 2.d0*zzs(nz) - zzs(nzm)

      yyl(1:nyr) = ycr(1:nyr)
      yyl(0) = 2.d0*yyl(1) - yyl(2)
      yys(1:ny) = yc(1:ny)
      yys(0) = 2.d0*yys(1) - yys(2)
      yys(-1) = 2.d0*yys(0) - yys(1)
      yys(ny+1) = 2.d0*yys(ny) - yys(nym)


      xxl(0) = 0.d0 ! xcr(1)
      xxl(1:nxmr) = xmr(1:nxmr)
      xxl(nxr) = alx3 ! xcr(nxr)
      xxs(0) = 0.d0 ! xc(1)
      xxs(1:nxm) = xm(1:nxm)
      xxs(nx) = alx3 ! xc(nx)

      !-- Get weights for gradient calculations
      !-- Now construct Hermite basis function in x | cxvy
      do ic=0,nxm
         if(ic.eq.0)then
            dlc = xxs(ic+1)-xxs(ic)
            dlp = xxs(ic+2)-xxs(ic+1)
            lxa = 1.d0/dlc
            do icr=max(irangs(ic),1),min(irangs(ic+1)-1,nxmr)
               lxm = (xxl(icr) - xxs(ic))*lxa
               lxp = 1.d0 - lxm
               h00=(1.d0+2.d0*lxm)*lxp*lxp
               h10=lxm*lxp*lxp
               h01=(1.d0+2.d0*lxp)*lxm*lxm
               h11=-lxp*lxm*lxm
               cxvy(1,icr)=0.d0
               cxvy(2,icr)=h00-h11*dlp/(dlp+dlc)-dble(inslws)*h10
               cxvy(3,icr)=h01+h11*(dlp-dlc)/dlp+dble(inslws)*h10
               cxvy(4,icr)=h11*dlc*dlc/dlp/(dlp+dlc)
            enddo
         elseif(ic.eq.nxm)then
            dlc = xxs(ic+1)-xxs(ic)
            dlm = xxs(ic)-xxs(ic-1)
            lxa = 1.d0/dlc
            do icr=max(irangs(ic),1),min(irangs(ic+1)-1,nxmr)
               lxm = (xxl(icr) - xxs(ic))*lxa
               lxp = 1.d0 - lxm
               h00=(1.d0+2.d0*lxm)*lxp*lxp
               h10=lxm*lxp*lxp
               h01=(1.d0+2.d0*lxp)*lxm*lxm
               h11=-lxp*lxm*lxm
               cxvy(1,icr)=-h10*dlc*dlc/dlm/(dlc+dlm)
               cxvy(2,icr)=h00+h10*(dlc-dlm)/dlm-dble(inslwn)*h11
               cxvy(3,icr)=h01+h10*dlm/(dlm+dlc)+dble(inslwn)*h11
               cxvy(4,icr)=0.d0
            enddo
         else
            dlc = xxs(ic+1)-xxs(ic)
            dlm = xxs(ic)-xxs(ic-1)
            dlp = xxs(ic+2)-xxs(ic+1)
            lxa = 1.d0/dlc
            do icr=max(irangs(ic),1),min(irangs(ic+1)-1,nxmr)
               lxm = (xxl(icr) - xxs(ic))*lxa
               lxp = 1.d0 - lxm
               h00=(1.d0+2.d0*lxm)*lxp*lxp
               h10=lxm*lxp*lxp
               h01=(1.d0+2.d0*lxp)*lxm*lxm
               h11=-lxp*lxm*lxm
               cxvy(1,icr)=-h10*dlc*dlc/dlm/(dlc+dlm)
               cxvy(2,icr)=h00-h11*dlp/(dlp+dlc)+h10*(dlc-dlm)/dlm
               cxvy(3,icr)=h01+h10*dlm/(dlm+dlc)+h11*(dlp-dlc)/dlp
               cxvy(4,icr)=h11*dlc*dlc/dlp/(dlp+dlc)
            enddo
         endif
      enddo
      !-- Now construct Hermite basis function in y | cyvy
      do jc=1,nym
         dlc = yys(jc+1)-yys(jc)
         dlm = yys(jc)-yys(jc-1)
         dlp = yys(jc+2)-yys(jc+1)
         lya=1.d0/dlc
         do jcr=max(jrangc(jc),1),min(jrangc(jc+1),nymr)
            lym = (yyl(jcr) - yys(jc))*lya
            lyp = 1.d0 - lym
            h00=(1.d0+2.d0*lym)*lyp*lyp
            h10=lym*lyp*lyp
            h01=(1.d0+2.d0*lyp)*lym*lym
            h11=-lyp*lym*lym
            cyvy(1,jcr)=-h10*dlc*dlc/dlm/(dlc+dlm)
            cyvy(2,jcr)=h00-h11*dlp/(dlp+dlc)+h10*(dlc-dlm)/dlm
            cyvy(3,jcr)=h01+h10*dlm/(dlm+dlc)+h11*(dlp-dlc)/dlp
            cyvy(4,jcr)=h11*dlc*dlc/dlp/(dlp+dlc)
         enddo
      enddo

      !-- Now construct Hermite basis function in z | czvy
      do kc=0,nzm
        dlc = zzs(kc+1)-zzs(kc)
        dlm = zzs(kc)-zzs(kc-1)
        dlp = zzs(kc+2)-zzs(kc+1)
        lza=1.d0/dlc
        do kcr=max(krangs(kc),1),min(krangs(kc+1)-1,nzmr)
         lzm = (zzl(kcr) - zzs(kc))*lza
         lzp = 1.d0 - lzm
         h00=(1.d0+2.d0*lzm)*lzp*lzp
         h10=lzm*lzp*lzp
         h01=(1.d0+2.d0*lzp)*lzm*lzm
         h11=-lzp*lzm*lzm
         czvy(1,kcr)=-h10*dlc*dlc/dlm/(dlc+dlm)
         czvy(2,kcr)=h00-h11*dlp/(dlp+dlc)+h10*(dlc-dlm)/dlm
         czvy(3,kcr)=h01+h10*dlm/(dlm+dlc)+h11*(dlp-dlc)/dlp
         czvy(4,kcr)=h11*dlc*dlc/dlp/(dlp+dlc)
        enddo
      enddo

      !=======================================
      !--- Interpolate Coefficients for Vz ---
      !=======================================
      !-- Set-up large (l) and small (s) arrays
      zzl(1:nzr) = zcr(1:nzr)
      zzl(0) = 2.d0*zzl(1) - zzl(2)
      zzs(1:nz) = zc(1:nz)
      zzs(0) = 2.d0*zzs(1) - zzs(2)
      zzs(-1) = 2.d0*zzs(0) - zzs(1)
      zzs(nz+1) = 2.d0*zzs(nz) - zzs(nzm)

      yyl(1:nymr) = ymr(1:nymr)
      yyl(0) = 2.d0*yyl(1) - yyl(2)
      yyl(nyr) = 2.d0*yyl(nymr) - yyl(nymr-1)
      yys(1:nym) = ym(1:nym)
      yys(0) = 2.d0*yys(1) - yys(2)
      yys(-1) = 2.d0*yys(0) - yys(1)
      yys(ny) = 2.d0*yys(nym) - yys(nym-1)
      yys(ny+1) = 2.d0*yys(ny) - yys(nym)

      xxl(0) = 0.d0
      xxl(1:nxmr) = xmr(1:nxmr)
      xxl(nxr) = alx3
      xxs(0) = 0.d0
      xxs(1:nxm) = xm(1:nxm)
      xxs(nx) = alx3

      !-- Get weights for gradient calculations
      !-- Now construct Hermite basis function in x | cxvz
      do ic=0,nxm
         if(ic.eq.0)then
            dlc = xxs(ic+1)-xxs(ic)
            dlp = xxs(ic+2)-xxs(ic+1)
            lxa = 1.d0/dlc
            do icr=max(irangs(ic),1),min(irangs(ic+1)-1,nxmr)
               lxm = (xxl(icr) - xxs(ic))*lxa
               lxp = 1.d0 - lxm
               h00=(1.d0+2.d0*lxm)*lxp*lxp
               h10=lxm*lxp*lxp
               h01=(1.d0+2.d0*lxp)*lxm*lxm
               h11=-lxp*lxm*lxm
               cxvz(1,icr)=0.d0
               cxvz(2,icr)=h00-h11*dlp/(dlp+dlc)-dble(inslws)*h10
               cxvz(3,icr)=h01+h11*(dlp-dlc)/dlp+dble(inslws)*h10
               cxvz(4,icr)=h11*dlc*dlc/dlp/(dlp+dlc)
            enddo
         elseif(ic.eq.nxm)then
            dlc = xxs(ic+1)-xxs(ic)
            dlm = xxs(ic)-xxs(ic-1)
            lxa = 1.d0/dlc
            do icr=max(irangs(ic),1),min(irangs(ic+1)-1,nxmr)
               lxm = (xxl(icr) - xxs(ic))*lxa
               lxp = 1.d0 - lxm
               h00=(1.d0+2.d0*lxm)*lxp*lxp
               h10=lxm*lxp*lxp
               h01=(1.d0+2.d0*lxp)*lxm*lxm
               h11=-lxp*lxm*lxm
               cxvz(1,icr)=-h10*dlc*dlc/dlm/(dlc+dlm)
               cxvz(2,icr)=h00+h10*(dlc-dlm)/dlm-dble(inslwn)*h11
               cxvz(3,icr)=h01+h10*dlm/(dlm+dlc)+dble(inslwn)*h11
               cxvz(4,icr)=0.d0
            enddo
         else
            dlc = xxs(ic+1)-xxs(ic)
            dlm = xxs(ic)-xxs(ic-1)
            dlp = xxs(ic+2)-xxs(ic+1)
            lxa = 1.d0/dlc
            do icr=max(irangs(ic),1),min(irangs(ic+1)-1,nxmr)
               lxm = (xxl(icr) - xxs(ic))*lxa
               lxp = 1.d0 - lxm
               h00=(1.d0+2.d0*lxm)*lxp*lxp
               h10=lxm*lxp*lxp
               h01=(1.d0+2.d0*lxp)*lxm*lxm
               h11=-lxp*lxm*lxm
               cxvz(1,icr)=-h10*dlc*dlc/dlm/(dlc+dlm)
               cxvz(2,icr)=h00-h11*dlp/(dlp+dlc)+h10*(dlc-dlm)/dlm
               cxvz(3,icr)=h01+h10*dlm/(dlm+dlc)+h11*(dlp-dlc)/dlp
               cxvz(4,icr)=h11*dlc*dlc/dlp/(dlp+dlc)
            enddo
         endif
      enddo
      !-- Now construct Hermite basis function in y | cyvz
      do jc=0,nym
         dlc = yys(jc+1)-yys(jc)
         dlm = yys(jc)-yys(jc-1)
         dlp = yys(jc+2)-yys(jc+1)
         lya=1.d0/dlc
         do jcr=max(jrangs(jc),1),min(jrangs(jc+1),nymr)
            lym = (yyl(jcr) - yys(jc))*lya
            lyp = 1.d0 - lym
            h00=(1.d0+2.d0*lym)*lyp*lyp
            h10=lym*lyp*lyp
            h01=(1.d0+2.d0*lyp)*lym*lym
            h11=-lyp*lym*lym
            cyvz(1,jcr)=-h10*dlc*dlc/dlm/(dlc+dlm)
            cyvz(2,jcr)=h00-h11*dlp/(dlp+dlc)+h10*(dlc-dlm)/dlm
            cyvz(3,jcr)=h01+h10*dlm/(dlm+dlc)+h11*(dlp-dlc)/dlp
            cyvz(4,jcr)=h11*dlc*dlc/dlp/(dlp+dlc)
         enddo
      enddo

      !-- Now construct Hermite basis function in z | czvz
      do kc=1,nzm
        dlc = zzs(kc+1)-zzs(kc)
        dlm = zzs(kc)-zzs(kc-1)
        dlp = zzs(kc+2)-zzs(kc+1)
        lza=1.d0/dlc
        do kcr=max(krangc(kc),1),min(krangc(kc+1)-1,nzmr)
         lzm = (zzl(kcr) - zzs(kc))*lza
         lzp = 1.d0 - lzm
         h00=(1.d0+2.d0*lzm)*lzp*lzp
         h10=lzm*lzp*lzp
         h01=(1.d0+2.d0*lzp)*lzm*lzm
         h11=-lzp*lzm*lzm
         czvz(1,kcr)=-h10*dlc*dlc/dlm/(dlc+dlm)
         czvz(2,kcr)=h00-h11*dlp/(dlp+dlc)+h10*(dlc-dlm)/dlm
         czvz(3,kcr)=h01+h10*dlm/(dlm+dlc)+h11*(dlp-dlc)/dlp
         czvz(4,kcr)=h11*dlc*dlc/dlp/(dlp+dlc)
        enddo
      enddo

      !===================================================
      !--- For second order interpolation of du_i/dx_i ---
      !    (Helps make the interpolated velocity
      !    field divergence free.)
      !===================================================
      !-- Set-up large (l) and small (s) arrays
      zzl(1:nzmr) = zmr(1:nzmr)
      zzl(0) = 2.d0*zzl(1) - zzl(2)
      zzl(nzr) = 2.d0*zzl(nzmr) - zzl(nzmr-1)
      zzs(1:nzm) = zm(1:nzm)
      zzs(0) = 2.d0*zzs(1) - zzs(2)
      zzs(-1) = 2.d0*zzs(0) - zzs(1)
      zzs(nz) = 2.d0*zzs(nzm) - zzs(nzm-1)
      zzs(nz+1) = 2.d0*zzs(nz) - zzs(nzm)

      yyl(1:nymr) = ymr(1:nymr)
      yyl(0) = 2.d0*yyl(1) - yyl(2)
      yyl(nyr) = 2.d0*yyl(nymr) - yyl(nymr-1)
      yys(1:nym) = ym(1:nym)
      yys(0) = 2.d0*yys(1) - yys(2)
      yys(-1) = 2.d0*yys(0) - yys(1)
      yys(ny) = 2.d0*yys(nym) - yys(nym-1)
      yys(ny+1) = 2.d0*yys(ny) - yys(nym)

      xxl(0) = xcr(1)
      xxl(1:nxmr) = xmr(1:nxmr)
      xxl(nxr) = xcr(nxr)
      xxs(0) = xc(1)
      xxs(1:nxm) = xm(1:nxm)
      xxs(nx) = xc(nx)

      !-- Get weights for gradient calculations
      !-- Now construct Hermite basis function in x | cxrs
      do ic=0,nxm
         if(ic.eq.0)then
            dlc = xxs(ic+1)-xxs(ic)
            dlp = xxs(ic+2)-xxs(ic+1)
            lxa = 1.d0/dlc
            do icr=max(irangs(ic),1),min(irangs(ic+1)-1,nxmr)
               lxm = (xxl(icr) - xxs(ic))*lxa
               lxp = 1.d0 - lxm
               h00=(1.d0+2.d0*lxm)*lxp*lxp
               h10=lxm*lxp*lxp
               h01=(1.d0+2.d0*lxp)*lxm*lxm
               h11=-lxp*lxm*lxm
               cxrs(1,icr)=0.d0
               cxrs(2,icr)=lzp  !h00-h10-h11*dlp/(dlp+dlc)
               cxrs(3,icr)=lzm  !h10+h01+h11*(dlp-dlc)/dlp
               cxrs(4,icr)=0.d0 !h11*dlc*dlc/dlp/(dlp+dlc)
            enddo
         elseif(ic.eq.nxm)then
            dlc = xxs(ic+1)-xxs(ic)
            dlm = xxs(ic)-xxs(ic-1)
            lxa = 1.d0/dlc
            do icr=max(irangs(ic),1),min(irangs(ic+1)-1,nxmr)
               lxm = (xxl(icr) - xxs(ic))*lxa
               lxp = 1.d0 - lxm
               h00=(1.d0+2.d0*lxm)*lxp*lxp
               h10=lxm*lxp*lxp
               h01=(1.d0+2.d0*lxp)*lxm*lxm
               h11=-lxp*lxm*lxm
               cxrs(1,icr)=0.d0 !-h10*dlc*dlc/dlm/(dlc+dlm)
               cxrs(2,icr)=lzp  !h00-h11+h10*(dlc-dlm)/dlm
               cxrs(3,icr)=lzm  !h01+h11+h10*dlm/(dlm+dlc)
               cxrs(4,icr)=0.d0
            enddo
         else
            dlc = xxs(ic+1)-xxs(ic)
            dlm = xxs(ic)-xxs(ic-1)
            dlp = xxs(ic+2)-xxs(ic+1)
            lxa = 1.d0/dlc
            do icr=max(irangs(ic),1),min(irangs(ic+1)-1,nxmr)
               lxm = (xxl(icr) - xxs(ic))*lxa
               lxp = 1.d0 - lxm
               h00=(1.d0+2.d0*lxm)*lxp*lxp
               h10=lxm*lxp*lxp
               h01=(1.d0+2.d0*lxp)*lxm*lxm
               h11=-lxp*lxm*lxm
               cxrs(1,icr)=-h10*dlc*dlc/dlm/(dlc+dlm)
               cxrs(2,icr)=h00-h11*dlp/(dlp+dlc)+h10*(dlc-dlm)/dlm
               cxrs(3,icr)=h01+h10*dlm/(dlm+dlc)+h11*(dlp-dlc)/dlp
               cxrs(4,icr)=h11*dlc*dlc/dlp/(dlp+dlc)
            enddo
         endif
      enddo
      !-- Now construct Hermite basis function in y | cyrs
      do jc=0,nym
         dlc = yys(jc+1)-yys(jc)
         dlm = yys(jc)-yys(jc-1)
         dlp = yys(jc+2)-yys(jc+1)
         lya=1.d0/dlc
         do jcr=max(jrangs(jc),1),min(jrangs(jc+1)-1,nymr)
            lym = (yyl(jcr) - yys(jc))*lya
            lyp = 1.d0 - lym
            h00=(1.d0+2.d0*lym)*lyp*lyp
            h10=lym*lyp*lyp
            h01=(1.d0+2.d0*lyp)*lym*lym
            h11=-lyp*lym*lym
            cyrs(1,jcr)=-h10*dlc*dlc/dlm/(dlc+dlm)
            cyrs(2,jcr)=h00-h11*dlp/(dlp+dlc)+h10*(dlc-dlm)/dlm
            cyrs(3,jcr)=h01+h10*dlm/(dlm+dlc)+h11*(dlp-dlc)/dlp
            cyrs(4,jcr)=h11*dlc*dlc/dlp/(dlp+dlc)
         enddo
      enddo
      !-- Now construct Hermite basis function in z | czrs
      do kc=0,nzm
         dlc = zzs(kc+1)-zzs(kc)
         dlm = zzs(kc)-zzs(kc-1)
         dlp = zzs(kc+2)-zzs(kc+1)
         lza=1.d0/dlc
         do kcr=max(krangs(kc),1),min(krangs(kc+1)-1,nzmr)
            lzm = (zzl(kcr) - zzs(kc))*lza
            lzp = 1.d0 - lzm
            h00=(1.d0+2.d0*lzm)*lzp*lzp
            h10=lzm*lzp*lzp
            h01=(1.d0+2.d0*lzp)*lzm*lzm
            h11=-lzp*lzm*lzm
            czrs(1,kcr)=-h10*dlc*dlc/dlm/(dlc+dlm)
            czrs(2,kcr)=h00-h11*dlp/(dlp+dlc)+h10*(dlc-dlm)/dlm
            czrs(3,kcr)=h01+h10*dlm/(dlm+dlc)+h11*(dlp-dlc)/dlp
            czrs(4,kcr)=h11*dlc*dlc/dlp/(dlp+dlc)
         enddo
      enddo

      !==============================================
      !--- Interpolate Coefficients for Coarse Vx ---
      !==============================================
      !-- Set-up large (l) and small (s) arrays
      zzl(1:nzmr) = zmr(1:nzmr)
      zzl(0) = 2.d0*zzl(1) - zzl(2)
      zzl(-1) = 2.d0*zzl(0) - zzl(1)
      zzl(nzr) = 2.d0*zzl(nzmr) - zzl(nzmr-1)
      zzs(1:nzm) = zm(1:nzm)
      zzs(0) = 2.d0*zzs(1) - zzs(2)
      zzs(-1) = 2.d0*zzs(0) - zzs(1)
      zzs(nz) = 2.d0*zzs(nzm) - zzs(nzm-1)
      zzs(nz+1) = 2.d0*zzs(nz) - zzs(nzm)

      yyl(1:nymr) = ymr(1:nymr)
      yyl(0) = 2.d0*yyl(1) - yyl(2)
      yyl(-1) = 2.d0*yyl(0) - yyl(1)
      yyl(nyr) = 2.d0*yyl(nymr) - yyl(nymr-1)
      yys(1:nym) = ym(1:nym)
      yys(0) = 2.d0*yys(1) - yys(2)
      yys(-1) = 2.d0*yys(0) - yys(1)
      yys(ny) = 2.d0*yys(nym) - yys(nym-1)
      yys(ny+1) = 2.d0*yys(ny) - yys(nym)

      xxl(1:nxr) = xcr(1:nxr)
      xxs(1:nx) = xc(1:nx)

      !-- Get weights for gradient calculations
      !-- Now construct Hermite basis function in x | cxvxc
      do ic=1,nxm
         if(ic.eq.1)then
            cxvxc(1,ic)=0.d0!-h10*dlc*dlc/dlm/(dlc+dlm)
            cxvxc(2,ic)=1.d0!h00-h11*dlp/(dlp+dlc)+h10*(dlc-dlm)/dlm
            cxvxc(3,ic)=0.d0!h01+h10*dlm/(dlm+dlc)+h11*(dlp-dlc)/dlp
            cxvxc(4,ic)=0.d0!h11*dlc*dlc/dlp/(dlp+dlc)
         elseif(ic.eq.nxm)then
            cxvxc(1,ic)=0.d0!-h10*dlc*dlc/dlm/(dlc+dlm)
            cxvxc(2,ic)=0.d0!h00-h11*dlp/(dlp+dlc)+h10*(dlc-dlm)/dlm
            cxvxc(3,ic)=1.d0!h01+h10*dlm/(dlm+dlc)+h11*(dlp-dlc)/dlp
            cxvxc(4,ic)=0.d0!h11*dlc*dlc/dlp/(dlp+dlc)
         else
            icr = irangc(ic)-1
            dlc = xxl(icr+1)-xxl(icr)
            dlm = xxl(icr)-xxl(icr-1)
            dlp = xxl(icr+2)-xxl(icr+1)
            lxa = 1.d0/dlc
            lxm = (xxs(ic) - xxl(icr))*lxa
            lxp = 1.d0 - lxm

            h00=(1.d0+2.d0*lxm)*lxp*lxp
            h10=lxm*lxp*lxp
            h01=(1.d0+2.d0*lxp)*lxm*lxm
            h11=-lxp*lxm*lxm
            cxvxc(1,ic)=-h10*dlc*dlc/dlm/(dlc+dlm)
            cxvxc(2,ic)=h00-h11*dlp/(dlp+dlc)+h10*(dlc-dlm)/dlm
            cxvxc(3,ic)=h01+h10*dlm/(dlm+dlc)+h11*(dlp-dlc)/dlp
            cxvxc(4,ic)=h11*dlc*dlc/dlp/(dlp+dlc)
         endif
      enddo
      !-- Now construct Hermite basis function in y | cyvxc
      do jc=0,nym
         jcr = jrangs(jc)-1
         dlc = yyl(jcr+1)-yyl(jcr)
         dlm = yyl(jcr)-yyl(jcr-1)
         dlp = yyl(jcr+2)-yyl(jcr+1)
         lya = 1.d0/dlc
         lym = (yys(jc) - yyl(jcr))*lya
         lyp = 1.d0 - lym

         h00=(1.d0+2.d0*lym)*lyp*lyp
         h10=lym*lyp*lyp
         h01=(1.d0+2.d0*lyp)*lym*lym
         h11=-lyp*lym*lym
         cyvxc(1,jc)=-h10*dlc*dlc/dlm/(dlc+dlm)
         cyvxc(2,jc)=h00-h11*dlp/(dlp+dlc)+h10*(dlc-dlm)/dlm
         cyvxc(3,jc)=h01+h10*dlm/(dlm+dlc)+h11*(dlp-dlc)/dlp
         cyvxc(4,jc)=h11*dlc*dlc/dlp/(dlp+dlc)
      enddo
      !-- Now construct Hermite basis function in z | czvxc
      do kc=0,nzm
         kcr = krangs(kc)-1
         dlc = zzl(kcr+1)-zzl(kcr)
         dlm = zzl(kcr)-zzl(kcr-1)
         dlp = zzl(kcr+2)-zzl(kcr+1)
         lza = 1.d0/dlc
         lzm = (zzs(kc) - zzl(kcr))*lza
         lzp = 1.d0 - lzm

         h00=(1.d0+2.d0*lzm)*lzp*lzp
         h10=lzm*lzp*lzp
         h01=(1.d0+2.d0*lzp)*lzm*lzm
         h11=-lzp*lzm*lzm
         czvxc(1,kc)=-h10*dlc*dlc/dlm/(dlc+dlm)
         czvxc(2,kc)=h00-h11*dlp/(dlp+dlc)+h10*(dlc-dlm)/dlm
         czvxc(3,kc)=h01+h10*dlm/(dlm+dlc)+h11*(dlp-dlc)/dlp
         czvxc(4,kc)=h11*dlc*dlc/dlp/(dlp+dlc)
      enddo

      return
      end subroutine CreateMgrdStencil

