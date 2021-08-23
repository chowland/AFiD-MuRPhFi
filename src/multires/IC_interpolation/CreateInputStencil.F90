!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: MgrdAuxiliaryRoutines.F90                      !
!    CONTAINS: subroutines CreateMgrdStencil,             !
!               InterpVelMgrd, Interp*Mgrd                !
!                                                         ! 
!    PURPOSE: Auxiliary routines for multigrid interp     !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine CreateInputStencil

    use param
    use input_grids
    use mgrd_arrays
    implicit none
    integer :: jc,kc,ic,i,j,k
    integer :: icr, jcr, kcr

    real xxl(-2:nx+2), yyl(-1:ny), zzl(-1:nz)  !CS xxl extended for str+uni
    real xxs(-1:nxo+1), yys(-1:nyo+1), zzs(-1:nzo+1)
    real h00, h01, h10, h11
    real lxm,lxp, lym,lyp, lzm,lzp, lxa,lya,lza
    real dlc, dlm, dlp

    irangs(:) = 0
    irangc(:) = 0
    jrangs(:) = 0
    jrangc(:) = 0
    krangs(:) = 0
    krangc(:) = 0

    !===================================================
    !  First create index look-up arrays for fine grid
    !===================================================
    !-- irangs, jrangs, krangs are at staggered location
    !-- irangc, jrangc, krangc are at cell-face location
    IF(nxmo.eq.nxm .AND. istr3o.eq.istr3) then

        irangs(0) = 1
        irangc(0) = 1
        do ic=1,nxmo  
            irangs(ic) = ic
            irangc(ic) = ic
        enddo
        irangs(nxo) = nxo
        irangc(nxo) = nxo

    ELSE

        irangs(0) = 1
        irangc(0) = 1
        do ic=1,nxmo
            do i=1,nxm
                if(xm(i).lt.xmo(ic) .and. xm(i+1).ge.xmo(ic))then
                    irangs(ic) = i+1
                endif
            enddo
            do i=1,nxm
                if(xc(i).lt.xco(ic) .and. xc(i+1).ge.xco(ic))then
                    irangc(ic) = i+1
                endif
            enddo
        enddo
        irangs(nxo) = nx
        irangc(1) = 1
        irangc(nxo) = nx

    ENDIF 

    IF(nymo.eq.nym) then
         
        jrangs(0) = 1
        jrangc(0) = 1
        do jc=1,nymo
         jrangs(jc) = jc
         jrangc(jc) = jc
        enddo
        jrangs(nyo) = ny
        jrangc(nyo) = ny
    
    ELSE

        jrangs(0) = 1
        jrangc(0) = 1
        do jc=1,nymo
          do j=0,ny
            if(ym(j).lt.ymo(jc) .and. ym(j+1).ge.ymo(jc))then
              jrangs(jc) = j+1
            endif
          enddo
          do j=1,nym
            if(yc(j).lt.yco(jc) .and. yc(j+1).ge.yco(jc))then
              jrangc(jc) = j+1
            endif
          enddo
        enddo
        jrangs(nyo) = ny
        jrangc(1) = 1
        jrangc(nyo) = ny

    ENDIF 

    IF(nzmo.eq.nzm) then

        krangs(0) = 1
        krangc(0) = 1
        do kc=1,nzmo
         krangs(kc) = kc
         krangc(kc) = kc
        enddo
        krangs(nzo) = nz
        krangc(nzo) = nz

    ELSE

        krangs(0) = 1
        krangc(0) = 1
        do kc=1,nzmo
          do k=1,nzm
            if(zm(k).lt.zmo(kc) .and. zm(k+1).ge.zmo(kc))then
              krangs(kc) = k+1
            endif
          enddo
          do k=1,nzm
            if(zc(k).lt.zco(kc) .and. zc(k+1).ge.zco(kc))then
              krangc(kc) = k+1
            endif
          enddo
        enddo
        krangs(nzo) = nz
        krangc(1) = 1
        krangc(nzo) = nz

    ENDIF

    !=======================================
    !--- Interpolate Coefficients for Vx ---
    !=======================================
    !-- Set-up large (l) and small (s) arrays
    zzl(1:nzm) = zm(1:nzm)
    zzl(0) = 2.d0*zzl(1) - zzl(2)
    zzl(nz) = 2.d0*zzl(nzm) - zzl(nzm-1)
    zzs(1:nzmo) = zmo(1:nzmo)
    zzs(0) = 2.d0*zzs(1) - zzs(2)
    zzs(-1) = 2.d0*zzs(0) - zzs(1)
    zzs(nzo) = 2.d0*zzs(nzmo) - zzs(nzmo-1)
    zzs(nzo+1) = 2.d0*zzs(nzo) - zzs(nzmo)

    yyl(1:nym) = ym(1:nym)
    yyl(0) = 2.d0*yyl(1) - yyl(2)
    yyl(ny) = 2.d0*yyl(nym) - yyl(nym-1)
    yys(1:nymo) = ymo(1:nymo)
    yys(0) = 2.d0*yys(1) - yys(2)
    yys(-1) = 2.d0*yys(0) - yys(1)
    yys(nyo) = 2.d0*yys(nymo) - yys(nymo-1)
    yys(nyo+1) = 2.d0*yys(nyo) - yys(nymo)

    xxl(1:nx) = xc(1:nx)
    xxs(1:nxo) = xco(1:nxo)

    !-- Get weights for gradient calculations
    !-- Now construct Hermite basis function in x | cxvx
    do ic=1,nxmo
        if(ic.eq.1)then
            dlc = xxs(ic+1)-xxs(ic)
            dlp = xxs(ic+2)-xxs(ic+1)
            lxa = 1.d0/dlc
            do icr=max(irangc(ic),1),min(irangc(ic+1)-1,nxm)
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
        elseif(ic.eq.nxmo)then
            dlc = xxs(ic+1)-xxs(ic)
            dlm = xxs(ic)-xxs(ic-1)
            lxa = 1.d0/dlc
            do icr=max(irangc(ic),1),min(irangc(ic+1)-1,nxm)
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
            do icr=max(irangc(ic),1),min(irangc(ic+1)-1,nxm)
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
    do jc=0,nymo
        dlc = yys(jc+1)-yys(jc)
        dlm = yys(jc)-yys(jc-1)
        dlp = yys(jc+2)-yys(jc+1)
        lya=1.d0/dlc
        do jcr=max(jrangs(jc),1),min(jrangs(jc+1)-1,nym)
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
    do kc=0,nzmo
        dlc = zzs(kc+1)-zzs(kc)
        dlm = zzs(kc)-zzs(kc-1)
        dlp = zzs(kc+2)-zzs(kc+1)
        lza=1.d0/dlc
        do kcr=max(krangs(kc),1),min(krangs(kc+1)-1,nzm)
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
    zzl(1:nzm) = zm(1:nzm)
    zzl(0) = 2.d0*zzl(1) - zzl(2)
    zzl(nz) = 2.d0*zzl(nzm) - zzl(nzm-1)
    zzs(1:nzmo) = zmo(1:nzmo)
    zzs(0) = 2.d0*zzs(1) - zzs(2)
    zzs(-1) = 2.d0*zzs(0) - zzs(1)
    zzs(nzo) = 2.d0*zzs(nzmo) - zzs(nzmo-1)
    zzs(nzo+1) = 2.d0*zzs(nzo) - zzs(nzmo)

    yyl(1:ny) = yc(1:ny)
    yyl(0) = 2.d0*yyl(1) - yyl(2)
    yys(1:nyo) = yco(1:nyo)
    yys(0) = 2.d0*yys(1) - yys(2)
    yys(-1) = 2.d0*yys(0) - yys(1)
    yys(nyo+1) = 2.d0*yys(nyo) - yys(nymo)


    xxl(0) = 0.d0 ! xcr(1)
    xxl(1:nxm) = xm(1:nxm)
    xxl(nx) = alx3 ! xcr(nxr)
    xxs(0) = 0.d0 ! xc(1)
    xxs(1:nxmo) = xmo(1:nxmo)
    xxs(nxo) = alx3 ! xc(nx)

    !-- Get weights for gradient calculations
    !-- Now construct Hermite basis function in x | cxvy
    do ic=0,nxmo
        if(ic.eq.0)then
            dlc = xxs(ic+1)-xxs(ic)
            dlp = xxs(ic+2)-xxs(ic+1)
            lxa = 1.d0/dlc
            do icr=max(irangs(ic),1),min(irangs(ic+1)-1,nxm)
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
        elseif(ic.eq.nxmo)then
            dlc = xxs(ic+1)-xxs(ic)
            dlm = xxs(ic)-xxs(ic-1)
            lxa = 1.d0/dlc
            do icr=max(irangs(ic),1),min(irangs(ic+1)-1,nxm)
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
            do icr=max(irangs(ic),1),min(irangs(ic+1)-1,nxm)
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
    do jc=1,nymo
        dlc = yys(jc+1)-yys(jc)
        dlm = yys(jc)-yys(jc-1)
        dlp = yys(jc+2)-yys(jc+1)
        lya=1.d0/dlc
        do jcr=max(jrangc(jc),1),min(jrangc(jc+1),nym)
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
    do kc=0,nzmo
        dlc = zzs(kc+1)-zzs(kc)
        dlm = zzs(kc)-zzs(kc-1)
        dlp = zzs(kc+2)-zzs(kc+1)
        lza=1.d0/dlc
        do kcr=max(krangs(kc),1),min(krangs(kc+1)-1,nzm)
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
    zzl(1:nz) = zc(1:nz)
    zzl(0) = 2.d0*zzl(1) - zzl(2)
    zzs(1:nzo) = zco(1:nzo)
    zzs(0) = 2.d0*zzs(1) - zzs(2)
    zzs(-1) = 2.d0*zzs(0) - zzs(1)
    zzs(nzo+1) = 2.d0*zzs(nzo) - zzs(nzmo)

    yyl(1:nym) = ym(1:nym)
    yyl(0) = 2.d0*yyl(1) - yyl(2)
    yyl(ny) = 2.d0*yyl(nym) - yyl(nym-1)
    yys(1:nymo) = ymo(1:nymo)
    yys(0) = 2.d0*yys(1) - yys(2)
    yys(-1) = 2.d0*yys(0) - yys(1)
    yys(nyo) = 2.d0*yys(nymo) - yys(nymo-1)
    yys(nyo+1) = 2.d0*yys(nyo) - yys(nymo)

    xxl(0) = 0.d0
    xxl(1:nxm) = xm(1:nxm)
    xxl(nx) = alx3
    xxs(0) = 0.d0
    xxs(1:nxmo) = xmo(1:nxmo)
    xxs(nxo) = alx3

    !-- Get weights for gradient calculations
    !-- Now construct Hermite basis function in x | cxvz
    do ic=0,nxmo
        if(ic.eq.0)then
            dlc = xxs(ic+1)-xxs(ic)
            dlp = xxs(ic+2)-xxs(ic+1)
            lxa = 1.d0/dlc
            do icr=max(irangs(ic),1),min(irangs(ic+1)-1,nxm)
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
        elseif(ic.eq.nxmo)then
            dlc = xxs(ic+1)-xxs(ic)
            dlm = xxs(ic)-xxs(ic-1)
            lxa = 1.d0/dlc
            do icr=max(irangs(ic),1),min(irangs(ic+1)-1,nxm)
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
            do icr=max(irangs(ic),1),min(irangs(ic+1)-1,nxm)
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
    do jc=0,nymo
        dlc = yys(jc+1)-yys(jc)
        dlm = yys(jc)-yys(jc-1)
        dlp = yys(jc+2)-yys(jc+1)
        lya=1.d0/dlc
        do jcr=max(jrangs(jc),1),min(jrangs(jc+1),nym)
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
    do kc=1,nzmo
        dlc = zzs(kc+1)-zzs(kc)
        dlm = zzs(kc)-zzs(kc-1)
        dlp = zzs(kc+2)-zzs(kc+1)
        lza=1.d0/dlc
        do kcr=max(krangc(kc),1),min(krangc(kc+1)-1,nzm)
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
    zzl(1:nzm) = zm(1:nzm)
    zzl(0) = 2.d0*zzl(1) - zzl(2)
    zzl(nz) = 2.d0*zzl(nzm) - zzl(nzm-1)
    zzs(1:nzmo) = zmo(1:nzmo)
    zzs(0) = 2.d0*zzs(1) - zzs(2)
    zzs(-1) = 2.d0*zzs(0) - zzs(1)
    zzs(nzo) = 2.d0*zzs(nzmo) - zzs(nzmo-1)
    zzs(nzo+1) = 2.d0*zzs(nzo) - zzs(nzmo)

    yyl(1:nym) = ym(1:nym)
    yyl(0) = 2.d0*yyl(1) - yyl(2)
    yyl(ny) = 2.d0*yyl(nym) - yyl(nym-1)
    yys(1:nymo) = ymo(1:nymo)
    yys(0) = 2.d0*yys(1) - yys(2)
    yys(-1) = 2.d0*yys(0) - yys(1)
    yys(nyo) = 2.d0*yys(nymo) - yys(nymo-1)
    yys(nyo+1) = 2.d0*yys(nyo) - yys(nymo)

    xxl(0) = xc(1)
    xxl(1:nxm) = xm(1:nxm)
    xxl(nx) = xc(nx)
    xxs(0) = xco(1)
    xxs(1:nxmo) = xmo(1:nxmo)
    xxs(nxo) = xco(nxo)

    !-- Get weights for gradient calculations
    !-- Now construct Hermite basis function in x | cxrs
    do ic=0,nxmo
        if(ic.eq.0)then
            dlc = xxs(ic+1)-xxs(ic)
            dlp = xxs(ic+2)-xxs(ic+1)
            lxa = 1.d0/dlc
            do icr=max(irangs(ic),1),min(irangs(ic+1)-1,nxm)
                lxm = (xxl(icr) - xxs(ic))*lxa
                lxp = 1.d0 - lxm
                h00=(1.d0+2.d0*lxm)*lxp*lxp
                h10=lxm*lxp*lxp
                h01=(1.d0+2.d0*lxp)*lxm*lxm
                h11=-lxp*lxm*lxm
                cxrs(1,icr)=0.d0
                cxrs(2,icr)=lxp  !h00-h10-h11*dlp/(dlp+dlc)
                cxrs(3,icr)=lxm  !h10+h01+h11*(dlp-dlc)/dlp
                cxrs(4,icr)=0.d0 !h11*dlc*dlc/dlp/(dlp+dlc)
            enddo
        elseif(ic.eq.nxmo)then
            dlc = xxs(ic+1)-xxs(ic)
            dlm = xxs(ic)-xxs(ic-1)
            lxa = 1.d0/dlc
            do icr=max(irangs(ic),1),min(irangs(ic+1)-1,nxm)
                lxm = (xxl(icr) - xxs(ic))*lxa
                lxp = 1.d0 - lxm
                h00=(1.d0+2.d0*lxm)*lxp*lxp
                h10=lxm*lxp*lxp
                h01=(1.d0+2.d0*lxp)*lxm*lxm
                h11=-lxp*lxm*lxm
                cxrs(1,icr)=0.d0 !-h10*dlc*dlc/dlm/(dlc+dlm)
                cxrs(2,icr)=lxp  !h00-h11+h10*(dlc-dlm)/dlm
                cxrs(3,icr)=lxm  !h01+h11+h10*dlm/(dlm+dlc)
                cxrs(4,icr)=0.d0
            enddo
        else
            dlc = xxs(ic+1)-xxs(ic)
            dlm = xxs(ic)-xxs(ic-1)
            dlp = xxs(ic+2)-xxs(ic+1)
            lxa = 1.d0/dlc
            do icr=max(irangs(ic),1),min(irangs(ic+1)-1,nxm)
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
    do jc=0,nymo
        dlc = yys(jc+1)-yys(jc)
        dlm = yys(jc)-yys(jc-1)
        dlp = yys(jc+2)-yys(jc+1)
        lya=1.d0/dlc
        do jcr=max(jrangs(jc),1),min(jrangs(jc+1)-1,nym)
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
    do kc=0,nzmo
        dlc = zzs(kc+1)-zzs(kc)
        dlm = zzs(kc)-zzs(kc-1)
        dlp = zzs(kc+2)-zzs(kc+1)
        lza=1.d0/dlc
        do kcr=max(krangs(kc),1),min(krangs(kc+1)-1,nzm)
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

    return
end subroutine CreateInputStencil