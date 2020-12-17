subroutine CreateSalStencil

    use param
    use input_grids
    use mgrd_arrays
    implicit none
    integer :: jc,kc,ic,i,j,k
    integer :: icr, jcr, kcr

    real xxl(-2:nxr+2), yyl(-1:nyr), zzl(-1:nzr)  !CS xxl extended for str+uni
    real xxs(-1:nxro+1), yys(-1:nyro+1), zzs(-1:nzro+1)
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
    IF(nxmro.eq.nxmr .AND. istr3ro.eq.istr3r) then

        irangs(0) = 1
        irangc(0) = 1
        do ic=1,nxmro  
            irangs(ic) = ic
            irangc(ic) = ic
        enddo
        irangs(nxro) = nxro
        irangc(nxro) = nxro

    ELSE

        irangs(0) = 1
        irangc(0) = 1
        do ic=1,nxmro
            do i=1,nxmr
                if(xmr(i).lt.xmro(ic) .and. xmr(i+1).ge.xmro(ic))then
                    irangs(ic) = i+1
                endif
            enddo
            do i=1,nxmr
                if(xcr(i).lt.xcro(ic) .and. xcr(i+1).ge.xcro(ic))then
                    irangc(ic) = i+1
                endif
            enddo
        enddo
        irangs(nxro) = nxr
        irangc(1) = 1
        irangc(nxro) = nxr

    ENDIF 

    IF(nymro.eq.nymr) then
         
        jrangs(0) = 1
        jrangc(0) = 1
        do jc=1,nymro
         jrangs(jc) = jc
         jrangc(jc) = jc
        enddo
        jrangs(nyro) = nyr
        jrangc(nyro) = nyr
    
    ELSE

        jrangs(0) = 1
        jrangc(0) = 1
        do jc=1,nymro
          do j=1,nymr
            if(ymr(j).lt.ymro(jc) .and. ymr(j+1).ge.ymro(jc))then
              jrangs(jc) = j+1
            endif
          enddo
          do j=1,nymr
            if(ycr(j).lt.ycro(jc) .and. ycr(j+1).ge.ycro(jc))then
              jrangc(jc) = j+1
            endif
          enddo
        enddo
        jrangs(nyro) = nyr
        jrangc(1) = 1
        jrangc(nyro) = nyr

    ENDIF 

    IF(nzmro.eq.nzmr) then

        krangs(0) = 1
        krangc(0) = 1
        do kc=1,nzmro
         krangs(kc) = kc
         krangc(kc) = kc
        enddo
        krangs(nzro) = nzr
        krangc(nzro) = nzr

    ELSE

        krangs(0) = 1
        krangc(0) = 1
        do kc=1,nzmro
          do k=1,nzmr
            if(zmr(k).lt.zmro(kc) .and. zmr(k+1).ge.zmro(kc))then
              krangs(kc) = k+1
            endif
          enddo
          do k=1,nzmr
            if(zcr(k).lt.zcro(kc) .and. zcr(k+1).ge.zcro(kc))then
              krangc(kc) = k+1
            endif
          enddo
        enddo
        krangs(nzro) = nzr
        krangc(1) = 1
        krangc(nzro) = nzr

    ENDIF

    !=======================================
    !--- Interpolate Coefficients for S ----
    !=======================================
    !-- Set-up large (l) and small (s) arrays
    zzl(1:nzmr) = zmr(1:nzmr)
    zzl(0) = 2.d0*zzl(1) - zzl(2)
    zzl(nzr) = 2.d0*zzl(nzmr) - zzl(nzmr-1)
    zzs(1:nzmro) = zmro(1:nzmro)
    zzs(0) = 2.d0*zzs(1) - zzs(2)
    zzs(-1) = 2.d0*zzs(0) - zzs(1)
    zzs(nzro) = 2.d0*zzs(nzmro) - zzs(nzmro-1)
    zzs(nzro+1) = 2.d0*zzs(nzro) - zzs(nzmro)

    yyl(1:nymr) = ymr(1:nymr)
    yyl(0) = 2.d0*yyl(1) - yyl(2)
    yyl(nyr) = 2.d0*yyl(nymr) - yyl(nymr-1)
    yys(1:nymro) = ymro(1:nymro)
    yys(0) = 2.d0*yys(1) - yys(2)
    yys(-1) = 2.d0*yys(0) - yys(1)
    yys(nyro) = 2.d0*yys(nymro) - yys(nymro-1)
    yys(nyro+1) = 2.d0*yys(nyro) - yys(nymro)

    xxl(0) = xcr(1)
    xxl(1:nxmr) = xmr(1:nxmr)
    xxl(nxr) = xcr(nxr)
    xxs(0) = xcro(1)
    xxs(1:nxmro) = xmro(1:nxmro)
    xxs(nxro) = xcro(nxro)

    !-- Get weights for gradient calculations
    !-- Now construct Hermite basis function in x | cxrs
    do ic=0,nxmro
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
                cxrs(2,icr)=lxp  !h00-h10-h11*dlp/(dlp+dlc)
                cxrs(3,icr)=lxm  !h10+h01+h11*(dlp-dlc)/dlp
                cxrs(4,icr)=0.d0 !h11*dlc*dlc/dlp/(dlp+dlc)
            enddo
        elseif(ic.eq.nxmro)then
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
                cxrs(2,icr)=lxp  !h00-h11+h10*(dlc-dlm)/dlm
                cxrs(3,icr)=lxm  !h01+h11+h10*dlm/(dlm+dlc)
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
    do jc=0,nymro
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
    do kc=0,nzmro
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

    return
end subroutine CreateSalStencil