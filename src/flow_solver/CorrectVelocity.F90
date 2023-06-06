!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: CorrectVelocity.F90                            !
!    CONTAINS: subroutine CorrectVelocity                 !
!                                                         ! 
!    PURPOSE: Update velocities with the pressure         !
!     correction to enforce incompresibility              !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine CorrectVelocity
    use param
    use local_arrays, only: vy,vx,dphhalo,vz,temp
    use afid_salinity, only: sal
    use decomp_2d, only: xstart,xend,xstartr,xendr
    use mpih
    implicit none
    integer :: jc,jm,kc,km,ic,im,kmid
    real    :: usukm,udy,udz,locdph
    real, dimension(nxm) :: vxbar, Gshape
    real, dimension(nxmr):: GshapeR
    real :: vybulk, vzbulk, Tbulk, Sbulk, idx, vz_target, lam, lfac

    udy = al*dt*dy
    udz = al*dt*dz

!$OMP  PARALLEL DO &
!$OMP   DEFAULT(none) &
!$OMP   SHARED(vz,vy,vx,dphhalo,udz,udy,udx3c) &
!$OMP   SHARED(xstart,xend,nxm,kmv,dt,al) &
!$OMP   PRIVATE(ic,jc,kc) &
!$OMP   PRIVATE(im,jm,km,usukm,locdph)
    do ic=xstart(3),xend(3)
        im=ic-1
        do jc=xstart(2),xend(2)
            jm=jc-1
            do kc=1,nxm
                km=kmv(kc)
                usukm = al*dt*udx3c(kc)
                locdph=dphhalo(kc,jc,ic)
                vx(kc,jc,ic)=vx(kc,jc,ic)- &
                    (locdph-dphhalo(km,jc,ic))*usukm
                vy(kc,jc,ic)=vy(kc,jc,ic)- &
                    (locdph-dphhalo(kc,jm,ic))*udy
                vz(kc,jc,ic)=vz(kc,jc,ic)- &
                    (locdph-dphhalo(kc,jc,im))*udz
            enddo
       enddo
    enddo
!$OMP END PARALLEL DO

    !CJH Prescribe mean volume flux
    !Treat dPdz input as a desired Re_b
    vz_target = dPdz/ren
    if (dPdz/=0) then
        vzbulk = 0.d0
        lam = sqrt(2.0*ren/al/dt)
        do ic=xstart(3),xend(3)
            do jc=xstart(2),xend(2)
                do kc=1,nxm
                    idx = 1/udx3m(kc)
                    vzbulk = vzbulk + vz(kc,jc,ic)*idx
                end do
            end do
        end do

        call MpiAllSumRealScalar(vzbulk)
        vzbulk = vzbulk/nym/nzm

        kmid = nxm/2
        lfac = (1.0 + 2.0/lam*(exp(-lam/2.0) - 1.0))
        do kc=1,kmid
            Gshape(kc) = (1.0 - exp(-lam*xm(kc)))/lfac
        end do
        do kc=kmid+1,nxm
            Gshape(kc) = (1.0 - exp(lam*(xm(kc) - 1.0)))/lfac
        end do
        do ic=xstart(3),xend(3)
            do jc=xstart(2),xend(2)
                do kc=1,nxm
                    vz(kc,jc,ic) = vz(kc,jc,ic) + (vz_target - vzbulk)*Gshape(kc)
                end do
            end do
        end do
    end if

    !CJH Remove mean mass flux
    if((.not.melt .and. .not.phasefield) .and. (gAxis==2 .and. inslwN==1)) then

        vxbar(:)=0.d0
        vybulk = 0.d0; vzbulk = 0.d0; Tbulk = 0.d0
        do ic=xstart(3),xend(3)
            do jc=xstart(2),xend(2)
                do kc=2,nxm
                    vxbar(kc) = vxbar(kc) + vx(kc,jc,ic)
                end do
                do kc=1,nxm
                    idx = 1/udx3m(kc)
                    vybulk = vybulk + vy(kc,jc,ic)*idx
                    vzbulk = vzbulk + vz(kc,jc,ic)*idx
                    Tbulk = Tbulk + temp(kc,jc,ic)*idx
                end do
            end do
        end do

        call MpiAllSumReal1D(vxbar,nxm)
        call MpiAllSumRealScalar(vybulk)
        call MpiAllSumRealScalar(vzbulk)
        call MpiAllSumRealScalar(Tbulk)

        if (salinity) then
            Sbulk = 0.d0
            do ic=xstartr(3),xendr(3)
                do jc=xstartr(2),xendr(2)
                    do kc=1,nxmr
                        idx = 1/udx3mr(kc)
                        Sbulk = Sbulk + sal(kc,jc,ic)*idx
                    end do
                end do
            end do

            call MpiAllSumRealScalar(Sbulk)

            Sbulk = Sbulk/nymr/nzmr
        end if
        Tbulk = Tbulk/nym/nzm
        vzbulk = vzbulk/nym/nzm
        vybulk = vybulk/nym/nzm
        do kc=1,nxm
            vxbar(kc) = vxbar(kc)/nym/nzm
        end do
        
        lam = sqrt(2.0*ren/al/dt)
        kmid = nxm/2
        lfac = (1.0 + 2.0/lam*(exp(-lam/2.0) - 1.0))
        do kc=1,kmid
            Gshape(kc) = (1.0 - exp(-lam*xm(kc)))/lfac
        end do
        do kc=kmid+1,nxm
            Gshape(kc) = (1.0 - exp(lam*(xm(kc) - 1.0)))/lfac
        end do
        kmid = nxmr/2
        do kc=1,kmid
            GshapeR(kc) = (1.0 - exp(-lam*xmr(kc)))/lfac
        end do
        do kc=kmid+1,nxmr
            GshapeR(kc) = (1.0 - exp(lam*(xmr(kc) - 1.0)))/lfac
        end do

        do ic=xstart(3),xend(3)
            do jc=xstart(2),xend(2)
                do kc=2,nxm
                    vx(kc,jc,ic) = vx(kc,jc,ic) - vxbar(kc)
                end do
                do kc=1,nxm
                    vy(kc,jc,ic) = vy(kc,jc,ic) - vybulk*Gshape(kc)
                end do
            end do
        end do
        if (.not.phasefield) then
            do ic=xstart(3),xend(3)
                do jc=xstart(2),xend(2)
                    do kc=1,nxm
                        temp(kc,jc,ic) = temp(kc,jc,ic) - Tbulk*Gshape(kc)
                    end do
                end do
            end do
        end if
        if (dPdz.eq.0) then
            do ic=xstart(3),xend(3)
                do jc=xstart(2),xend(2)
                    do kc=1,nxm
                        vz(kc,jc,ic) = vz(kc,jc,ic) + (vz_target - vzbulk)*Gshape(kc)
                    end do
                end do
            end do
        end if
        if (salinity) then
            do ic=xstartr(3),xendr(3)
                do jc=xstartr(2),xendr(2)
                    do kc=1,nxmr
                        sal(kc,jc,ic) = sal(kc,jc,ic) - Sbulk*GshapeR(kc)
                    end do
                end do
            end do
        end if

    end if

    return
end