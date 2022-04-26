!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: ExplicitTermsSal.F90                           !
!    CONTAINS: subroutine ExplicitTermsSal                !
!                                                         ! 
!    PURPOSE: Compute the non-linear terms associated to  !
!     the salinity.                                       !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ExplicitTermsSal
    use param
    use mgrd_arrays, only: vxr,vyr,vzr,sal,hsal,phi,rhsr
    use ibm_param, only: solidr
    use decomp_2d, only: xstartr,xendr
    implicit none
    integer :: jc,kc,ic
    integer :: km,kp,jm,jp,im,ip
    real    :: hsx,hsy,hsz,udyr,udzr
    real    :: udzrq,udyrq
    real    :: dyys,dzzs,pf_delta
    real, dimension(1:nxmr) :: sdx
    real    :: dpdsx, dpdsy, dpdsz, phfrac, aldt

    udzr=dzr*0.5d0
    udyr=dyr*0.5d0
    udzrq=dzqr/pecs
    udyrq=dyqr/pecs
    pf_delta = 1e-6     !0.01

    do kc=1,nxmr
        sdx(kc) = 0.5*dxr/g3rmr(kc)
    end do
    aldt = 1.0/al/dt

!$OMP  PARALLEL DO &
!$OMP   DEFAULT(none) &
!$OMP   SHARED(xstartr,xendr,vzr,vyr,vxr,nxmr) &
!$OMP   SHARED(kmv,kpv,am3sk,ac3sk,ap3sk,udzr) &
!$OMP   SHARED(udyr,udzrq,udyrq,udx3cr,sal,hsal) &
!$OMP   PRIVATE(ic,jc,kc,im,ip,km,kp,jm,jp) &
!$OMP   PRIVATE(hsx,hsy,hsz,dyys,dzzs)
    do ic=xstartr(3),xendr(3)
        im=ic-1
        ip=ic+1
        do jc=xstartr(2),xendr(2)
            jm=jc-1
            jp=jc+1

            kc = 1
            kp = 2
            hsx=(vxr(kp,jc,ic)*(sal(kp,jc,ic)+sal(kc,jc,ic))- &
                    vxr(kc,jc,ic)*2.d0*salbp(1,jc,ic))*udx3mr(kc)*0.5d0
            hsz=(vzr(kc,jc,ip)*(sal(kc,jc,ip)+sal(kc,jc,ic))- &
                    vzr(kc,jc,ic)*(sal(kc,jc,ic)+sal(kc,jc,im)) &
                )*udzr
            hsy=(vyr(kc,jp,ic)*(sal(kc,jp,ic)+sal(kc,jc,ic))- &
                    vyr(kc,jc,ic)*(sal(kc,jc,ic)+sal(kc,jm,ic)) &
                )*udyr
            ! dzzs=(sal(kc,jc,ip) - 2.0*sal(kc,jc,ic) + sal(kc,jc,im))*udzrq
            ! dyys=(sal(kc,jp,ic) - 2.0*sal(kc,jc,ic) + sal(kc,jm,ic))*udyrq
            if (solidr(kc,jc,ip)) then
                dzzs = (sal(kc,jc,im) - sal(kc,jc,ic))*udzrq
            elseif (solidr(kc,jc,im)) then
                dzzs = (sal(kc,jc,ip) - sal(kc,jc,ic))*udzrq
            else
                dzzs = (sal(kc,jc,ip) - 2.0*sal(kc,jc,ic) + sal(kc,jc,im))*udzrq
            end if
            if (solidr(kc,jp,ic)) then
                dyys = (sal(kc,jm,ic) - sal(kc,jc,ic))*udyrq
            elseif (solidr(kc,jm,ic)) then
                dyys = (sal(kc,jp,ic) - sal(kc,jc,ic))*udyrq
            else
                dyys = (sal(kc,jp,ic) - 2.0*sal(kc,jc,ic) + sal(kc,jm,ic))*udyrq
            end if
            hsal(kc,jc,ic) = -(hsx+hsy+hsz)+dyys+dzzs

            do kc=2,nxmr-1
                km=kc-1
                kp=kc+1
                !
                !    sal vxr term
                !
                !
                !                 d  sal q_x 
                !                -----------
                !                 d   x      
                !
                hsx=(vxr(kp,jc,ic)*(sal(kp,jc,ic)+sal(kc,jc,ic))- &
                    vxr(kc,jc,ic)*(sal(kc,jc,ic)+sal(km,jc,ic)) &
                    )*udx3mr(kc)*0.5d0
        !
        !
        !    sal vz term
        !
        !
        !                d  sal q_z
        !             -----------
        !                d   z      
        !
                hsz=(vzr(kc,jc,ip)*(sal(kc,jc,ip)+sal(kc,jc,ic))- &
                    vzr(kc,jc,ic)*(sal(kc,jc,ic)+sal(kc,jc,im)) &
                    )*udzr
        !
        !
        !    sal vyr term
        !
        !
        !                d  sal q_y 
        !             -----------
        !                d   y      
        !
                hsy=(vyr(kc,jp,ic)*(sal(kc,jp,ic)+sal(kc,jc,ic))- &
                    vyr(kc,jc,ic)*(sal(kc,jc,ic)+sal(kc,jm,ic)) &
                    )*udyr
        !
        !
        !   zz second derivatives of sal
        !
                ! dzzs=(sal(kc,jc,ip) &
                ! -2.0*sal(kc,jc,ic) &
                !     +sal(kc,jc,im))*udzrq
                if (solidr(kc,jc,ip)) then
                    dzzs = (sal(kc,jc,im) - sal(kc,jc,ic))*udzrq
                elseif (solidr(kc,jc,im)) then
                    dzzs = (sal(kc,jc,ip) - sal(kc,jc,ic))*udzrq
                else
                    dzzs = (sal(kc,jc,ip) - 2.0*sal(kc,jc,ic) + sal(kc,jc,im))*udzrq
                end if
        !
        !   yy second derivatives of sal
        !
                ! dyys=(sal(kc,jp,ic) &
                ! -2.0*sal(kc,jc,ic) &
                !     +sal(kc,jm,ic))*udyrq
                if (solidr(kc,jp,ic)) then
                    dyys = (sal(kc,jm,ic) - sal(kc,jc,ic))*udyrq
                elseif (solidr(kc,jm,ic)) then
                    dyys = (sal(kc,jp,ic) - sal(kc,jc,ic))*udyrq
                else
                    dyys = (sal(kc,jp,ic) - 2.0*sal(kc,jc,ic) + sal(kc,jm,ic))*udyrq
                end if

                hsal(kc,jc,ic) = -(hsx+hsy+hsz)+dyys+dzzs
            enddo
            
            kc = nxmr
            kp = nxr
            km = nxmr - 1
            hsx=(vxr(kp,jc,ic)*2.d0*saltp(1,jc,ic)- &
                vxr(kc,jc,ic)*(sal(kc,jc,ic)+sal(km,jc,ic)) &
                )*udx3mr(kc)*0.5d0
            hsz=(vzr(kc,jc,ip)*(sal(kc,jc,ip)+sal(kc,jc,ic))- &
                vzr(kc,jc,ic)*(sal(kc,jc,ic)+sal(kc,jc,im)) &
                )*udzr
            hsy=(vyr(kc,jp,ic)*(sal(kc,jp,ic)+sal(kc,jc,ic))- &
                vyr(kc,jc,ic)*(sal(kc,jc,ic)+sal(kc,jm,ic)) &
                )*udyr
            if (solidr(kc,jc,ip)) then
                dzzs = (sal(kc,jc,im) - sal(kc,jc,ic))*udzrq
            elseif (solidr(kc,jc,im)) then
                dzzs = (sal(kc,jc,ip) - sal(kc,jc,ic))*udzrq
            else
                dzzs = (sal(kc,jc,ip) - 2.0*sal(kc,jc,ic) + sal(kc,jc,im))*udzrq
            end if
            if (solidr(kc,jp,ic)) then
                dyys = (sal(kc,jm,ic) - sal(kc,jc,ic))*udyrq
            elseif (solidr(kc,jm,ic)) then
                dyys = (sal(kc,jp,ic) - sal(kc,jc,ic))*udyrq
            else
                dyys = (sal(kc,jp,ic) - 2.0*sal(kc,jc,ic) + sal(kc,jm,ic))*udyrq
            end if
            ! dzzs=(sal(kc,jc,ip) - 2.0*sal(kc,jc,ic) + sal(kc,jc,im))*udzrq
            ! dyys=(sal(kc,jp,ic) - 2.0*sal(kc,jc,ic) + sal(kc,jm,ic))*udyrq
            hsal(kc,jc,ic) = -(hsx+hsy+hsz)+dyys+dzzs

        enddo
    enddo
!$OMP  END PARALLEL DO

    if (phasefield) then
        do ic=xstartr(3),xendr(3)
            im = ic - 1
            ip = ic + 1
            do jc=xstartr(2),xendr(2)
                jm = jc - 1
                jp = jc + 1

                kc = 1
                kp = 2
                dpdsx = sdx(kc)**2*(phi(kp,jc,ic) - phi(kc,jc,ic))&
                        *(sal(kp,jc,ic) - sal(kc,jc,ic) &
                            + 2.0*SfixS*(sal(kc,jc,ic) - salbp(1,jc,ic)))
                dpdsy = udyr**2*(phi(kc,jp,ic) - phi(kc,jm,ic))&
                        *(sal(kc,jp,ic) - sal(kc,jm,ic))
                dpdsz = udzr**2*(phi(kc,jc,ip) - phi(kc,jc,im))&
                        *(sal(kc,jc,ip) - sal(kc,jc,im))
                phfrac = 1.0/(1.0 - phi(kc,jc,ic) + pf_delta)
                hsal(kc,jc,ic) = hsal(kc,jc,ic) &
                        + phfrac*(sal(kc,jc,ic)*rhsr(kc,jc,ic)*aldt &
                            - (dpdsx + dpdsy + dpdsz)/pecs)

                do kc=2,nxmr-1
                    km = kc - 1
                    kp = kc + 1
                    dpdsx = sdx(kc)**2*(phi(kp,jc,ic) - phi(km,jc,ic))&
                            *(sal(kp,jc,ic) - sal(km,jc,ic))
                    dpdsy = udyr**2*(phi(kc,jp,ic) - phi(kc,jm,ic))&
                            *(sal(kc,jp,ic) - sal(kc,jm,ic))
                    dpdsz = udzr**2*(phi(kc,jc,ip) - phi(kc,jc,im))&
                            *(sal(kc,jc,ip) - sal(kc,jc,im))
                    phfrac = 1.0/(1.0 - phi(kc,jc,ic) + pf_delta)
                    hsal(kc,jc,ic) = hsal(kc,jc,ic) &
                        + phfrac*(sal(kc,jc,ic)*rhsr(kc,jc,ic)*aldt &
                            - (dpdsx + dpdsy + dpdsz)/pecs)
                end do

                kc = nxmr
                km = nxmr - 1
                dpdsx = sdx(kc)**2*(phi(kc,jc,ic) - phi(km,jc,ic))&
                        *(sal(kc,jc,ic) - sal(km,jc,ic) &
                            + 2.0*SfixN*(saltp(1,jc,ic) - sal(kc,jc,ic)))
                dpdsy = udyr**2*(phi(kc,jp,ic) - phi(kc,jm,ic))&
                        *(sal(kc,jp,ic) - sal(kc,jm,ic))
                dpdsz = udzr**2*(phi(kc,jc,ip) - phi(kc,jc,im))&
                        *(sal(kc,jc,ip) - sal(kc,jc,im))
                phfrac = 1.0/(1.0 - phi(kc,jc,ic) + pf_delta)
                hsal(kc,jc,ic) = hsal(kc,jc,ic) &
                    + phfrac*(sal(kc,jc,ic)*rhsr(kc,jc,ic)*aldt &
                        - (dpdsx + dpdsy + dpdsz)/pecs)
            end do
        end do
    end if

    return

end subroutine ExplicitTermsSal