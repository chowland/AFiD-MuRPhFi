!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: ExplicitTermsVZ.F90                            !
!    CONTAINS: subroutine ExplicitTermsVZ                 !
!                                                         ! 
!    PURPOSE: Compute the non-linear terms associated to  !
!     the velocity in the z (horizontal) dimension        !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ExplicitTermsVZ
    use param
    use local_arrays, only: vx,vy,vz,dq,temp
    use mgrd_arrays, only: salc,phic
    use decomp_2d, only: xstart,xend
    implicit none
    integer :: kc,kp,jpp,jmm,jc,ic,imm,ipp
    integer :: kmm,kpp
    real    :: hzx,hzy,hzz,udy,udz
    real    :: udyq,udzq
    real    :: dzzvz,dyyvz,pf_eta
    real    :: tempit, salit,volpen

    pf_eta = ren*(1.51044385*pf_eps)**2

    udyq=dyq/ren
    udzq=dzq/ren

    udy=dy*0.25
    udz=dz*0.25

!$OMP  PARALLEL DO &
!$OMP  DEFAULT(none) &
!$OMP  SHARED(xstart,xend,vz,vy,vx,dz,dy,udx3m) &
!$OMP  SHARED(kmv,kpv,am3sk,ac3sk,ap3sk,udz) &
!$OMP  SHARED(udy,udzq,udyq,dq,nxm) &
!$OMP  PRIVATE(ic,jc,kc,imm,ipp,kmm,kp,kpp) &
!$OMP  PRIVATE(jmm,jpp) &
!$OMP  PRIVATE(hzz,hzy,hzx,dzzvz,dyyvz)

    do ic=xstart(3),xend(3)
        imm=ic-1
        ipp=ic+1
        do jc=xstart(2),xend(2)
            jmm=jc-1
            jpp=jc+1
            do kc=1,nxm
                kmm=kmv(kc)
                kpp=kpv(kc)
                kp=kc+1
        
!     vz vz term
!
!
!                 d  q_t q_t 
!                ------------
!                 d   t      
!
                hzz=( (vz(kc,jc,ipp)+vz(kc,jc,ic)) &
                    *(vz(kc,jc,ipp)+vz(kc,jc,ic)) &
                    -(vz(kc,jc,imm)+vz(kc,jc,ic)) &
                    *(vz(kc,jc,imm)+vz(kc,jc,ic)) &
                    )*udz

!     vz vy term
!
!
!                 d  q_t q_r 
!                ------------
!                 d   r      
!
                hzy=( (vy(kc,jpp,ic)+vy(kc,jpp,imm)) &
                    *(vz(kc,jpp,ic)+vz(kc,jc,ic)) &
                    -(vy(kc,jc,ic)+vy(kc,jc,imm)) &
                    *(vz(kc,jc,ic)+vz(kc,jmm,ic)) &
                    )*udy
!
!     vz vx term
!
!
!                 d  q_t q_x 
!                -----------
!                 d   x      
!
                hzx=((vx(kp,jc,ic)+vx(kp,jc,imm))*(vz(kpp,jc,ic)+vz(kc,jc,ic)) &
                    -(vx(kc,jc,ic)+vx(kc,jc,imm))*(vz(kc,jc,ic)+vz(kmm,jc,ic)) &
                    )*udx3m(kc)*0.25d0
!
!
!
!   11 second derivative of vz
!
                dzzvz=(vz(kc,jc,ipp) &
                    -2.0*vz(kc,jc,ic) &
                    +vz(kc,jc,imm))*udzq
!
!   22 second derivative of vz
!
                dyyvz=(vz(kc,jpp,ic) &
                    -2.0*vz(kc,jc,ic) &
                    +vz(kc,jmm,ic))*udyq

!
                dq(kc,jc,ic)=-(hzx+hzy+hzz)+dyyvz+dzzvz - dPdz
!
            enddo
        enddo
    enddo

    !CJH Add the buoyancy term if z is chosen gAxis
    if (gAxis.eq.3) then
        do ic=xstart(3),xend(3)
            imm=ic-1
            do jc=xstart(2),xend(2)
                do kc=1,nxm
                    tempit=active_T*0.5d0*(temp(kc,jc,ic)+temp(kc,jc,imm))
                    dq(kc,jc,ic) = dq(kc,jc,ic) + byct*tempit
                end do
            end do
        end do

        !CJH Add salinity component of buoyancy if used
        if (salinity) then
            do ic=xstart(3),xend(3)
                imm=ic-1
                do jc=xstart(2),xend(2)
                    do kc=1,nxm
                        salit =active_S*0.5d0*(salc(kc,jc,ic)+salc(kc,jc,imm))
                        dq(kc,jc,ic) = dq(kc,jc,ic) - bycs*salit
                    end do
                end do
            end do
        end if
    end if
!$OMP END PARALLEL DO

    if (phasefield .and. .not.IBM) then
        do ic=xstart(3),xend(3)
            imm=ic-1
            do jc=xstart(2),xend(2)
                do kc=1,nxm
                volpen = 0.5d0*(phic(kc,jc,ic) + phic(kc,jc,imm))* &
                            vz(kc,jc,ic)/pf_eta
                dq(kc,jc,ic) = dq(kc,jc,ic) - volpen
                end do
            end do
        end do
    end if

    return
end subroutine ExplicitTermsVZ