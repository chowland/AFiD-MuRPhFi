!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: ExplicitTermsVX.F90                            !
!    CONTAINS: subroutine ExplicitTermsVX                 !
!                                                         ! 
!    PURPOSE: Compute the non-linear terms associated to  !
!     the velocity in the x (vertical) dimension          !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ExplicitTermsVX
    use param
    use local_arrays, only: vz,vy,vx,temp,qcap
    use mgrd_arrays, only: salc,phic
    use decomp_2d, only: xstart,xend
    implicit none
    integer :: jc,kc
    integer :: km,kp,jm,jp,ic,im,ip
    real    :: hxx,hxy,hxz
    real    :: udz,udy,tempit,salit,volpen
    real    :: udzq,udyq
    real    :: dzzvx,dyyvx,pf_eta
    real, dimension(nxm) :: idx
    
    pf_eta = ren*(1.51044385*pf_eps)**2
    
    udy=dy*0.25
    udz=dz*0.25
    
    udyq=dyq/ren
    udzq=dzq/ren

    do kc=1,nxm
        idx(kc) = 1.0/udx3m(kc)
    end do
    
    !$OMP  PARALLEL DO &
    !$OMP   DEFAULT(none) &
    !$OMP   SHARED(xstart,xend,nxm,vz,vy,vx,dz,dy) &
    !$OMP   SHARED(kmv,kpv,am3sk,ac3sk,ap3sk,udz) &
    !$OMP   SHARED(udy,udzq,udyq,udx3c,qcap,temp) &
    !$OMP   PRIVATE(ic,jc,kc,imm,ipp,km,kp) &
    !$OMP   PRIVATE(jmm,jpp,tempit) &
    !$OMP   PRIVATE(hxx,hxy,hxz,dzzvx,dyyvx)
    do ic=xstart(3),xend(3)
        im=ic-1
        ip=ic+1
        do jc=xstart(2),xend(2)
            jm=jc-1
            jp=jc+1
            do kc=2,nxm
                km=kc-1
                kp=kc+1
                !
                !    vx vz term
                !
                !
                !                d  q_x q_t 
                !             -----------
                !                d   t      
                !
                !
                ! hxz=(((vz(kc,jc,ipp)+vz(km,jc,ipp)) &
                ! *(vx(kc,jc,ipp)+vx(kc,jc,ic))) &
                ! -((vz(kc,jc,ic)+vz(km,jc,ic)) &
                ! *(vx(kc,jc,ic)+vx(kc,jc,imm))))*udz
                hxz = udx3c(kc)*udz*( &
                    (idx(km)*vz(km,jc,ip) + idx(kc)*vz(kc,jc,ip))*(vx(kc,jc,ip) + vx(kc,jc,ic)) &
                    - (idx(km)*vz(km,jc,ic) + idx(kc)*vz(kc,jc,ic))*(vx(kc,jc,ic) + vx(kc,jc,im)) &
                )
                !
                !    vx vy term
                !
                !
                !                d  q_x q_r 
                !             -----------
                !                d   r      
                !
                ! hxy=(((vy(kc,jpp,ic)+vy(km,jpp,ic)) &
                ! *(vx(kc,jpp,ic)+vx(kc,jc,ic))) &
                ! -((vy(kc,jc,ic)+vy(km,jc,ic)) &
                ! *(vx(kc,jc,ic)+vx(kc,jmm,ic))))*udy
                hxy = udx3c(kc)*udy*( &
                    (idx(km)*vy(km,jp,ic) + idx(kc)*vy(kc,jp,ic))*(vx(kc,jp,ic) + vx(kc,jc,ic)) &
                    - (idx(km)*vy(km,jc,ic) + idx(kc)*vy(kc,jc,ic))*(vx(kc,jc,ic) + vx(kc,jm,ic)) &
                )
                !
                !    vx vx term
                !
                !
                !                 d  q_x q_x 
                !                -----------
                !                 d   x      
                !
                hxx=((vx(kp,jc,ic)+vx(kc,jc,ic))*(vx(kp,jc,ic)+vx(kc,jc,ic)) &
                -(vx(kc,jc,ic)+vx(km,jc,ic))*(vx(kc,jc,ic)+vx(km,jc,ic)) &
                )*udx3c(kc)*0.25d0
                
                !
                !   z second derivatives of vx
                !
                dzzvx=(vx(kc,jc,im) &
                -2.0*vx(kc,jc,ic) &
                +vx(kc,jc,ip))*udzq
                !
                !   y second derivatives of vx
                !
                dyyvx=(vx(kc,jm,ic) &
                -2.0*vx(kc,jc,ic) &
                +vx(kc,jp,ic))*udyq
                
                
                qcap(kc,jc,ic) =-(hxx+hxy+hxz)+dyyvx+dzzvx
                
            end do
        end do
    end do
    
    !CJH Add the buoyancy term if x is chosen gAxis
    if (gAxis.eq.1) then
        do ic=xstart(3),xend(3)
            do jc=xstart(2),xend(2)
                do kc=2,nxm
                    km=kc-1
                    tempit=active_T*0.5d0*(temp(kc,jc,ic)+temp(km,jc,ic))
                    qcap(kc,jc,ic) = qcap(kc,jc,ic) + byct*tempit
                end do
            end do
        end do
        
        !CJH Add salinity component of buoyancy if used
        if (salinity) then
            do ic=xstart(3),xend(3)
                do jc=xstart(2),xend(2)
                    do kc=2,nxm
                        km=kc-1
                        salit =active_S*0.5d0*(salc(kc,jc,ic)+salc(km,jc,ic))
                        qcap(kc,jc,ic) = qcap(kc,jc,ic) - bycs*salit
                    end do
                end do
            end do

            if (melt) then
                !CJH For melting boundary, remove far-field buoyancy
                do ic=xstart(3),xend(3)
                    do jc=xstart(2),xend(2)
                        do kc=2,nxm
                            qcap(kc,jc,ic) = qcap(kc,jc,ic) + bycs - byct
                        end do
                    end do
                end do
            end if
        end if
    end if
    
    if (phasefield .and. .not.IBM) then
        do ic=xstart(3),xend(3)
            do jc=xstart(2),xend(2)
                do kc=2,nxm
                    km = kc - 1
                    volpen = 0.5d0*(phic(kc,jc,ic) + phic(km,jc,ic))* &
                    vx(kc,jc,ic)/pf_eta
                    ! volpen = 0.5*(1.0+tanh((xc(kc)-0.5)/2/pf_eps))* &
                    !           vx(kc,jc,ic)/pf_gamma
                    qcap(kc,jc,ic) = qcap(kc,jc,ic) - volpen
                end do
            end do
        end do
    end if
    !$OMP END PARALLEL DO
    
    return
end
!
