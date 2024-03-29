!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: ExplicitTermsVY.F90                            !
!    CONTAINS: subroutine ExplicitTermsVY                 !
!                                                         ! 
!    PURPOSE: Compute the non-linear terms associated to  !
!     the velocity in the y (horizontal) dimension        !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ExplicitTermsVY
    use param
    use local_arrays, only: vx,vy,vz,dph,temp
    use decomp_2d, only: xstart,xend
    implicit none
    integer :: kc,kp,jpp,jmm,jc,ic,imm,ipp
    integer :: kpp,kmm
    real    :: udzq,udyq
    real    :: udy,udz,hyx,hyy,hyz 
    real    :: dyyvy, dzzvy
    real    :: tempit, Gy

    udyq=dyq/ren
    udzq=dzq/ren

    udy=dy*0.25
    udz=dz*0.25

    !CJH Treat dPdy input variable as Re_tau
    if (inslwn==0 .or. inslws==0) then
        Gy = dPdy**2/ren**2
    else
        Gy = 8.0*dPdy**2/ren**2
    end if

!$OMP  PARALLEL DO &
!$OMP  DEFAULT(none) &
!$OMP  SHARED(xstart,xend,vz,vy,vx,dz,dy) &
!$OMP  SHARED(kmv,kpv,am3sk,ac3sk,ap3sk,udz) &
!$OMP  SHARED(udy,udzq,udyq,udx3m,dph,nxm) &
!$OMP  PRIVATE(ic,jc,kc,imm,ipp,kmm,kp,kpp) &
!$OMP  PRIVATE(jmm,jpp) &
!$OMP  PRIVATE(hyx,hyy,hyz,dyyvy,dzzvy)
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
!
!     vy vx term
!
!                 d  q_x q_r 
!                -----------
!                 d   x      
!
                hyx=((vx(kp,jc,ic)+vx(kp,jmm,ic))*(vy(kpp,jc,ic)+vy(kc,jc,ic)) &
                    -(vx(kc,jc,ic)+vx(kc,jmm,ic))*(vy(kc,jc,ic)+vy(kmm,jc,ic)) &
                    )*udx3m(kc)*0.25d0

!     
!     vy vy term
!
!                 d  q_r q_r 
!                ------------
!                 d   r      
!
                hyy=( (vy(kc,jpp,ic)+vy(kc,jc,ic)) &
                    *(vy(kc,jpp,ic)+vy(kc,jc,ic)) &
                    -(vy(kc,jmm,ic)+vy(kc,jc,ic)) &
                    *(vy(kc,jmm,ic)+vy(kc,jc,ic)) &
                    )*udy
!
!     vz vy term
!
!                 d  q_t q_r 
!                ------------
!                 d   t      
!
                hyz=( (vy(kc,jc,ipp)+vy(kc,jc,ic)) &
                    *(vz(kc,jc,ipp)+vz(kc,jmm,ipp)) &
                    -(vy(kc,jc,ic)+vy(kc,jc,imm)) &
                    *(vz(kc,jc,ic)+vz(kc,jmm,ic)) &
                    )*udz

!
!   yy second derivative of vy
!
                dyyvy=(vy(kc,jpp,ic) &
                    -2.0*vy(kc,jc,ic) &
                    +vy(kc,jmm,ic))*udyq

!
!   zz second derivative of vy
!
                dzzvy=(vy(kc,jc,ipp) &
                    -2.0*vy(kc,jc,ic) &
                    +vy(kc,jc,imm))*udzq


                dph(kc,jc,ic)=-(hyx+hyy+hyz)+dyyvy+dzzvy + Gy

            enddo
        enddo
    enddo

    !CJH Add the buoyancy term if y is chosen gAxis
    if (gAxis.eq.2) then
        do ic=xstart(3),xend(3)
            do jc=xstart(2),xend(2)
                jmm=jc-1
                do kc=1,nxm
                    tempit=active_T*0.5d0*(temp(kc,jc,ic)+temp(kc,jmm,ic))
                    dph(kc,jc,ic) = dph(kc,jc,ic) + byct*tempit
                end do
            end do
        end do
    end if

!$OMP END PARALLEL DO

    return
end subroutine ExplicitTermsVY