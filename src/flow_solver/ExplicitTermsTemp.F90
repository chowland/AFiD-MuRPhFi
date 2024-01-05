!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: ExplicitTermsTemp.F90                          !
!    CONTAINS: subroutine ExplicitTermsTemp               !
!                                                         ! 
!    PURPOSE: Compute the non-linear terms associated to  !
!     the temperature.                                    !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ExplicitTermsTemp
    use param
    use local_arrays, only: vy,vx,temp,vz,hro
    use decomp_2d, only: xstart,xend
    implicit none
    integer :: jc,kc,ic
    integer :: jm,jp,im,ip
    real    :: htx,hty,htz,udy,udz
    real    :: udzq,udyq
    real    :: dyyt,dzzt
    
    udz=dz*0.5d0
    udy=dy*0.5d0
    udzq=dzq/pect
    udyq=dyq/pect
    
    !$OMP  PARALLEL DO &
    !$OMP   DEFAULT(none) &
    !$OMP   SHARED(xstart,xend,vz,vy,vx,nxm) &
    !$OMP   SHARED(kmv,kpv,am3sk,ac3sk,ap3sk,udz) &
    !$OMP   SHARED(udy,udzq,udyq,udx3c,temp,hro) &
    !$OMP   PRIVATE(ic,jc,kc,im,ip,km,kp,jm,jp) &
    !$OMP   PRIVATE(htx,hty,htz,dyyt,dzzt)
    do ic=xstart(3),xend(3)
        im=ic-1
        ip=ic+1
        do jc=xstart(2),xend(2)
            jm=jc-1
            jp=jc+1
            do kc=1,nxm

                ! ! x-advection d/dx(vx * T)
                if (kc==1) then
                    htx = ( &
                          vx(kc+1,jc,ic)*(temp(kc+1,jc,ic) + temp(kc,jc,ic)) &
                        - vx(kc  ,jc,ic)*2.0*tempbp(1,jc,ic) &
                    )*0.5*udx3m(kc)
                elseif (kc==nxm) then
                    htx = ( &
                          vx(kc+1,jc,ic)*2.0*temptp(1,jc,ic) &
                        - vx(kc  ,jc,ic)*(temp(kc,jc,ic) + temp(kc-1,jc,ic)) &
                    )*0.5*udx3m(kc)
                else
                    htx = ( &
                          vx(kc+1,jc,ic)*(temp(kc+1,jc,ic) + temp(kc  ,jc,ic)) &
                        - vx(kc  ,jc,ic)*(temp(kc  ,jc,ic) + temp(kc-1,jc,ic)) &
                    )*0.5*udx3m(kc)
                end if

                ! ! z-advection d/dx(vz * T)
                htz = ( &
                      vz(kc,jc,ip)*(temp(kc,jc,ip)+temp(kc,jc,ic)) &
                    - vz(kc,jc,ic)*(temp(kc,jc,ic)+temp(kc,jc,im)) &
                )*udz
                
                ! ! y-advection d/dx(vy * T)
                hty=( &
                      vy(kc,jp,ic)*(temp(kc,jp,ic)+temp(kc,jc,ic)) &
                    - vy(kc,jc,ic)*(temp(kc,jc,ic)+temp(kc,jm,ic)) &
                )*udy

                !   zz second derivatives of temp
                dzzt=(temp(kc,jc,ip) - 2.0*temp(kc,jc,ic) + temp(kc,jc,im))*udzq
                !   yy second derivatives of temp
                dyyt=(temp(kc,jp,ic) - 2.0*temp(kc,jc,ic) + temp(kc,jm,ic))*udyq

                hro(kc,jc,ic) = -(htx+hty+htz)+dyyt+dzzt
                !   hro(kc,jc,ic) = dyyt+dzzt
            end do
            
        end do
    end do
    !$OMP  END PARALLEL DO
    
    return
end