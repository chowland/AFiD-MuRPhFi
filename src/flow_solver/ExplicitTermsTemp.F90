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
end  subroutine ExplicitTermsTemp

subroutine ExplicitTermsTempPhi
    use param
    use local_arrays, only: vy,vx,temp,vz,hro
    use afid_phasefield, only:  phic
    use decomp_2d, only: xstart,xend
    implicit none
    integer :: jc,kc,ic
    integer :: km,kp,jm,jp,im,ip
    real    :: htx,hty,htz,udy,udz
    real    :: udzq,udyq
    real    :: dyyt,dzzt
    real    :: i_midp,i_midm,j_midp,j_midm,k_midp,k_midm,tmp1,tmp2
    real    :: rho,Cp
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
			    km=kc-1
				kp=kc+1
                ! ! x-advection d/dx(vx * T)
                ! ! x-advection d/dx(vx * T)
				tmp1=vx(kc+1,jc,ic)
				tmp2=vx(kc  ,jc,ic)
				if(tmp1.ge.0)then
					k_midp=kc
				else
					k_midp=kp
				endif
				if(tmp2.ge.0)then
					k_midm=km
				else
					k_midm=kc
				endif 
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
                    ! htx = ( &
                    !         vx(kc+1,jc,ic)*((temp(kc+1,jc,ic) + temp(kc  ,jc,ic))*0.5 &
                    !           -E_quick/3.0*(temp(k_midp-1,jc,ic)-2.0*temp(k_midp,jc,ic)+temp(k_midp+1,jc,ic))) &
                    !       - vx(kc  ,jc,ic)*((temp(kc  ,jc,ic) + temp(kc-1,jc,ic))*0.5 &
                    !           -E_quick/3.0*(temp(k_midm-1,jc,ic)-2.0*temp(k_midm,jc,ic)+temp(k_midm+1,jc,ic))) &
                    !       )*dx
                end if

                ! ! z-advection d/dx(vz * T)
                tmp1=vz(kc,jc,ip)
				tmp2=vz(kc,jc,ic)
				if(tmp1.ge.0)then
					i_midp=ic
				else
					i_midp=ip
				endif
				if(tmp2.ge.0)then
					i_midm=im
				else
					i_midm=ic
				endif 
                htz = ( &
                      vz(kc,jc,ip)*(temp(kc,jc,ip)+temp(kc,jc,ic)) &
                    - vz(kc,jc,ic)*(temp(kc,jc,ic)+temp(kc,jc,im)) &
                )*udz
                ! htz = ( &
                !          vz(kc,jc,ip)*((temp(kc,jc,ip)+temp(kc,jc,ic))*0.5 &
                !               -E_quick/3.0*(temp(kc,jc,i_midp-1)-2.0*temp(kc,jc,i_midp)+temp(kc,jc,i_midp+1))) &
                !        - vz(kc,jc,ic)*((temp(kc,jc,ic)+temp(kc,jc,im))*0.5 &
                !               -E_quick/3.0*(temp(kc,jc,i_midm-1)-2.0*temp(kc,jc,i_midm)+temp(kc,jc,i_midm+1))) &
                !  )*dz
                ! ! y-advection d/dx(vy * T)
                 tmp1=vy(kc,jp,ic)
                 tmp2=vy(kc,jc,ic)
                 if(tmp1.ge.0)then
                     j_midp=jc
                 else
                     j_midp=jp
                 endif
                 if(tmp2.ge.0)then
                     j_midm=jm
                 else
                     j_midm=jc
                 endif
                hty=( &
                      vy(kc,jp,ic)*(temp(kc,jp,ic)+temp(kc,jc,ic)) &
                    - vy(kc,jc,ic)*(temp(kc,jc,ic)+temp(kc,jm,ic)) &
                )*udy
                 
                !  hty=( &
                !          vy(kc,jp,ic)*((temp(kc,jp,ic)+temp(kc,jc,ic))*0.5&
                !                -E_quick/3.0*(temp(kc,j_midp-1,ic)-2.0*temp(kc,j_midp,ic)+temp(kc,j_midp+1,ic))) &
                !         - vy(kc,jc,ic)*((temp(kc,jc,ic)+temp(kc,jm,ic))*0.5&
                !                -E_quick/3.0*(temp(kc,j_midm-1,ic)-2.0*temp(kc,j_midm,ic)+temp(kc,j_midm+1,ic))) &
                !     )*dy
                !   zz second derivatives of temp
                rho = (1.0-phic(kc,jc,ic))*(1.0-rd_sl) +rd_sl
                Cp  = (1.0-phic(kc,jc,ic))*(1.0-rCp_sl)+rCp_sl
                !dzzt=(temp(kc,jc,ip) - 2.0*temp(kc,jc,ic) + temp(kc,jc,im))*udzq
                dzzt= (&
                        (temp(kc,jc,ip) - temp(kc,jc,ic))*((2.0-(phic(kc,jc,ip)+phic(kc,jc,ic)))*0.5*(1.0-k_c)+k_c)&
                       -(temp(kc,jc,ic) - temp(kc,jc,im))*((2.0-(phic(kc,jc,im)+phic(kc,jc,ic)))*0.5*(1.0-k_c)+k_c)&
                        )*udzq
                !   yy second derivatives of temp
                !dyyt=(temp(kc,jp,ic) - 2.0*temp(kc,jc,ic) + temp(kc,jm,ic))*udyq
                dyyt= (&
                        (temp(kc,jp,ic) - temp(kc,jc,ic))*((2.0-(phic(kc,jp,ic)+phic(kc,jc,ic)))*0.5*(1.0-k_c)+k_c)&
                       -(temp(kc,jc,ic) - temp(kc,jm,ic))*((2.0-(phic(kc,jm,ic)+phic(kc,jc,ic)))*0.5*(1.0-k_c)+k_c)&
                        )*udyq
                hro(kc,jc,ic) = -(htx+hty+htz)+(dyyt+dzzt)/rho/Cp
                !   hro(kc,jc,ic) = dyyt+dzzt
            end do
            
        end do
    end do
    !$OMP  END PARALLEL DO
    
    return
end  subroutine ExplicitTermsTempPhi