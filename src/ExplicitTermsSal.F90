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
      use mgrd_arrays, only: vxr,vyr,vzr
      use local_arrays, only: sal,hsal
      use decomp_2d, only: xstartr,xendr
      implicit none
      integer :: jc,kc,ic
      integer :: km,kp,jm,jp,im,ip
      real    :: hsx,hsy,hsz,udyr,udzr
      real    :: udzrq,udyrq
      real    :: dyys,dzzs

      udzr=dzr*0.25
      udyr=dyr*0.25
      udzrq=dzq/pecs
      udyrq=dyq/pecs

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
        do kc=2,nxmr
         km=kc-1
         kp=kc+1
!
!
!    sal vz term
!
!
!                d  sal q_z
!             -----------
!                d   z      
!
      hsz=((vzr(km,jc,ip)+vzr(kc,jc,ip))*(sal(kc,jc,ip)+sal(kc,jc,ic))- &
           (vzr(km,jc,ic)+vzr(kc,jc,ic))*(sal(kc,jc,ic)+sal(kc,jc,im)) &
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
      hsy=((vyr(kc,jp,ic)+vyr(km,jp,ic))*(sal(kc,jp,ic)+sal(kc,jc,ic))- &
           (vyr(kc,jc,ic)+vyr(km,jc,ic))*(sal(kc,jc,ic)+sal(kc,jm,ic)) &
          )*udyr
!
!    sal vxr term
!
!
!                 d  sal q_x 
!                -----------
!                 d   x      
!
      hsx=((vxr(kp,jc,ic)+vxr(kc,jc,ic))*(sal(kp,jc,ic)+sal(kc,jc,ic))- &
           (vxr(kc,jc,ic)+vxr(km,jc,ic))*(sal(kc,jc,ic)+sal(km,jc,ic)) &
          )*udx3cr(kc)*0.25d0
!
!
!   zz second derivatives of sal
!
            dzzs=(sal(kc,jc,ip) &
             -2.0*sal(kc,jc,ic) &
                 +sal(kc,jc,im))*udzrq
      
!
!   yy second derivatives of sal
!
            dyys=(sal(kc,jp,ic) &
             -2.0*sal(kc,jc,ic) &
                 +sal(kc,jm,ic))*udyrq
!
            hsal(kc,jc,ic) = -(hsx+hsy+hsz)+dyys+dzzs
      enddo
      enddo
      enddo
!$OMP  END PARALLEL DO

      return
      end
