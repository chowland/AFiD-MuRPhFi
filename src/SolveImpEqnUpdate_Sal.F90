!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: SolveImpEqnUpdate_Sal.F90                      !
!    CONTAINS: subroutine SolveImpEqnUpdate_Sal           !
!                                                         ! 
!    PURPOSE: Inverts the implicit equation for           !
!     salinity, and updates it to time t+dt               !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine SolveImpEqnUpdate_Sal
      use param
      use local_arrays, only: sal,rhsr
      use decomp_2d, only: xstartr,xendr
      implicit none
      real, dimension(nxr) :: amkl,apkl,ackl
      integer :: jc,kc,info,ipkv(nxr),ic,nrhs
      real :: betadx,ackl_b
      real :: amkT(nxr-1),ackT(nxr),apkT(nxr-1),appk(nxr-2)

!     Calculate the coefficients of the tridiagonal matrix
!     The coefficients are normalized to prevent floating
!     point errors.

      betadx=0.5d0*al*dt/pec !TODO

      amkl(1)=0.d0
      apkl(1)=0.d0
      ackl(1)=1.d0
      do kc=2,nxmr
        ackl_b=1.0d0/(1.-ac3sskr(kc)*betadx)
        amkl(kc)=-am3sskr(kc)*betadx*ackl_b
        ackl(kc)=1.0d0
        apkl(kc)=-ap3sskr(kc)*betadx*ackl_b
      enddo
      amkl(nxr)=0.d0
      apkl(nxr)=0.d0
      ackl(nxr)=1.d0

      amkT=amkl(2:nxr)
      apkT=apkl(1:nxmr)
      ackT=ackl(1:nxr)

!     Call to LAPACK library to factor tridiagonal matrix.
!     No solving is done in this call.

      call dgttrf(nxr,amkT,ackT,apkT,appk,ipkv,info)
     
      nrhs=(xendr(3)-xstartr(3)+1)*(xendr(2)-xstartr(2)+1)
      do ic=xstartr(3),xendr(3)
        do jc=xstartr(2),xendr(2)
           do kc=2,nxmr
              ackl_b=1.0/(1.0-ac3sskr(kc)*betadx)
              rhsr(kc,jc,ic)=rhsr(kc,jc,ic)*ackl_b
           end do
        end do
      end do
      
      call dgttrs('N',nxr,nrhs,amkT,ackT,apkT,appk,ipkv,rhsr,nxr,info)

       do ic=xstartr(3),xendr(3)
         do jc=xstartr(2),xendr(2)
            do kc=2,nxmr
              sal(kc,jc,ic)=sal(kc,jc,ic) + rhsr(kc,jc,ic)
             end do
          end do
      end do

      return
      end
