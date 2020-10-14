!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: CreateInitialConditions.F90                    !
!    CONTAINS: subroutine CreateInitialConditions         !
!                                                         ! 
!    PURPOSE: Initialization routine. Sets initial        !
!     conditions for velocity and temperature             !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine CreateInitialConditions
      use param
      use local_arrays, only: vy,vx,temp,vz,sal
      use decomp_2d, only: xstart,xend,xstartr,xendr
      use mpih
      implicit none
      integer :: j,k,i
      real :: xxx,yyy,eps,varptb

      call random_seed()
      eps=0.5d0
      do i=xstart(3),xend(3)
      do j=xstart(2),xend(2)
      do k=1,nxm
           vz(k,j,i)=0.0d0
           yyy=xm(k) 
           xxx=yc(j)            
           call random_number(varptb)
           varptb=0.d0
           vy(k,j,i)=-yyy*(1.d0-yyy)*(2.d0*yyy-1.d0)*sin(pi*xxx)
          !  vy(k,j,i)=(xm(k))+0.3*(2.d0*varptb-1.d0)+ &
!           vy(k,j,i)=+0.3*(2.d0*varptb-1.d0)+ &
!     &                  (2.0d0*yyy-6.0d0*yyy**2+4.0d0*yyy**3) &
!     &                  *sin(3*xxx)*eps

           yyy=xc(k)          
           xxx=ym(j)
           vx(k,j,i)=-yyy**2*(1.0d0-yyy)**2*cos(pi*xxx)*eps
!           vx(k,j,i)=1.d0 ! -yyy**2*(1.0d0-yyy)**2*cos(3.1*xxx)*eps

         enddo
        enddo
      enddo

      !assign linear temperature profile in the nodes k=1 to k=nxm
      do i=xstart(3),xend(3)
      do j=xstart(2),xend(2)
      do k=1,nxm
      temp(k,j,i)=tempbp(j,i)+(temptp(j,i)-tempbp(j,i))*xm(k)/xc(nx)
      enddo
      enddo
      enddo

      !assign the boundary conditions at k=1 and k=nx
      ! do i=xstart(3),xend(3)
      ! do j=xstart(2),xend(2)
      ! temp(1 ,j,i) = tempbp(j,i)
      ! temp(nx,j,i) = temptp(j,i)
      ! end do
      ! end do

      !assign linear salinity profile in the nodes k=1 to k=nxmr
      do i=xstartr(3),xendr(3)
      do j=xstartr(2),xendr(2)
      do k=1,nxmr
      sal(k,j,i)=salbp(j,i)-(salbp(j,i)-saltp(j,i))*xmr(k)/xcr(nxr)
      enddo
      enddo
      enddo

      !assign the boundary conditions at k=1 and k=nx
      ! do i=xstartr(3),xendr(3)
      ! do j=xstartr(2),xendr(2)
      ! sal(1  ,j,i) = salbp(j,i)
      ! sal(nxr,j,i) = saltp(j,i)
      ! end do
      ! end do

      return
      end
