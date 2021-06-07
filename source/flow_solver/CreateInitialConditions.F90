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
      use local_arrays, only: vy,vx,temp,vz
      ! use mgrd_arrays, only: sal
      use decomp_2d, only: xstart,xend!,xstartr,xendr
      use mpih
      implicit none
      integer :: j,k,i,kmid
      real :: xxx,yyy,eps,varptb,amp

      call random_seed()
      eps=1d-2
      do i=xstart(3),xend(3)
        do j=xstart(2),xend(2)
          do k=1,nxm
            call random_number(varptb)
            vz(k,j,i)=eps*(2.d0*varptb - 1.d0)
            yyy=xm(k) 
            xxx=yc(j)            
            call random_number(varptb)
          !  varptb=0.d0
          !  vy(k,j,i)=-yyy*(1.d0-yyy)*(2.d0*yyy-1.d0)*sin(2.d0*pi*xxx/ylen)/pi
          !  vy(k,j,i)=(xm(k))+0.3*(2.d0*varptb-1.d0)+ &
!           vy(k,j,i)=+0.3*(2.d0*varptb-1.d0)+ &
!     &                  (2.0d0*yyy-6.0d0*yyy**2+4.0d0*yyy**3) &
!     &                  *sin(3*xxx)*eps

          !  yyy=xc(k)          
          !  xxx=ym(j)
          !  vx(k,j,i)=-yyy**2*(1.0d0-yyy)**2*cos(2.d0*pi*xxx/ylen)*eps
!           vx(k,j,i)=1.d0 ! -yyy**2*(1.0d0-yyy)**2*cos(3.1*xxx)*eps
            !CJH Laminar vertical convection as Batchelor (54) + noise
            vy(k,j,i) = min(ren,320.0)/12.d0*yyy*(2*yyy-1)*(yyy-1) + eps*(2.d0*varptb - 1.d0)
            call random_number(varptb)
            vx(k,j,i) = eps*(2.d0*varptb - 1.d0)
          enddo
          vx(1,j,i) = 0.d0   !CJH Enforce no-slip BC at x=0
        enddo
      enddo

      if (melt .or. phasefield) then
        kmid = nxm/2
        do i=xstart(3),xend(3)
          do j=xstart(2),xend(2)
            do k=1,nxm
              vx(k,j,i) = 0.d0
              vy(k,j,i) = 0.d0
              vz(k,j,i) = 0.d0
            end do
            ! call random_number(varptb)
            ! vy(1,j,i) = eps*(2.d0*varptb - 1.d0)
            ! do k=2,kmid
            !   call random_number(varptb)
            !   vx(k,j,i) = eps*(2.d0*varptb - 1.d0)
            !   call random_number(varptb)
            !   vy(k,j,i) = eps*(2.d0*varptb - 1.d0)
            ! end do
          end do
        end do
      end if

      !assign linear temperature profile in the nodes k=1 to k=nxm
      eps=0.1 !0.d0!5d-2
      kmid = nxm/2
      do i=xstart(3),xend(3)
      do j=xstart(2),xend(2)
      do k=1,nxm
        call random_number(varptb)
        temp(k,j,i)=tempbp(1,j,i)+(temptp(1,j,i)-tempbp(1,j,i))*xm(k)/xc(nx)
        ! if (abs(xm(k)-0.5) + eps > 0.5) then
        !   amp = 0.5 - abs(xm(k)-0.5) ! CJH Prevent values of |T| exceeding 0.5
        !   temp(k,j,i) = temp(k,j,i) + amp*(2.d0*varptb - 1.d0)
        ! else
        !   temp(k,j,i) = temp(k,j,i) + eps*(2.d0*varptb - 1.d0)
        ! end if
      end do
      !CJH Favier et al (2019) Appendix A3 Validation Case
      do k=1,kmid
        temp(k,j,i) = temp(k,j,i) &
              + eps*sin(4*pi*ym(j))*cos(4*pi*zm(j))*sin(2*pi*xm(k))**2
      enddo
      enddo
      enddo

      if (melt) then
        do i=xstart(3),xend(3)
          do j=xstart(2),xend(2)
            do k=1,nxm
              call random_number(varptb)
              temp(k,j,i) = eps*(2.d0*varptb - 1.d0) * exp(-xm(k)/0.1)
            end do
          end do
        end do
      end if

      ! FAVIER ET AL. (2019) APPENDIX A1 VALIDATION CASE
      ! if (phasefield) then
      !   ! kmid = nxm/2
      !   eps = 8.041
      !   do i=xstart(3),xend(3)
      !     do j=xstart(2),xend(2)
      !       do k=1,nxm
      !         temp(k,j,i) = (exp(-eps*(xm(k) - 1.0)) - 1.0)/(exp(eps) - 1.0)
      !       end do
      !       ! do k=1,kmid
      !       !   temp(k,j,i) = 1.0 - 2.0*xm(k)
      !       ! end do
      !       ! do k=kmid+1,nxm
      !       !   temp(k,j,i) = 0.1 - 0.2*xm(k)
      !       ! end do
      !     end do
      !   end do
      ! end if

      !assign the boundary conditions at k=1 and k=nx
      ! do i=xstart(3),xend(3)
      ! do j=xstart(2),xend(2)
      ! temp(1 ,j,i) = tempbp(1,j,i)
      ! temp(nx,j,i) = temptp(1,j,i)
      ! end do
      ! end do

      !assign linear salinity profile in the nodes k=1 to k=nxmr
      ! do i=xstartr(3),xendr(3)
      ! do j=xstartr(2),xendr(2)
      ! do k=1,nxmr
      !   call random_number(varptb)
      !   sal(k,j,i)=salbp(1,j,i)-(salbp(1,j,i)-saltp(1,j,i))*xmr(k)/xcr(nxr)
      !   if (abs(xmr(k)-0.5) + eps > 0.5) then
      !     amp = 0.5 - abs(xmr(k)-0.5) ! CJH Prevent values of |S| exceeding 0.5
      !     sal(k,j,i) = sal(k,j,i) + amp*(2.d0*varptb - 1.d0)
      !   else
      !     sal(k,j,i) = sal(k,j,i) + eps*(2.d0*varptb - 1.d0)
      !   end if
      ! enddo
      ! enddo
      ! enddo

      ! if (melt) then
      !   do i=xstartr(3),xendr(3)
      !     do j=xstartr(2),xendr(2)
      !       do k=1,nxmr
      !         call random_number(varptb)
      !         sal(k,j,i) = eps*(2.d0*varptb - 1.d0) * exp(-xmr(k)/0.1)
      !       end do
      !     end do
      !   end do
      ! end if

      !assign the boundary conditions at k=1 and k=nx
      ! do i=xstartr(3),xendr(3)
      ! do j=xstartr(2),xendr(2)
      ! sal(1  ,j,i) = salbp(1,j,i)
      ! sal(nxr,j,i) = saltp(1,j,i)
      ! end do
      ! end do

      return
      end
