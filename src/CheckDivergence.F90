!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: CheckDivergence.F90                            !
!    CONTAINS: subroutine CheckDivergence                 !
!                                                         ! 
!    PURPOSE: Check the maximum divergence of velocity    !
!     in the domain                                       !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine CheckDivergence(qmax,qmaxr)
      use param
      use local_arrays, only: vy,vx,vz
      use mgrd_arrays, only: vyr,vxr,vzr
      use mpih
      use decomp_2d, only: xstart,xend,xstartr,xendr
      implicit none
      real,intent(out) :: qmax,qmaxr
      integer :: jc,kc,kp,jp,ic,ip
      real    :: dqcap
        
      qmax =-huge(0.0)

!$OMP  PARALLEL DO &
!$OMP   DEFAULT(none) &
!$OMP   SHARED(xstart,xend,nxm,vz,vy,vx,dz,dy,udx3m) &
!$OMP   PRIVATE(ic,jc,kc,ip,jp,kp) &
!$OMP   PRIVATE(dqcap) &
!$OMP   REDUCTION(max:qmax)
      do ic=xstart(3),xend(3)
        ip=ic+1
        do jc=xstart(2),xend(2)
          jp=jc+1
            do kc=1,nxm
            kp=kc+1
              dqcap= (vz(kc,jc,ip)-vz(kc,jc,ic))*dz &
                    +(vy(kc,jp,ic)-vy(kc,jc,ic))*dy &
                    +(vx(kp,jc,ic)-vx(kc,jc,ic))*udx3m(kc)
              qmax = max(abs(dqcap),qmax)          
      enddo
      enddo
      enddo
!$OMP END PARALLEL DO

      call MpiMaxRealScalar(qmax)

!------------------------------------
      qmaxr =-huge(0.0)
! !$OMP  PARALLEL DO &
! !$OMP   DEFAULT(none) &
! !$OMP   SHARED(xstartr,xendr,nxmr,vzr,vyr,vxr,dzr,dyr,udx3mr) &
! !$OMP   PRIVATE(ic,jc,kc,ip,jp,kp) &
! !$OMP   PRIVATE(dqcap) &
! !$OMP   REDUCTION(max:qmaxr)
      do ic=xstartr(3),xendr(3)
        ip=ic+1
        do jc=xstartr(2),xendr(2)
          jp=jc+1
            do kc=1,nxmr
            kp=kc+1
              dqcap= (vzr(kc,jc,ip)-vzr(kc,jc,ic))*dzr &
                    +(vyr(kc,jp,ic)-vyr(kc,jc,ic))*dyr &
                    +(vxr(kp,jc,ic)-vxr(kc,jc,ic))*udx3mr(kc)
              qmaxr = max(abs(dqcap),qmaxr)
      enddo
      enddo
      enddo
! !$OMP END PARALLEL DO

      call MpiMaxRealScalar(qmaxr)
    
      
      return     
      end         
