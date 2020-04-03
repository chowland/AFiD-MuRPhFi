!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: CalcPlateCf.F90                                !
!    CONTAINS: subroutine CalcPlateCf                     !
!                                                         ! 
!    PURPOSE: Calculate the skin friction at the top      !
!     and bottom plates and output to a file.             !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine CalcPlateCf
      use param
      use local_arrays, only: vy
      use mpih
      use decomp_2d, only: xstart,xend
      implicit none
      integer :: j,i
      real ::  cflow, cfupp
      real :: del,deln
  

      cflow = 0.d0
      cfupp = 0.d0
      del  = 1.0/(xm(1)-xc(1))
      deln = 1.0/(xc(nx)-xm(nxm))

!$OMP  PARALLEL DO &
!$OMP   DEFAULT(none) &
!$OMP   SHARED(xstart,xend,vy,del,deln) &
!$OMP   SHARED(nxm,nx) &
!$OMP   PRIVATE(i,j) &
!$OMP   REDUCTION(+:cflow) &
!$OMP   REDUCTION(+:cfupp)
      do i=xstart(3),xend(3)
         do j=xstart(2),xend(2)
           cflow = cflow + vy(1,j,i)
           cfupp = cfupp + (1.d0-vy(nxm,j,i))
        enddo
      end do
!$OMP END PARALLEL DO

      !-- Cf = 2*u_\tau^2 = 2*(\nu*du/dz|_w)
      cflow = 2.d0 * (cflow / float(nzm*nym)) * del  / ren
      cfupp = 2.d0 * (cfupp / float(nzm*nym)) * deln / ren

      call MpiSumRealScalar(cflow)
      call MpiSumRealScalar(cfupp)

      if(ismaster) then
       open(98,file="outputdir/cf_plate.out",status='unknown', &
        access='sequential',position='append')
       write(98,546) time, cflow, cfupp, dsqrt(cflow/2.d0)*0.5d0*ren, dsqrt(cfupp/2.d0)*0.5d0*ren
 546   format(5(1x,e14.6))
       close(98)
      endif

      return         
      end                                                               
