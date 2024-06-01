!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: CalcPlateNu.F90                                !
!    CONTAINS: subroutine CalcPlateNu                     !
!                                                         ! 
!    PURPOSE: Calculate the Nusselt number at the top     !
!     and bottom plates and output to a file.             !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine CalcPlateNu
      use param
      use local_arrays, only: temp
      use mpih
      use decomp_2d, only: xstart,xend
      implicit none
      integer :: j,i
      real ::  nuslow, nusupp
      real :: del,deln
  

      nuslow = 0.d0
      nusupp = 0.d0
      del  = 1.0/(xm(1)-xc(1))
      deln = 1.0/(xc(nx)-xm(nxm))

!$OMP  PARALLEL DO &
!$OMP   DEFAULT(none) &
!$OMP   SHARED(xstart,xend,temp,del,deln) &
!$OMP   SHARED(nxm,nx) &
!$OMP   PRIVATE(i,j) &
!$OMP   REDUCTION(+:nuslow) &
!$OMP   REDUCTION(+:nusupp)
      do i=xstart(3),xend(3)
            do j=xstart(2),xend(2)
              
            if  (FixValueBCRegion_Length/=0 .and.FixValueBCRegion_Nord_or_Sud==0) then
                  if (yc(jc) < 0.01 * FixValueBCRegion_Length * YLEN
                       nusCol = nusCol + (temptp(1,j,i)-temp(nxm,j,i))*deln
                  else if ( yc(jc) > YLEN - 0.01 * FixValueBCRegion_Length * YLEN) then
                        nusHot = nusHot + (temptp(1,j,i)-temp(nxm,j,i))*deln   
                   end if 
            else if  (FixValueBCRegion_Length/=0 .and.FixValueBCRegion_Nord_or_Sud==1) then
                  if (yc(jc) < 0.01 * FixValueBCRegion_Length * YLEN
                      nusCol = nusCol + (temp(1,j,i)-tempbp(1,j,i))*del
                  else if ( yc(jc) > YLEN - 0.01 * FixValueBCRegion_Length * YLEN) then
                       nusHot = nusHot + (temp(1,j,i)-tempbp(1,j,i))*del

            end if

        enddo
      end do
!$OMP END PARALLEL DO

      nuslow = nuslow / (nzm*nym)
      nusupp = nusupp / (nzm*nym)

      call MpiSumRealScalar(nuslow)
      call MpiSumRealScalar(nusupp)

      if(ismaster) then
       open(97,file="outputdir/nu_plate.out",status='unknown', &
        access='sequential',position='append')
       write(97,546) time, nuslow, nusupp
 546   format(4(1x,e14.6))
       close(97)
      endif

      return         
      end                                                               
