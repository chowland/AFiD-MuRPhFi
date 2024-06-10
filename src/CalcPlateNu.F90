
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
      use local_arrays, only: temp,vx,vy,vz,temp
      use mpih
      use decomp_2d, only: xstart,xend,xstartr,xendr
      use afid_salinity
      implicit none
      integer :: j,i,ip,jp,k,kp
      real ::  nu_T_Col, nu_T_Hot,  nu_S_Col, nu_S_Hot,Urms
      real, dimension(nxm) :: vxrms, vyrms, vzrms,u_rms
      real:: del,deln,delr,delnr
      real :: inym, inzm
      integer :: unit

      real ::count
      logical :: file_exists
      character(len=1024) :: buffer
      vxrms(:)=0.0;   vyrms(:)=0.0;   vzrms(:)=0.0

      del  = 1.0/(xm(1)-xc(1))
      deln = 1.0/(xc(nx)-xm(nxm))
      if(salinity)then
      delr  = 1.0/(xmr(1)-xcr(1))
      delnr = 1.0/(xcr(nx)-xmr(nxm))
      endif 
    
   

      inym = 1.d0/nym
      inzm = 1.d0/nzm

!$OMP  PARALLEL DO &
!$OMP   DEFAULT(none) &
!$OMP   SHARED(xstart,xend,temp,del,deln) &
!$OMP   SHARED(nxm,nx) &
!$OMP   PRIVATE(i,j) &
!$OMP   REDUCTION(+:nuslow) &
!$OMP   REDUCTION(+:nusupp)
      !do j=1,nym
      !end do
      do i=xstart(3),xend(3)
            ip = i + 1
            do j=xstart(2),xend(2) 
                  jp = j + 1
                  do k=1,nxm
                        kp = k + 1
                         vxrms(k) = vxrms(k) + 0.5*(vx(k,j,i)**2+vx(kp,j,i)**2)
                         vyrms(k) = vyrms(k) + 0.5*(vy(k,j,i)**2+vy(k,jp,i)**2)
                         vzrms(k) = vzrms(k) + 0.5*(vz(k,j,i)**2+vz(k,j,ip)**2)
                  end do
            if  (FixValueBCRegion_Length/=0 .and.FixValueBCRegion_Nord_or_Sud==1) then
                  if (ym(j) <= 0.01 * FixValueBCRegion_Length * YLEN) then
                        nu_T_Hot = nu_T_Hot + (temptp(1,j,i)-temp(nxm,j,i))*deln  

                  else if ( ym(j) >= YLEN - 0.01 * FixValueBCRegion_Length * YLEN) then
                        nu_T_Col = nu_T_Col + (temptp(1,j,i)-temp(nxm,j,i))*deln
                  end if 
            else if  (FixValueBCRegion_Length/=0 .and.FixValueBCRegion_Nord_or_Sud==0) then
                  if (ym(j) <= 0.01 * FixValueBCRegion_Length * YLEN) then
                        nu_T_Hot = nu_T_Hot + (temp(1,j,i)-tempbp(1,j,i))*del
                  else if ( ym(j) >= YLEN - 0.01 * FixValueBCRegion_Length * YLEN) then
                        nu_T_Col = nu_T_Col + (temp(1,j,i)-tempbp(1,j,i))*del
                  end if
            end if

        enddo
      end do


      nu_T_Hot = nu_T_Hot / (nzm*ny_Hot)
      nu_T_Col = nu_T_Col / (nzm*ny_Cold)

      call MpiSumRealScalar(nu_T_Hot)
      call MpiSumRealScalar(nu_T_Col)
 if(salinity) then 
            do i=xstartr(3),xendr(3)
                  do j=xstartr(2),xendr(2) 

                  if  (FixValueBCRegion_Length/=0 .and.FixValueBCRegion_Nord_or_Sud==1) then
                        if (ymr(j) <= 0.01 * FixValueBCRegion_Length * YLEN) then
                              nu_S_Hot = nu_S_Hot + (saltp(1,j,i)-sal(nxmr,j,i))*delnr
                        else if ( ymr(j) >= YLEN - 0.01 * FixValueBCRegion_Length * YLEN) then
                              nu_S_Col = nu_S_Col + (saltp(1,j,i)-sal(nxmr,j,i))*delnr

                        end if 
                  else if  (FixValueBCRegion_Length/=0 .and.FixValueBCRegion_Nord_or_Sud==0) then
                        if (ymr(j) <= 0.01 * FixValueBCRegion_Length * YLEN) then
                              nu_S_Hot = nu_S_Hot + (sal(1,j,i)-salbp(1,j,i))*delr

                        else if ( ymr(j) >= YLEN - 0.01 * FixValueBCRegion_Length * YLEN) then
                              nu_S_Col = nu_S_Col + (sal(1,j,i)-salbp(1,j,i))*delr

                        end if
                  end if
      
              enddo
            end do


      nu_S_Hot = nu_S_Hot / (nzm*nyr_Hot)
      nu_S_Col = nu_S_Col / (nzm*nyr_Cold)

      call MpiSumRealScalar(nu_S_Hot)
      call MpiSumRealScalar(nu_S_Col)

end if 


do k=1,nxm
      vxrms(k) = sqrt(vxrms(k)*inym*inzm)
      vyrms(k) = sqrt(vyrms(k)*inym*inzm)
      vzrms(k) = sqrt(vzrms(k)*inym*inzm)
  end do

  u_rms = vxrms**2 + vyrms**2 +  vyrms**2
  u_rms = u_rms/nxm
  Urms = sum(u_rms)

!$OMP END PARALLEL DO

  if (ismaster) then
      ! Verifica se il file esiste
      inquire(file="outputdir/nu_plate.csv", exist=file_exists)
      if (.not. file_exists) then
          ! Se il file non esiste, crealo e scrivi la legenda
          open(newunit=unit, file="outputdir/nu_plate.csv", status='new', &
               access='sequential')
               if (salinity) then
                  write(unit, '(A)') "           time,           Urms,      nu_T_Hot,      nu_T_Col,      nu_S_Hot,      nu_S_Col"
              else
                  write(unit, '(A)') "           time,           Urms,      nu_T_Hot,      nu_T_Col"
              end if
              close(unit)
      end if
      ! Apri il file in modalità append
      open(newunit=unit, file="outputdir/nu_plate.csv", status='old', &
           access='sequential', position='append')

      if (salinity) then
          write(unit, '(ES22.14,1x,ES22.14,1x,ES22.14,1x,ES22.14,1x,ES22.14,1x,ES22.14)') &
              time, Urms, nu_T_Hot, nu_T_Col, nu_S_Hot, nu_S_Col
      else
          write(unit, '(ES22.14,1x,ES22.14,1x,ES22.14,1x,ES22.14)') &
              time, Urms, nu_T_Hot, nu_T_Col
      end if

      close(unit)
  end if
  

  
  
      return         
      end                                                               
