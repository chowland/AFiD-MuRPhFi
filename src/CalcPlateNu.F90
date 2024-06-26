
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
      use local_arrays, only: temp,vx,vy,vz
      use mpih
      use decomp_2d, only: xstart,xend,xstartr,xendr
      use afid_salinity, only: sal,salbp,&
                              &saltp,alpha_Sal
      use afid_Termperature_Fine
      implicit none
      integer :: j,i,ip,jp,k,kp
      real ::  nu_T_Col, nu_T_Hot,  nu_S_Col, nu_S_Hot,Urms
      real ::  nuTslow, nuTsupp, nuSslow, nuSsupp
      real :: Ghost_var
      real, dimension(nx-1) :: dx_mean
      real, dimension(nxm) :: vxrms, vyrms, vzrms,u_rms
      real:: del,deln,delr,delnr
      real :: inym, inzm
      integer :: unit
      real :: Ghost

      integer :: case_number
      real ::count
      logical :: file_exists
      character(len=1024) :: buffer
  


      vxrms(:)=0.0;   vyrms(:)=0.0;   vzrms(:)=0.0
      del  = 1.0/(xm(1)-xc(1))
      deln = 1.0/(xc(nx)-xm(nxm))
      if(salinity .or. multiRes_Temp)then
      delr  = 1.0/(xmr(1)-xcr(1))
      delnr = 1.0/(xcr(nxr)-xmr(nxmr))
      endif 

      if (FixValueBCRegion_Length /= 0) then
            case_number = 1
        else
            case_number = 2
        end if

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

      select case (case_number)
      case (1)
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
     if(.not.multiRes_Temp)then
                  if(Robin==1)then
                        if (ym(j) <= YLEN/2 .and. alpha_Temp(j)/=0) then
                              call T_G(alpha_Temp(j), deln, temp(nxm,j,i), temptp(1,j,i), Ghost)
                              nu_T_Hot =nu_T_Hot+(Ghost-temp(nxm,j,i)) *deln/2 
                              

                        elseif (ym(j) >= YLEN/2 .and. alpha_Temp(j)/=0) then
                              call T_G(alpha_Temp(j), deln, temp(nxm,j,i), temptp(1,j,i), Ghost)
                              nu_T_Col =nu_T_Col+(Ghost-temp(nxm,j,i)) *deln/2 

                                          
                        end if

                  else 
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
                  end if 
            end if 
        enddo
      end do

      if(.not.multiRes_Temp)then
      nu_T_Hot = nu_T_Hot / (nzm*ny_Hot)
      nu_T_Col = nu_T_Col / (nzm*ny_Cold)

      call MpiSumRealScalar(nu_T_Hot)
      call MpiSumRealScalar(nu_T_Col)
      end if 

 if(salinity) then 
            do i=xstartr(3),xendr(3)
                  do j=xstartr(2),xendr(2) 
                        if(Robin==1)then
                        if (ymr(j) <= YLEN/2 .and. alpha_Sal(j)/=0) then
                              call T_G(alpha_Sal(j), delnr, sal(nxmr,j,i), saltp(1,j,i), Ghost)
                              nu_S_Hot =nu_S_Hot+(Ghost-sal(nxmr,j,i)) *delnr/2 

                        elseif (ymr(j) >= YLEN/2 .and. alpha_Sal(j)/=0)then
                              call T_G(alpha_Sal(j), delnr, sal(nxmr,j,i), saltp(1,j,i), Ghost)
                              nu_S_Col =nu_S_Col+(Ghost-sal(nxmr,j,i)) *delnr/2 
                                          
                  end if
            end if 
            if(multiRes_Temp)then
                  if  (FixValueBCRegion_Length/=0 .and.FixValueBCRegion_Nord_or_Sud==1) then
                        if (ymr(j) <= 0.01 * FixValueBCRegion_Length * YLEN) then
                              nu_T_Hot = nu_T_Hot + (temp_fine_tp(1,j,i)-temp_fine(nxmr,j,i))*delnr
                        else if ( ymr(j) >= YLEN - 0.01 * FixValueBCRegion_Length * YLEN) then
                              nu_T_Col = nu_T_Col + (temp_fine_tp(1,j,i)-temp_fine(nxmr,j,i))*delnr

                        end if 
                  else if  (FixValueBCRegion_Length/=0 .and.FixValueBCRegion_Nord_or_Sud==0) then
                        if (ymr(j) <= 0.01 * FixValueBCRegion_Length * YLEN) then
                              nu_T_Hot = nu_T_Hot + (temp_fine_bp(1,j,i)-salbp(1,j,i))*delr

                        else if ( ymr(j) >= YLEN - 0.01 * FixValueBCRegion_Length * YLEN) then
                              nu_T_Col = nu_T_Col + (temp_fine_bp(1,j,i)-salbp(1,j,i))*delr

                        end if
                  end if
            end if 

              enddo
            end do


      nu_S_Hot = nu_S_Hot / (nzm*nyr_Hot)
      nu_S_Col = nu_S_Col / (nzm*nyr_Cold)

      call MpiSumRealScalar(nu_S_Hot)
      call MpiSumRealScalar(nu_S_Col)

      if(multiRes_Temp)then
            nu_T_Hot = nu_T_Hot / (nzm*nyr_Hot)
            nu_T_Col = nu_T_Col / (nzm*nyr_Cold)
      
            call MpiSumRealScalar(nu_T_Hot)
            call MpiSumRealScalar(nu_T_Col)

      end if 

end if 


do k=1,nxm
      vxrms(k) = sqrt(vxrms(k)*inym*inzm)
      vyrms(k) = sqrt(vyrms(k)*inym*inzm)
      vzrms(k) = sqrt(vzrms(k)*inym*inzm)
  end do

  u_rms = vxrms**2 + vyrms**2 +  vyrms**2
  do i = 1, nx-1
      dx_mean(i) = xc(i+1) - xc(i)
      Urms = Urms + u_rms(i) * dx_mean(i)
end do
write(*,*)Urms

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
case (2)
      do i=xstart(3),xend(3)
            ip = i + 1
            do j=xstart(2),xend(2)
                  jp = j + 1
                  do k=1,nxm
                         kp = k + 1
                               vxrms(k) = vxrms(k) + 0.5*(vx(k,j,i)**2+vx(kp,j,i)**2)
                               vyrms(k) = vyrms(k) + 0.5*(vy(k,j,i)**2+vy(k,jp,i)**2)
                               vzrms(k) = vzrms(k) + 0.5*(vz(k,j,i)**2+vz(k,j,ip)**2)
            if(.not. multiRes_Temp)then
              nuTslow = nuTslow + (temp(1,j,i)-tempbp(1,j,i))*del
              nuTsupp = nuTsupp + (temptp(1,j,i)-temp(nxm,j,i))*deln
            end if 
                  enddo
           enddo
         enddo
         if(.not. multiRes_Temp)then

         nuTslow = nuTslow / (nzm*nym)
         nuTsupp = nuTsupp / (nzm*nym)
   
         call MpiSumRealScalar(nuTslow)
         call MpiSumRealScalar(nuTsupp)
         end if 

      if(salinity) then
         do i=xstartr(3),xendr(3)
            do j=xstartr(2),xendr(2)
              nuSslow = nuSslow + (sal(1,j,i)-salbp(1,j,i))*delr
              nuSsupp = nuSsupp + (saltp(1,j,i)-sal(nxmr,j,i))*delnr

           enddo
         end do
         nuSslow = nuSslow / (nzmr*nymr)
         nuSsupp = nuSsupp / (nzmr*nymr)
   
         call MpiSumRealScalar(nuSslow)
         call MpiSumRealScalar(nuSsupp)
      end if 
         if(multiRes_Temp)then
            do i=xstartr(3),xendr(3)
               do j=xstartr(2),xendr(2)
                     nuTslow = nuTslow + (temp_fine(1,j,i)-temp_fine_bp(1,j,i))*delr
                     nuTsupp = nuTsupp + (temp_fine_tp(1,j,i)-temp_fine(nxmr,j,i))*delnr
              enddo
            end do

            nuTslow = nuTslow / (nzmr*nymr)
            nuTsupp = nuTsupp / (nzmr*nymr)
      
            call MpiSumRealScalar(nuTslow)
            call MpiSumRealScalar(nuTsupp)
         

      end if 

      if (ismaster) then
            ! Verifica se il file esiste
            inquire(file="outputdir/nu.csv", exist=file_exists)
            if (.not. file_exists) then
                ! Se il file non esiste, crealo e scrivi la legenda
                open(newunit=unit, file="outputdir/nu.csv", status='new', &
                     access='sequential')
                     if (salinity) then
                        write(unit, '(A)') &
                        "           time,           Urms,      nu_T_Up,      nu_T_Low,      nu_S_Up,      nu_S_Low"
                       else
                        write(unit, '(A)') "           time,           Urms,      nu_T_Up,      nu_T_Low"
                    end if
                    close(unit)
            end if
            ! Apri il file in modalità append
            open(newunit=unit, file="outputdir/nu.csv", status='old', &
                 access='sequential', position='append')
      
            if (salinity) then
                write(unit, '(ES22.14,1x,ES22.14,1x,ES22.14,1x,ES22.14,1x,ES22.14,1x,ES22.14)') &
                    time, Urms, nuTsupp, nuTslow,  nuSsupp, nuSslow
            else
                write(unit, '(ES22.14,1x,ES22.14,1x,ES22.14,1x,ES22.14)') &
                    time, Urms,  nuTsupp, nuTslow
            end if
      
            close(unit)
        end if
end select
  
  
      return         
      end   
      
      subroutine T_G(alpha, del, Var_nx, delta_Var, Ghost)
            implicit none
            real, intent(in):: alpha, Var_nx, delta_Var,del
            real, intent(out):: Ghost
            real :: term1, term2,delta_x

            delta_x = 1/del
            term1 = (-2.0 * alpha / (alpha + (1.0 - alpha) / delta_x)) * (Var_nx / 2.0 - delta_Var)
            term2 = ((1.0 - alpha) / (alpha + (1.0 - alpha) / delta_x)) * (Var_nx / delta_x)
            
            Ghost = term1 + term2  
        end subroutine T_G
        
