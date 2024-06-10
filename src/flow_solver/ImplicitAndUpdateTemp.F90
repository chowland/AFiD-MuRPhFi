!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: ImplicitAndUpdateTemp.F90                      !
!    CONTAINS: subroutine ImplicitAndUpdateTemp           !
!                                                         ! 
!    PURPOSE: Compute the linear terms associated to      !
!     the temperature and call the implicit solver.       !
!     After this routine, the temperature has been        !
!     updated to the new timestep                         !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ImplicitAndUpdateTemp
    use param
    use local_arrays, only: temp,hro,rutemp,rhs
    use decomp_2d, only: xstart,xend
    use ibm_param
    implicit none
    integer :: jc,kc,ic,ii
    real    :: alpec,dxxt, FlagBC_Nord,FlagBC_Sud
    alpec=al/pect
!$OMP  PARALLEL DO &
!$OMP   DEFAULT(none) &
!$OMP   SHARED(xstart,xend,nxm,temp) &
!$OMP   SHARED(kmv,kpv,am3ck,ac3ck,ap3ck) &
!$OMP   SHARED(ga,ro,alpec,dt) &
!$OMP   SHARED(rhs,rutemp,hro) &
!$OMP   PRIVATE(ic,jc,kc,km,kp) &
!$OMP   PRIVATE(amm,acc,app) &
!$OMP   PRIVATE(dxxt)
    if (FixValueBCRegion_Length==0) then
        FlagBC_Sud = TfixS
        FlagBC_Nord = TfixN
    else if  (FixValueBCRegion_Length/=0 .and. FixValueBCRegion_Nord_or_Sud==0) then
        FlagBC_Nord = TfixN
    else if  (FixValueBCRegion_Length/=0 .and. FixValueBCRegion_Nord_or_Sud==1) then
        FlagBC_Sud = TfixS
    end if

    do ic=xstart(3),xend(3)
        do jc=xstart(2),xend(2)
            do kc=1,nxm
                     if (FixValueBCRegion_Length/=0) then
                        if (FixValueBCRegion_Nord_or_Sud==0) then
                            if (ym(jc) <= 0.01 * FixValueBCRegion_Length * YLEN .or. &
                                ym(jc) >= YLEN - 0.01 * FixValueBCRegion_Length * YLEN) then
                                ii = 1
                                FlagBC_Sud = 1
                             else 
                                ii = 2
                                FlagBC_Sud = 0
                             end if
                        else if (FixValueBCRegion_Nord_or_Sud == 1) then
                            if (ym(jc) <= 0.01 * FixValueBCRegion_Length * YLEN .or. &
                                ym(jc) >= YLEN - 0.01 * FixValueBCRegion_Length * YLEN) then
                                
                                ii = 1
                                FlagBC_Nord = 1
                             else 
                                ii = 2
                                FlagBC_Nord = 0
                             end if
                        end if
                    else 
                    ii = 1
                    end if
!   Calculate second derivative of temperature in the x-direction.
!   This is the only term calculated implicitly for temperature.
                if (kc.eq.1) then       !CJH Apply lower BC
                    dxxt = temp(kc+1,jc,ic)*ap3ssk(kc,ii) &
                        + temp(kc,jc,ic)*ac3ssk(kc,ii) &
                        - (ap3ssk(kc,ii)+ac3ssk(kc,ii))*tempbp(1,jc,ic)*FlagBC_Sud
                       

                elseif(kc.eq.nxm) then  !CJH Apply upper BC
                    dxxt = temp(kc,jc,ic)*ac3ssk(kc,ii) &
                        + temp(kc-1,jc,ic)*am3ssk(kc,ii) &
                        - (am3ssk(kc,ii)+ac3ssk(kc,ii))*temptp(1,jc,ic)*FlagBC_Nord
                else
                    dxxt = temp(kc+1,jc,ic)*ap3ssk(kc,ii) &
                        + temp(kc  ,jc,ic)*ac3ssk(kc,ii) &
                        + temp(kc-1,jc,ic)*am3ssk(kc,ii)
                end if
   
!    Calculate right hand side of Eq. 5 (VO96)

                rhs(kc,jc,ic)=(ga*hro(kc,jc,ic)+ro*rutemp(kc,jc,ic) &
                        +alpec*dxxt)*dt

!    Store the non-linear terms for the calculation of 
!    the next timestep

                rutemp(kc,jc,ic)=hro(kc,jc,ic)

            enddo
        enddo
    enddo
!$OMP END PARALLEL DO

    if (IBM .and. .not. phasefield) then
        call SolveImpEqnUpdate_Temp_ibm
    else
!  Solve equation and update temperature

        call SolveImpEqnUpdate_Temp
    end if

    return
end subroutine ImplicitAndUpdateTemp