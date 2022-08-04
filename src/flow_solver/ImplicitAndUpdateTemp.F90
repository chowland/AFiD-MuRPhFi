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
    integer :: jc,kc,ic
    integer :: km,kp,ke
    real    :: alpec,dxxt
    real    :: app,acc,amm
    real    :: usaldto,dense

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
    do ic=xstart(3),xend(3)
        do jc=xstart(2),xend(2)
            do kc=1,nxm

!   Calculate second derivative of temperature in the x-direction.
!   This is the only term calculated implicitly for temperature.
                if (kc.eq.1) then       !CJH Apply lower BC
                    dxxt = temp(kc+1,jc,ic)*ap3ssk(kc) &
                        + temp(kc,jc,ic)*ac3ssk(kc) &
                        - (ap3ssk(kc)+ac3ssk(kc))*tempbp(1,jc,ic)*TfixS
                elseif(kc.eq.nxm) then  !CJH Apply upper BC
                    dxxt = temp(kc,jc,ic)*ac3ssk(kc) &
                        + temp(kc-1,jc,ic)*am3ssk(kc) &
                        - (am3ssk(kc)+ac3ssk(kc))*temptp(1,jc,ic)*TfixN
                else
                    dxxt = temp(kc+1,jc,ic)*ap3ssk(kc) &
                        + temp(kc  ,jc,ic)*ac3ssk(kc) &
                        + temp(kc-1,jc,ic)*am3ssk(kc)
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

    if (IBM) then
        forclo = 1.d0
        usaldto = 1.0/aldto
        do n=1,npunte
            ic = indgeot(n,1)
            jc = indgeot(n,2)
            kc = indgeot(n,3)
            forclo(kc,jc,ic) = 0.d0
            ke = indgeoet(n,3)
            ! dense = ((al*dt + aldto)*temp(ke,jc,ic) - al*dt*densb(n))*usaldto
            ! rhs(kc,jc,ic) = -temp(kc,jc,ic) + dense*distbt(n) &
            !                 + (1.0 - distbt(n))*temb(n)
            rhs(kc,jc,ic) = distbt(n)*temp(ke,jc,ic) - temp(kc,jc,ic) &
                            + (1.0 - distbt(n))*temb(n)
            hro(kc,jc,ic) = distbt(n)
            ! densb(n)= temp(ke,jc,ic)
        end do
        call SolveImpEqnUpdate_Temp_ibm
    else
!  Solve equation and update temperature

        call SolveImpEqnUpdate_Temp
    end if

    return
end subroutine ImplicitAndUpdateTemp