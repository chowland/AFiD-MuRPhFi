!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: ImplicitAndUpdateSal.F90                       !
!    CONTAINS: subroutine ImplicitAndUpdateSal            !
!                                                         ! 
!    PURPOSE: Compute the linear terms associated to      !
!     the salinity and call the implicit solver.          !
!     After this routine, the salinity has been           !
!     updated to the new timestep                         !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ImplicitAndUpdateSal
    use param
    use mgrd_arrays, only: sal,hsal,rhsr,rusal
    use decomp_2d, only: xstartr,xendr
    use ibm_param
    implicit none
    integer :: jc,kc,ic
    integer :: km,kp,ke
    real    :: alpec,dxxs
    real    :: app,acc,amm
    real    :: usaldto,sale

    alpec=al/pecs

!$OMP  PARALLEL DO &
!$OMP   DEFAULT(none) &
!$OMP   SHARED(xstartr,xendr,nxmr,sal) &
!$OMP   SHARED(kmv,kpv,am3ckr,ac3ckr,ap3ckr) &
!$OMP   SHARED(ga,ro,alpec,dt) &
!$OMP   SHARED(rhsr,rusal,hsal) &
!$OMP   PRIVATE(ic,jc,kc,km,kp) &
!$OMP   PRIVATE(amm,acc,app) &
!$OMP   PRIVATE(dxxs)
    do ic=xstartr(3),xendr(3)
        do jc=xstartr(2),xendr(2)
            do kc=1,nxmr

!   Calculate second derivative of salinty in the x-direction.
!   This is the only term calculated implicitly for salinity.
                if (kc.eq.1) then       !CJH Apply lower BC
                    dxxs = sal(kc+1,jc,ic)*ap3sskr(kc) &
                            + sal(kc,jc,ic)*ac3sskr(kc) &
                            - (ap3sskr(kc) + ac3sskr(kc))*salbp(1,jc,ic)*SfixS
                elseif(kc.eq.nxmr) then !CJH Apply upper BC
                    dxxs = sal(kc,jc,ic)*ac3sskr(kc) &
                            + sal(kc-1,jc,ic)*am3sskr(kc) &
                            - (am3sskr(kc) + ac3sskr(kc))*saltp(1,jc,ic)*SfixN
                else
                    dxxs = sal(kc+1,jc,ic)*ap3sskr(kc) &
                            + sal(kc  ,jc,ic)*ac3sskr(kc) &
                            + sal(kc-1,jc,ic)*am3sskr(kc)
                end if


!    Calculate right hand side of Eq. 5 (VO96)

                rhsr(kc,jc,ic) = (ga*hsal(kc,jc,ic) + ro*rusal(kc,jc,ic) &
                        + alpec*dxxs)*dt

!    Store the non-linear terms for the calculation of 
!    the next timestep

                rusal(kc,jc,ic) = hsal(kc,jc,ic)

            enddo
        enddo
    enddo
!$OMP END PARALLEL DO


!  Solve equation and update salinity

    if (IBM) then
        forclor = 1.d0
        usaldto = 1.0/aldto
        do n=1,npuntr
            ic = indgeor(n,1)
            jc = indgeor(n,2)
            kc = indgeor(n,3)
            forclor(kc,jc,ic) = 0.d0
            ke = indgeoer(n,3)
            ! sale = ((al*dt + aldto)*sal(ke,jc,ic) - al*dt*salb(n))*usaldto
            ! rhsr(kc,jc,ic) = -sal(kc,jc,ic) + sale*distbr(n) !&
                            ! + (1.0 - distbr(n))*salfix(n)
            rhsr(kc,jc,ic) = sal(ke,jc,ic) - sal(kc,jc,ic)
            ! salb(n) = sal(ke,jc,ic)
        end do
        call SolveImpEqnUpdate_Sal_ibm
    else
        call SolveImpEqnUpdate_Sal
    end if

    return
end subroutine ImplicitAndUpdateSal