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
      use local_arrays, only: sal,hsal,rhsr,rusal
      use decomp_2d, only: xstartr,xendr
      implicit none
      integer :: jc,kc,ic
      integer :: km,kp
      real    :: alpec,dxxs
      real    :: app,acc,amm

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
                        +sal(kc  ,jc,ic)*ac3sskr(kc) &
                        -(ap3sskr(kc)+ac3sskr(kc))*salbp(1,jc,ic)*SfixS
            elseif(kc.eq.nxmr) then !CJH Apply upper BC
                  dxxs = sal(kc  ,jc,ic)*ac3sskr(kc) &
                        +sal(kc-1,jc,ic)*am3sskr(kc) &
                        -(am3sskr(kc)+ac3sskr(kc))*saltp(1,jc,ic)*SfixN
            else
               dxxs= sal(kc+1,jc,ic)*ap3sskr(kc) &
                    +sal(kc  ,jc,ic)*ac3sskr(kc) &
                    +sal(kc-1,jc,ic)*am3sskr(kc)
            end if


!    Calculate right hand side of Eq. 5 (VO96)

            rhsr(kc,jc,ic)=(ga*hsal(kc,jc,ic)+ro*rusal(kc,jc,ic) &
                    +alpec*dxxs)*dt

!    Store the non-linear terms for the calculation of 
!    the next timestep

            rusal(kc,jc,ic)=hsal(kc,jc,ic)

        enddo
       enddo
      enddo
!$OMP END PARALLEL DO


!  Solve equation and update salinity

      call SolveImpEqnUpdate_Sal

!  Set boundary conditions on the salinity field at top
!  and bottom plates. This seems necessary.

      !  sal(1,xstartr(2):xendr(2),xstartr(3):xendr(3)) &
      !     = salbp(1,xstartr(2):xendr(2),xstartr(3):xendr(3))

      !  sal(nxr,xstartr(2):xendr(2),xstartr(3):xendr(3)) &
      !     = saltp(1,xstartr(2):xendr(2),xstartr(3):xendr(3))


      return
      end
