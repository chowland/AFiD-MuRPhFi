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

      alpec=al/pec

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
      do kc=2,nxmr

!   Calculate second derivative of salinty in the x-direction.
!   This is the only term calculated implicitly for salinity.

               dxxs= sal(kc+1,jc,ic)*ap3ckr(kc) &
                    +sal(kc  ,jc,ic)*ac3ckr(kc) &
                    +sal(kc-1,jc,ic)*am3ckr(kc)


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

       sal(1,xstartr(2):xendr(2),xstartr(3):xendr(3)) &
          = salbp(xstartr(2):xendr(2),xstartr(3):xendr(3))

       sal(nx,xstartr(2):xendr(2),xstartr(3):xendr(3)) &
          = saltp(xstartr(2):xendr(2),xstartr(3):xendr(3))


      return
      end
