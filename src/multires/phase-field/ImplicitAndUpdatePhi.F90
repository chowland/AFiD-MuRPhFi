!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: ImplicitAndUpdatePhi.F90                       !
!    CONTAINS: subroutine ImplicitAndUpdatePhi            !
!                                                         ! 
!    PURPOSE: Compute the linear terms associated to      !
!     the phase-field and call the implicit solver.       !
!     After this routine, the phase-field has been        !
!     updated to the new timestep                         !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ImplicitAndUpdatePhi
    use param
    use mgrd_arrays, only: phi,hphi,rhsr,ruphi
    use decomp_2d, only: xstartr,xendr,nrank
    use mpih
    implicit none
    integer :: jc,kc,ic
    integer :: km,kp
    real    :: alpec,dxxp
    real    :: app,acc,amm

    alpec=al*pf_D

    do ic=xstartr(3),xendr(3)
        do jc=xstartr(2),xendr(2)
            do kc=1,nxmr

!   Calculate second derivative of phase field variable in the x-direction.
!   This is the only term calculated implicitly for phi.
                if (kc.eq.1) then       !CJH Apply lower BC d/dx(phi)=0
                    dxxp = phi(kc+1,jc, ic)*ap3spkr(kc) &
                          +phi(kc  ,jc, ic)*ac3spkr(kc)
                elseif(kc.eq.nxmr) then !CJH Apply upper BC d/dx(phi)=0
                    dxxp = phi(kc  ,jc,ic)*ac3spkr(kc) &
                          +phi(kc-1,jc,ic)*am3spkr(kc)
                else
                    dxxp= phi(kc+1,jc,ic)*ap3spkr(kc) &
                         +phi(kc  ,jc,ic)*ac3spkr(kc) &
                         +phi(kc-1,jc,ic)*am3spkr(kc)
                end if


!    Calculate right hand side of Eq. 5 (VO96)

                rhsr(kc,jc,ic)=(ga*hphi(kc,jc,ic) + ro*ruphi(kc,jc,ic) &
                                + alpec*dxxp)*dt

!    Store the non-linear terms for the calculation of 
!    the next timestep

                ruphi(kc,jc,ic)=hphi(kc,jc,ic)

            enddo
        enddo
    enddo
    
    ! iF(ANY(IsNaN(rhsr))) then
    !     write(*,*)nrank,'NaN in rhsr pre-solve'
    !     write(*,*)'Please try to reduce dtmax in bou.in'
    !     call MPI_Abort(MPI_COMM_WORLD,1,ierr)
    ! end if

!  Solve equation and update salinity

    call SolveImpEqnUpdate_Phi

    !Prevent phase field going outside of realistic values
    ! do ic=xstartr(3),xendr(3)
    !     do jc=xstartr(2),xendr(2)
    !         do kc=1,nxmr
    !             if (phi(kc,jc,ic).lt.0.0) then
    !                 phi(kc,jc,ic) = 0.0
    !             elseif (phi(kc,jc,ic).gt.1.0) then
    !                 phi(kc,jc,ic) = 1.0
    !             end if
    !         end do
    !     end do
    ! end do

    return
end subroutine ImplicitAndUpdatePhi