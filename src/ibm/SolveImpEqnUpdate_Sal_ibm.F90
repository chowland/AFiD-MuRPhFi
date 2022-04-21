!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: SolveImpEqnUpdate_Sal.F90                      !
!    CONTAINS: subroutine SolveImpEqnUpdate_Sal           !
!                                                         ! 
!    PURPOSE: Inverts the implicit equation for           !
!     salinity, and updates it to time t+dt               !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine SolveImpEqnUpdate_Sal_ibm
    use param
    use mgrd_arrays, only: sal,rhsr
    use decomp_2d, only: xstartr,xendr
    use ibm_param, only: forclor
    ! use param_particle      ! SL
    implicit none
    real, dimension(nxr) :: amkl,apkl,ackl,fkl
    integer :: jc,kc,info,ipkv(nxmr),ic,nrhs
    real :: betadx,ackl_b
    real :: amkT(nxmr-1),ackT(nxmr),apkT(nxmr-1),appk(nxmr-2)

!   Calculate the coefficients of the tridiagonal matrix
!   The coefficients are normalized to prevent floating
!   point errors.

    betadx=0.5d0*al*dt/pecs

    do ic=xstartr(3),xendr(3)
        do jc=xstartr(2),xendr(2)
            do kc=1,nxmr
                ackl_b = 1.0d0/(1.-ac3sskr(kc)*betadx*forclor(kc,jc,ic))
                amkl(kc) = -am3sskr(kc)*betadx*ackl_b
                ackl(kc) = 1.0d0
                apkl(kc) = -ap3sskr(kc)*betadx*ackl_b
                fkl(kc) = rhsr(kc,jc,ic)*ackl_b
            end do
            fkl(nxr) = 0.d0
            amkT = amkl(2:nxmr)
            apkT = apkl(1:(nxmr-1))
            ackT = ackl(1:nxmr)
            call dgttrf(nxmr,amkT,ackT,apkT,appk,ipkv,info)
            call dgttrs('N',nxmr,1,amkT,ackT,apkT,appk,ipkv,fkl,nxr,info)
            do kc=1,nxmr
                sal(kc,jc,ic) = sal(kc,jc,ic) + fkl(kc)
            end do
        end do
    end do
    return
end subroutine SolveImpEqnUpdate_Sal_ibm