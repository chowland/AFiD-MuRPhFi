!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: SolveImpEqnUpdate_Phi.F90                      !
!    CONTAINS: subroutine SolveImpEqnUpdate_Phi           !
!                                                         ! 
!    PURPOSE: Inverts the implicit equation for           !
!     salinity, and updates it to time t+dt               !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine SolveImpEqnUpdate_Phi
    use param
    use mgrd_arrays, only: phi,rhsr
    use decomp_2d, only: xstartr,xendr,nrank
    implicit none
    real, dimension(nxr) :: amkl,apkl,ackl
    integer :: jc,kc,info,ipkv(nxr),ic,nrhs
    real :: betadx,ackl_b
    real :: amkT(nxmr-1),ackT(nxmr),apkT(nxmr-1),appk(nxmr-2)

!   Calculate the coefficients of the tridiagonal matrix
!   The coefficients are normalized to prevent floating
!   point errors.

    betadx=0.5d0*al*dt*pf_D

    do kc=1,nxmr
        ackl_b=1.0d0/(1.-ac3sskr(kc)*betadx)
        amkl(kc)=-am3sskr(kc)*betadx*ackl_b
        ackl(kc)=1.0d0
        apkl(kc)=-ap3sskr(kc)*betadx*ackl_b
    end do

    amkT=amkl(2:nxmr)
    apkT=apkl(1:(nxmr-1))
    ackT=ackl(1:nxmr)

!   Call to LAPACK library to factor tridiagonal matrix.
!   No solving is done in this call.

    call dgttrf(nxmr,amkT,ackT,apkT,appk,ipkv,info)
    
    nrhs=(xendr(3)-xstartr(3)+1)*(xendr(2)-xstartr(2)+1)
    do ic=xstartr(3),xendr(3)
        do jc=xstartr(2),xendr(2)
            do kc=1,nxmr
                ackl_b=1.0/(1.0-ac3sskr(kc)*betadx)
                rhsr(kc,jc,ic)=rhsr(kc,jc,ic)*ackl_b
            end do
        end do
    end do

    call dgttrs('N',nxmr,nrhs,amkT,ackT,apkT,appk,ipkv,rhsr,nxmr,info)

    do ic=xstartr(3),xendr(3)
        do jc=xstartr(2),xendr(2)
            do kc=1,nxmr
                phi(kc,jc,ic) = phi(kc,jc,ic) + rhsr(kc,jc,ic)
            end do
        end do
    end do

    return
end subroutine SolveImpEqnUpdate_Phi
