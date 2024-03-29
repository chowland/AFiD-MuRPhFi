!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: SolveImpEqnUpdate_Temp.F90                     !
!    CONTAINS: subroutine SolveImpEqnUpdate_Temp          !
!                                                         ! 
!    PURPOSE: Inverts the implicit equation for           !
!     temperature, and updates it to time t+dt            !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine SolveImpEqnUpdate_Temp
    use param
    use local_arrays, only : temp,rhs
    use decomp_2d, only: xstart,xend
    implicit none
    real, dimension(nx) :: amkl,apkl,ackl
    integer :: jc,kc,info,ipkv(nxm),ic,nrhs
    real :: betadx,ackl_b
    real :: amkT(nxm-1),ackT(nxm),apkT(nxm-1),appk(nxm-2)

!     Calculate the coefficients of the tridiagonal matrix
!     The coefficients are normalized to prevent floating
!     point errors.

    betadx=0.5d0*al*dt/pect

    do kc=1,nxm
        ackl_b=1.0d0/(1.0d0-ac3ssk(kc)*betadx)
        amkl(kc)=-am3ssk(kc)*betadx*ackl_b
        ackl(kc)=1.0d0
        apkl(kc)=-ap3ssk(kc)*betadx*ackl_b
    end do

    amkT=amkl(2:nxm)
    apkT=apkl(1:(nxm-1))
    ackT=ackl(1:nxm)

!     Call to LAPACK library to factor tridiagonal matrix.
!     No solving is done in this call.

    call dgttrf(nxm,amkT,ackT,apkT,appk,ipkv,info)
    
    nrhs=(xend(3)-xstart(3)+1)*(xend(2)-xstart(2)+1)
    do ic=xstart(3),xend(3)
        do jc=xstart(2),xend(2)
            do kc=1,nxm
                ackl_b=1.0/(1.0-ac3ssk(kc)*betadx)
                rhs(kc,jc,ic)=rhs(kc,jc,ic)*ackl_b
            end do
        end do
    end do
      
    call dgttrs('N',nxm,nrhs,amkT,ackT,apkT,appk,ipkv,rhs(1:nxm,:,:),nxm,info)

    do ic=xstart(3),xend(3)
        do jc=xstart(2),xend(2)
            do kc=1,nxm
                temp(kc,jc,ic)=temp(kc,jc,ic) + rhs(kc,jc,ic)
            end do
        end do
    end do

    return
end subroutine SolveImpEqnUpdate_Temp