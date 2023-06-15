!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: SolveImpEqnUpdate_X.F90                        !
!    CONTAINS: subroutine SolveImpEqnUpdate_X             !
!                                                         ! 
!    PURPOSE: Inverts the implicit equation for velocity  !
!     in any the vertical direction, and updates it to    !
!     time t+dt                                           !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine SolveImpEqnUpdate_X_ibm
    use param
    use local_arrays, only : vx,rhs
    use decomp_2d, only: xstart,xend
    use ibm_param, only: ibmaskx, distx
    implicit none
    real, dimension(nx) :: amkl,apkl,ackl, fkl
    real :: amkT(nx-1),apkT(nx-1)
    real :: appk(nx-2)
    real :: ackT(nx)
    integer :: jc,kc,info,ic,km,kp,n
    integer :: ipkv(nx)
    real :: betadx,ackl_b

    betadx=beta*al

    amkl(1)=0.d0
    apkl(1)=0.d0
    
    amkl(nx)=0.d0
    apkl(nx)=0.d0
    ackl(:)=1.d0

    n = 1

    do ic=xstart(3),xend(3)
        do jc=xstart(2),xend(2)
            fkl(1)= 0.d0
            do kc=1,nxm
                km = max(1,kc - 1)
                kp = kc + 1
                if (ibmaskx(kc,jc,ic) == 2) then ! Liquid phase
                    ackl_b=1.0d0/(1.0d0-ac3ck(kc)*betadx)
                    amkl(kc)=-am3ck(kc)*betadx*ackl_b
                    ackl(kc)=1.d0
                    apkl(kc)=-ap3ck(kc)*betadx*ackl_b
                    fkl(kc) = rhs(kc,jc,ic)*ackl_b
                elseif (ibmaskx(kc,jc,ic) == 0) then ! Solid phase
                    amkl(kc) = 0.d0
                    ackl(kc) = 1.d0
                    apkl(kc) = 0.d0
                    fkl(kc) = -vx(kc,jc,ic)
                elseif (ibmaskx(kc,jc,ic) == 1) then ! Upper boundary points
                    amkl(kc) = 0.d0
                    ackl(kc) = 1.d0
                    apkl(kc) = -distx(n)
                    fkl(kc) = distx(n)*vx(kp,jc,ic) - vx(kc,jc,ic)
                    n = n + 1
                elseif (ibmaskx(kc,jc,ic) == -1) then ! Lower boundary points
                    amkl(kc) = -distx(n)
                    ackl(kc) = 1.d0
                    apkl(kc) = 0.d0
                    fkl(kc) = distx(n)*vx(km,jc,ic) - vx(kc,jc,ic)
                    n = n + 1
                end if
            end do
            amkl(1) = 0.d0
            ackl(1) = 1.d0
            apkl(1) = 0.d0
            fkl(1 )= -vx(1, jc,ic)
            fkl(nx)= -vx(nx,jc,ic)
         
            amkT=amkl(2:nx)
            apkT=apkl(1:(nx-1))
            ackT=ackl(1:nx) 
            call dgttrf(nx,amkT,ackT,apkT,appk,ipkv,info)
            call dgttrs('N',nx,1,amkT,ackT,apkT,appk,ipkv,fkl,nx,info)

            do kc=2,nxm
                vx(kc,jc,ic)=vx(kc,jc,ic) + fkl(kc)
            end do
        end do
    end do

    return
end subroutine SolveImpEqnUpdate_X_ibm