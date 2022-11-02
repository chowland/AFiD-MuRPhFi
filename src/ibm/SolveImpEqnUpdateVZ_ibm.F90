!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: SolveImpEqnUpdate_VY.F90                       !
!    CONTAINS: subroutine SolveImpEqnUpdate_VY            !
!                                                         ! 
!    PURPOSE: Inverts the implicit equation for velocity  !
!     in any of the horizontal directions, and updates    !
!     it to time t+dt                                     !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
subroutine SolveImpEqnUpdate_VZ_ibm
    use param
    use local_arrays, only: vz, rhs
    use decomp_2d, only: xstart,xend
    use ibm_param, only: ibmaskz, distz
    implicit none
    real, dimension(nx) :: amkl,apkl,ackl,fkl
    integer :: jc,kc,info,ipkv(nxm),ic,km,kp,n
    real :: betadx,ackl_b
    real :: amkT(nxm-1),ackT(nxm),apkT(nxm-1),appk(nx-3)

    betadx = beta*al
    ackl(:) = 1.0d0

    do ic=xstart(3),xend(3)
        do jc=xstart(2),xend(2)
            do kc=1,nxm
                km = max(1,kc - 1)
                kp = min(kc + 1,nxm)
                if (ibmaskz(kc,jc,ic) == 2) then ! Liquid phase
                    ackl_b = 1.0d0/(1.0d0 - ac3sk(kc)*betadx)
                    amkl(kc) = -am3sk(kc)*betadx*ackl_b
                    ackl(kc) = 1.d0
                    apkl(kc) = -ap3sk(kc)*betadx*ackl_b
                    fkl(kc) = rhs(kc,jc,ic)*ackl_b
                elseif (ibmaskz(kc,jc,ic) == 0) then ! Solid phase
                    amkl(kc) = 0.d0
                    ackl(kc) = 1.d0
                    apkl(kc) = 0.d0
                    fkl(kc) = -vz(kc,jc,ic)
                elseif (ibmaskz(kc,jc,ic) == 1) then ! Upper boundary points
                    amkl(kc) = 0.d0
                    ackl(kc) = 1.d0
                    apkl(kc) = -distz(n)
                    fkl(kc) = distz(n)*vz(kp,jc,ic) - vz(kc,jc,ic)
                    n = n + 1
                elseif (ibmaskz(kc,jc,ic) == -1) then ! Lower boundary points
                    amkl(kc) = -distz(n)
                    ackl(kc) = 1.d0
                    apkl(kc) = 0.d0
                    fkl(kc) = distz(n)*vz(km,jc,ic) - vz(kc,jc,ic)
                    n = n + 1
                end if
            end do
            amkT = amkl(2:nxm)
            apkT = apkl(1:(nxm-1))
            ackT = ackl(1:nxm)

            call dgttrf(nxm,amkT,ackT,apkT,appk,ipkv,info)
            call dgttrs('N',nxm,1,amkT,ackT,apkT,appk,ipkv,fkl,nxm,info)

            do kc=1,nxm
                vz(kc,jc,ic) = vz(kc,jc,ic) + fkl(kc)
            end do
        end do
    end do

    return
end subroutine SolveImpEqnUpdate_VZ_ibm