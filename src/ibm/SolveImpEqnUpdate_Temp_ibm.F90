!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: SolveImpEqnUpdate_Temp.F90                     !
!    CONTAINS: subroutine SolveImpEqnUpdate_Temp          !
!                                                         ! 
!    PURPOSE: Inverts the implicit equation for           !
!     temperature, and updates it to time t+dt            !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine SolveImpEqnUpdate_Temp_ibm
    use param
    use local_arrays, only : temp,rhs,hro
    use decomp_2d, only: xstart,xend
    use ibm_param, only: forclo
    ! use param_particle      ! SL
    implicit none
    real, dimension(nx) :: amkl,apkl,ackl,fkl
    integer :: jc,kc,info,ipkv(nxm),ic,nrhs,km,kp
    real :: betadx,ackl_b
    real :: amkT(nxm-1),ackT(nxm),apkT(nxm-1),appk(nxm-2)
    ! real :: kcpp,kcpm       ! SL

!   Calculate the coefficients of the tridiagonal matrix
!   The coefficients are normalized to prevent floating
!   point errors.

    betadx=0.5d0*al*dt/pect

!   Call to LAPACK library to factor tridiagonal matrix.
!   No solving is done in this call.

    do ic=xstart(3),xend(3)
        do jc=xstart(2),xend(2)

!   Normalize RHS of equation

            do kc=1,nxm
! SL ==================================================   
                ! if(ifparticle.eq.0) then
                km = kc - 1
                kp = kc + 1
                ackl_b=1.0d0/(1.-ac3ssk(kc)*forclo(kc,jc,ic)*betadx)
                amkl(kc)=-am3ssk(kc)*forclo(kc,jc,ic)*betadx*ackl_b &
                        - (1.0 - forclo(kc,jc,ic))*forclo(km,jc,ic)*hro(kc,jc,ic)
                ackl(kc)=1.d0
                apkl(kc)=-ap3ssk(kc)*forclo(kc,jc,ic)*betadx*ackl_b &
                        - (1.0 - forclo(kc,jc,ic))*forclo(kp,jc,ic)*hro(kc,jc,ic)
                ! end if
        
                ! if(ifparticle.ne.0) then    ! SL  consider kcp, rhocpcp
                    
                !     kcpp = (kcp(kc+1,jc,ic)+kcp(kc,jc,ic))/2.0
                !     kcpm = (kcp(kc,jc,ic)+kcp(kc-1,jc,ic))/2.0

                !     ackl_b=1.0d0/(1.-ac3ssk(kc)*forclo(kc,jc,ic)*betadx*(kcpp+kcpm)/2.0/rhocpcp(kc,jc,ic))
                !     amkl(kc)=-am3ssk(kc)*forclo(kc,jc,ic)*betadx*ackl_b*kcpm/rhocpcp(kc,jc,ic)
                !     ackl(kc)=1.d0
                !     apkl(kc)=-ap3ssk(kc)*forclo(kc,jc,ic)*betadx*ackl_b*kcpp/rhocpcp(kc,jc,ic)
                ! end if
! ======================================================
                fkl(kc)=rhs(kc,jc,ic)*ackl_b
            end do
            fkl(nx)= 0.d0
            amkT=amkl(2:nxm)
            apkT=apkl(1:(nxm-1))
            ackT=ackl(1:nxm) 
!     Solve equation using LAPACK library
            call dgttrf(nxm,amkT,ackT,apkT,appk,ipkv,info)
            call dgttrs('N',nxm,1,amkT,ackT,apkT,appk,ipkv,fkl,nx,info)
          
!      Update temperature field

            do kc=1,nxm
                temp(kc,jc,ic) = temp(kc,jc,ic) + fkl(kc)
            end do

        end do
    end do

    return
end subroutine SolveImpEqnUpdate_Temp_ibm