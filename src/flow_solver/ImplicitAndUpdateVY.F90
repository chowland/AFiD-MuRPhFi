!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: ImplicitAndUpdateVY.F90                        !
!    CONTAINS: subroutine ImplicitAndUpdateVY             !
!                                                         ! 
!    PURPOSE: Compute the linear terms associated to      !
!     the velocity in the y (horizontal) dimension        !
!     and call the implicit solver                        !
!     After this routine, the velocity field in y has     !
!     been updated to the new timestep                    !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ImplicitAndUpdateVY
    use param
    use local_arrays, only: vy,ruy,pr,rhs,dph
    use decomp_2d, only: xstart,xend
    use ibm_param
    implicit none
    integer :: kc,jmm,jc,ic
    integer :: kpp,kmm,ke
    real    :: alre,udy
    real    :: amm,acc,app
    real    :: dyp,dxxvy
    real    :: usaldto,q2e


    alre=al/ren
    udy=dy*al

!$OMP  PARALLEL DO &
!$OMP   DEFAULT(none) &
!$OMP   SHARED(xstart,xend,nxm,vy,pr,xminusU,xplusU) &
!$OMP   SHARED(kmv,kpv,am3sk,ac3sk,ap3sk) &
!$OMP   SHARED(dy,al,ga,ro,alre,dt,dph) &
!$OMP   SHARED(udy,udx3m,rhs,ruy) &
!$OMP   PRIVATE(ic,jc,kc,kmm,kpp,jmm) &
!$OMP   PRIVATE(amm,acc,app) &
!$OMP   PRIVATE(dyp,dxxvy)
    do ic=xstart(3),xend(3)
        do jc=xstart(2),xend(2)
            jmm = jc-1
            do kc=1,nxm
                kmm = kmv(kc)
                kpp = kpv(kc)
                amm = am3sk(kc)
                acc = ac3sk(kc)
                app = ap3sk(kc)


!   Second derivative in x-direction of vy
!
!
                if(kc.eq.1) then
                    dxxvy = vy(kpp,jc,ic)*app &
                            + vy(kc,jc,ic)*acc &
                            - (acc + app)*xminusU*inslws
                elseif(kc.eq.nxm) then
                    dxxvy = -(acc + amm)*xplusU*inslwn &
                            + vy(kc,jc,ic)*acc &
                            + vy(kmm,jc,ic)*amm
                else
                    dxxvy = vy(kpp,jc,ic)*app &
                            + vy(kc,jc,ic)*acc &
                            + vy(kmm,jc,ic)*amm
                end if

!   component of grad(pr) along y direction
!
                dyp = (pr(kc,jc,ic) - pr(kc,jmm,ic))*udy

!    Calculate right hand side of Eq. 5 (VO96)
!
                rhs(kc,jc,ic) = (ga*dph(kc,jc,ic) + ro*ruy(kc,jc,ic) &
                            + alre*dxxvy - dyp)*dt

!    Store the non-linear terms for the calculation of 
!    the next timestep

                ruy(kc,jc,ic)=dph(kc,jc,ic)
            enddo
        enddo
    enddo

!  Solve equation and update velocity

    if (IBM) then
        forclo = 1.d0
        usaldto = 1.0/aldto
        do n=1,npuny
            ic = indgeo(2,n,1)
            jc = indgeo(2,n,2)
            kc = indgeo(2,n,3)
            forclo(kc,jc,ic) = 0.d0
            ke = indgeoe(2,n,3)
            q2e = ((al*dt + aldto)*vy(ke,jc,ic) - al*dt*q2bo(n))*usaldto
            rhs(kc,jc,ic) = -vy(kc,jc,ic) + q2e*distb(2,n)
            q2bo(n) = vy(ke,jc,ic)
        end do
        call SolveImpEqnUpdate_YZ_ibm(vy,rhs,forclo)
    else
        call SolveImpEqnUpdate_YZ(vy,rhs)
    end if
      
    return
end subroutine ImplicitAndUpdateVY