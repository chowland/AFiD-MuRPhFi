!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: ImplicitAndUpdateVZ.F90                        !
!    CONTAINS: subroutine ImplicitAndUpdateVZ             !
!                                                         ! 
!    PURPOSE: Compute the linear terms associated to      !
!     the velocity in the z (horizontal) dimension        !
!     and call the implicit solver.                       !
!     After this routine, the velocity field in z has     !
!     been updated to the new timestep                    !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ImplicitAndUpdateVZ
    use param
    use local_arrays, only: vz,dq,ruz,rhs,pr
    use decomp_2d, only: xstart,xend
    use ibm_param
    implicit none
    integer :: kc,jc,ic,imm
    integer :: kmm,kpp,ke
    real    :: alre,amm,acc,app,udz
    real    :: dxxvz,dzp
    real    :: usaldto,q1e

    alre=al/ren
    udz=dz*al

!$OMP  PARALLEL DO &
!$OMP   DEFAULT(none) &
!$OMP   SHARED(xstart,xend,nxm,vz,pr) &
!$OMP   SHARED(kmv,kpv,am3sk,ac3sk,ap3sk) &
!$OMP   SHARED(dz,al,ga,ro,alre,dt,dq) &
!$OMP   SHARED(udx3m,rhs,ruz) &
!$OMP   PRIVATE(ic,jc,kc,imm,kmm,kpp) &
!$OMP   PRIVATE(amm,acc,app) &
!$OMP   PRIVATE(dxxvz,dzp)
    do ic=xstart(3),xend(3)
        imm=ic-1
        do jc=xstart(2),xend(2)
            do kc=1,nxm
                kmm=kmv(kc)
                kpp=kpv(kc)
                amm=am3sk(kc)
                acc=ac3sk(kc)
                app=ap3sk(kc)

!   Second derivative in x-direction of vz
!
                dxxvz=vz(kpp,jc,ic)*app &
                    +vz(kc,jc,ic)*acc &
                    +vz(kmm,jc,ic)*amm
      
!   component of grad(pr) along z direction
!
                dzp=(pr(kc,jc,ic)-pr(kc,jc,imm))*dz*al

!    Calculate right hand side of Eq. 5 (VO96)
!
                rhs(kc,jc,ic)=(ga*dq(kc,jc,ic)+ro*ruz(kc,jc,ic) &
                            +alre*dxxvz-dzp)*dt

!    Store the non-linear terms for the calculation of 
!    the next timestep

                ruz(kc,jc,ic)=dq(kc,jc,ic)
            enddo
        enddo
    enddo
!$OMP END PARALLEL DO

    if (IBM) then
        forclo = 1.d0
        usaldto = 1.0/aldto
        do n=1,npunz
            ic = indgeo(1,n,1)
            jc = indgeo(1,n,2)
            kc = indgeo(1,n,3)
            forclo(kc,jc,ic) = 0.d0
            ke = indgeoe(1,n,3)
            q1e = ((al*dt + aldto)*vz(ke,jc,ic) - al*dt*q1bo(n))*usaldto
            rhs(kc,jc,ic) = -vz(kc,jc,ic) + q1e*distb(1,n)
            q1bo(n) = vz(ke,jc,ic)
        end do
        call SolveImpEqnUpdate_YZ_ibm(vz,rhs,forclo)
    else
        call SolveImpEqnUpdate_YZ(vz,rhs)
    end if

    return
end subroutine ImplicitAndUpdateVZ