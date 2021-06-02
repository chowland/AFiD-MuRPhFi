!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: ExplicitTermsPhi.F90                           !
!    CONTAINS: subroutine ExplicitTermsPhi                !
!                                                         ! 
!    PURPOSE: Compute the non-linear terms associated to  !
!     the phase-field variable.                           !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ExplicitTermsPhi
    use param
    use mgrd_arrays, only: phi,hphi,tempr
    use decomp_2d, only: xstartr,xendr
    implicit none
    integer :: jc,kc,ic
    integer :: km,kp,jm,jp,im,ip
    real    :: nlphi
    real    :: udzrq,udyrq
    real    :: dyyp,dzzp
    real    :: pf_B

    pf_B = pf_A / (pf_eps)**2

    udzrq=dzqr
    udyrq=dyqr

    do ic=xstartr(3),xendr(3)
        im=ic-1
        ip=ic+1
        do jc=xstartr(2),xendr(2)
            jm=jc-1
            jp=jc+1
            do kc=1,nxmr
                km=kc-1
                kp=kc+1
                ! yy second derivative of phi
                dyyp = (phi(kc,jp,ic) - 2.0*phi(kc,jc,ic) + phi(kc,jm,ic))*udyrq
                ! zz second derivative of phi
                dzzp = (phi(kc,jc,ip) - 2.0*phi(kc,jc,ic) + phi(kc,jc,im))*udzrq
                ! Extra nonlinear terms
                nlphi = pf_B*phi(kc,jc,ic)*(1.0 - phi(kc,jc,ic)) &
                        *(1.0 - 2.0*phi(kc,jc,ic) + pf_C*(tempr(kc,jc,ic) - 0.2))

                hphi(kc,jc,ic) = pf_A*(dyyp + dzzp) - nlphi
            end do
        end do
    end do

    return
end subroutine ExplicitTermsPhi