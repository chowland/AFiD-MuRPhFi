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
    use mgrd_arrays, only: phi,hphi,tempr,sal
    use decomp_2d, only: xstartr,xendr
    implicit none
    integer :: jc,kc,ic
    integer :: jm,jp,im,ip
    real    :: nlphi
    real    :: udzrq,udyrq
    real    :: dyyp,dzzp
    real    :: pf_B, bcl

    pf_B = pf_D / (pf_eps)**2

    udzrq=dzqr
    udyrq=dyqr

    do ic=xstartr(3),xendr(3)
        im=ic-1
        ip=ic+1
        do jc=xstartr(2),xendr(2)
            jm=jc-1
            jp=jc+1
            do kc=1,nxmr
                ! yy second derivative of phi
                dyyp = (phi(kc,jp,ic) - 2.0*phi(kc,jc,ic) + phi(kc,jm,ic))*udyrq
                ! zz second derivative of phi
                dzzp = (phi(kc,jc,ip) - 2.0*phi(kc,jc,ic) + phi(kc,jc,im))*udzrq
                ! Extra nonlinear terms
                nlphi = pf_B*phi(kc,jc,ic)*(1.0 - phi(kc,jc,ic)) &
                        *(1.0 - 2.0*phi(kc,jc,ic) + pf_A*(tempr(kc,jc,ic) - pf_Tm))

                hphi(kc,jc,ic) = pf_D*(dyyp + dzzp) - nlphi
            end do
        end do
    end do

    if (salinity) then
        bcl = pf_B*pf_A*pf_Lambda
        do ic=xstartr(3),xendr(3)
            im=ic-1
            ip=ic+1
            do jc=xstartr(2),xendr(2)
                jm=jc-1
                jp=jc+1
                do kc=1,nxmr
                    hphi(kc,jc,ic) = hphi(kc,jc,ic) - &
                        bcl*sal(kc,jc,ic)*phi(kc,jc,ic)*(1.0 - phi(kc,jc,ic))
                end do
            end do
        end do
    end if

    return
end subroutine ExplicitTermsPhi