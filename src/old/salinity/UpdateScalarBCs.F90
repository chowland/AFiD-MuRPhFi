!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: UpdateScalarBCs.F90                            !
!    CONTAINS: subroutine UpdateBCs                       !
!                                                         ! 
!    PURPOSE: Updates the temperature and salinity        !
!     boundary conditions to simulate a melting wall      !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine UpdateBCs
    use param
    use local_arrays, only: temp
    use mgrd_arrays
    use decomp_2d
    implicit none
    integer :: ic, jc, icr, jcr
    real, dimension(4,4) :: qv2
    real, dimension(4) :: qv1
    real :: dxb, dxbr, Sb, aa, bb, cc!, m

    ! Interpolate (T(1,:,:) to refined grid in y and z)
    do ic=xstart(3)-1,xend(3)
        do jc=xstart(2)-1,xend(2)
            qv2(:,:) = temp(1,jc-1:jc+2,ic-1:ic+2)
            do icr=max(krangs(ic),xstartr(3)),min(krangs(ic+1)-1,xendr(3))
                qv1(:) = qv2(:,1)*czrs(1,icr) + qv2(:,2)*czrs(2,icr) &
                        +qv2(:,3)*czrs(3,icr) + qv2(:,4)*czrs(4,icr)
                do jcr=max(jrangs(jc),xstartr(2)),min(jrangs(jc+1)-1,xendr(2))
                    Tplaner(1,jcr,icr) = sum(qv1(1:4)*cyrs(1:4,jcr))
                end do
            end do
        end do
    end do

    ! Distance to boundary
    dxb = xm(1) - xc(1)
    dxbr = xmr(1) - xcr(1)
    ! Calculate coefficients for equations
    ! (PeT)^-1 * dT/dx = Stefan * m
    ! (PeS)^-1 * dS/dx = S * m
    ! T + Lambda * S = 0
    ! and update boundary values
    do icr=xstartr(3),xendr(3)
        do jcr=xstartr(2),xendr(2)
            aa = dxbr * pecs * pf_Lambda
            bb = dxbr*pecs*Tplaner(1,jcr,icr) + dxb*pect*pf_S
            cc = -pf_S*dxb*pect*sal(1,jcr,icr)
            Sb = (-bb + sqrt(bb**2 - 4*aa*cc))/2/aa
            salbp(1,jcr,icr) = Sb
            saltp(1,jcr,icr) = 0.d0
            Tplaner(1,jcr,icr) =  -pf_Lambda*Sb
            ! m = (sal(1,jcr,icr) - Sb)/dxbr/pecs/Sb
        end do
    end do

    !CJH Update halos for accurate interpolation
    call update_halo(Tplaner,lvlhalo)
    call update_halo(salbp,lvlhalo)
    call update_halo(saltp,lvlhalo)

    ! Interpolate the temperature BC back to the coarse field
    do icr=xstartr(3)-1,xendr(3)
        do jcr=xstartr(2)-1,xendr(2)
            
            qv2 = Tplaner(1,jcr-1:jcr+2,icr-1:icr+2)

            do ic=max(krangr(icr),xstart(3)),min(krangr(icr+1)-1,xend(3))
                qv1(:) = qv2(:,1)*czsalc(1,ic) + qv2(:,2)*czsalc(2,ic) &
                       + qv2(:,3)*czsalc(3,ic) + qv2(:,4)*czsalc(4,ic)
                do jc=max(jrangr(jcr),xstart(2)),min(jrangr(jcr+1)-1,xend(2))
                    tempbp(1,jc,ic) = sum(qv1(1:4) * cysalc(1:4,jc))
                    temptp(1,jc,ic) = 0.d0
                end do
            end do
        end do
    end do

    return
end subroutine UpdateBCs