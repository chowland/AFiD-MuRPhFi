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
    real :: dxb, dxbr, delT, a, b, c1, c2, aa, bb, cc, m, S0

    ! Interpolate (T(1,:,:) to refined grid in y and z)
    do ic=xstart(3)-1,xend(3)
        do jc=xstart(2)-1,xend(2)
            qv2(:,:) = temp(1,jc-1:jc+2,ic-1:ic+2)
            do icr=max(krangs(ic),1),min(krangs(ic+1)-1,nzmr)
                qv1(:) = qv2(:,1)*czrs(1,icr) + qv2(:,2)*czrs(2,icr) &
                        +qv2(:,3)*czrs(3,icr) + qv2(:,4)*czrs(4,icr)
                do jcr=max(jrangs(jc),1),min(jrangs(jc+1)-1,nymr)
                    tempr(1,jcr,icr) = sum(qv1(1:4)*cyrs(1:4,jcr))
                end do
            end do
        end do
    end do

    ! Distance to boundary
    dxb = xm(1) - xc(1)
    dxbr = xmr(1) - xcr(1)
    ! Calculate coefficients for equations
    ! dT/dx = a * m
    ! dS/dx = b * S * m
    ! T = c_1 * (S+S0) + c_2
    ! and update boundary values
    ! delT = 35.0*7.86e-4/3.87e-5*Rrho
    a = 25.2    !920.0/1025.0*3974.0/3.35e5/delT*pect
    b = 160.3    !920.0/1025.0*pecs*25 ! 25 factor to get Sc=2500
    c1 = -0.38      !-5.73e-2 * 35.0 /delT
    c2 = -1.28      !c1 + (8.32e-2 + 7.61e-2)/delT - 1.0
    S0 = 1.75
    do icr=xstartr(3),xendr(3)
        do jcr=xstartr(2),xendr(2)
            aa = dxb * dxbr * a * b
            bb = dxb*a + dxbr*b*c2 - dxbr*b*tempr(1,jcr,icr) - c1*dxbr*b*S0
            cc = c1*sal(1,jcr,icr) + c2 - tempr(1,jcr,icr)
            m = (-bb + sqrt(bb**2 - 4*aa*cc))/2/aa
            salbp(1,jcr,icr) = (sal(1,jcr,icr)-dxbr*b*m*S0)/(1+dxbr*b*m)
            saltp(1,jcr,icr) = 0.d0
            tempr(1,jcr,icr) =  tempr(1,jcr,icr) - dxb*a*m
        end do
    end do

    !CJH Update halos for accurate interpolation
    call update_halo(tempr,lvlhalo)
    call update_halo(salbp,lvlhalo)
    call update_halo(saltp,lvlhalo)

    ! Interpolate the temperature BC back to the coarse field
    do ic=xstart(3),xend(3)
        icr = krangs(ic) - 1
        do jc=xstart(2),xend(2)
            jcr = jrangs(jc) - 1
            qv2(:,:) = tempr(1,jcr-1:jcr+2,icr-1:icr+2)
            qv1(:) = qv2(:,1)*czsalc(1,ic) + qv2(:,2)*czsalc(2,ic) &
                    +qv2(:,3)*czsalc(3,ic) + qv2(:,4)*czsalc(4,ic)
            tempbp(1,jc,ic) = sum(qv1(1:4) * cysalc(1:4,jc))
            temptp(1,jc,ic) = 0.d0
        end do
    end do

    return
end subroutine UpdateBCs