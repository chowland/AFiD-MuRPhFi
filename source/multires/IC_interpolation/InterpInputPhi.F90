subroutine InterpInputPhi

    use param
    use input_grids
    use mgrd_arrays
    use mpih
    use decomp_2d
    use AuxiliaryRoutines
    implicit none

    integer  :: ic,jc,kc, ip,jp,kp, icr,jcr,kcr

    real,dimension(4,4,4) :: qv3 
    real,dimension(4,4) :: qv2
    real,dimension(4) :: qv1

    real, allocatable, dimension(:,:,:) :: phio

    phi(:,:,:) = 0.d0

!=========================================================
!     Interpolation of phi

    ! Allocate and read in old phi
    call AllocateReal3DArray(phio, -1, nxro+1, xstarto(2)-lvlhalo, xendo(2)+lvlhalo, xstarto(3)-lvlhalo, xendo(3)+lvlhalo)

    phio(:,:,:) = 0.d0

    call HdfReadContinua(nzro, nyro, nxro, xstarto(2), xendo(2), xstarto(3), xendo(3), 6, &
            phio(1:nxro, xstarto(2)-lvlhalo:xendo(2)+lvlhalo, xstarto(3)-lvlhalo:xendo(3)+lvlhalo))

    do ic=xstarto(3),xendo(3) ! BC: dphi/dx=0
        do jc=xstarto(2),xendo(2)
            phio(0,jc,ic) = phio(1,jc,ic)
            phio(nxro,jc,ic) = phio(nxmro,jc,ic)
        end do
    end do

    call update_halo(phio, lvlhalo)

    !-- Interpolate phase-field to refined grid
    do ic=xstarto(3)-1,xendo(3)
        do jc=xstarto(2)-1,xendo(2)
            do kc=0,nxmro

                qv3=phio(kc-1:kc+2,jc-1:jc+2,ic-1:ic+2)

                do icr=max(krangsr(ic),1),min(krangsr(ic+1)-1,nzmr)
                    qv2(:,:) = qv3(:,:,1)*czrs(1,icr)+qv3(:,:,2)*czrs(2,icr) &
                                +qv3(:,:,3)*czrs(3,icr)+qv3(:,:,4)*czrs(4,icr)
                    do jcr=max(jrangsr(jc),1),min(jrangsr(jc+1)-1,nymr)
                        qv1(:) = qv2(:,1)*cyrs(1,jcr)+qv2(:,2)*cyrs(2,jcr) &
                                +qv2(:,3)*cyrs(3,jcr)+qv2(:,4)*cyrs(4,jcr)
                        do kcr=max(irangsr(kc),1),min(irangsr(kc+1)-1,nxmr)
                            phi(kcr,jcr,icr) = sum(qv1*cxrs(:,kcr))
                        end do
                    end do
                end do

            end do
        end do
    end do

    call DestroyReal3DArray(phio)

    return

end subroutine InterpInputPhi