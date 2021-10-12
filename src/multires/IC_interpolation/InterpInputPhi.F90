subroutine InterpInputPhi

    use param
    use input_grids
    use mgrd_arrays
    use mpih
    use decomp_2d
    use AuxiliaryRoutines
    use HermiteInterpolations, only: interpolate_xyz_old_to_new_ref
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

    call interpolate_xyz_old_to_new_ref(phio, phi(1:nxmr,:,:))

    call DestroyReal3DArray(phio)

    return

end subroutine InterpInputPhi