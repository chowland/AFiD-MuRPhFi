subroutine InterpInputPhi

    use param
    use input_grids
    ! use mgrd_arrays
    use afid_phasefield, only: phi
    use mpih
    use decomp_2d
    use AuxiliaryRoutines
    use HermiteInterpolations, only: interpolate_xyz_old_to_new_ref
    use afid_sides
    implicit none

    integer  :: ic,jc
    real :: dyo, dzo

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

    !! Add side boundary conditions if using
    if (sidewall) then
        dyo = ycro(2) - ycro(1)
        dzo = zcro(2) - zcro(1)
        if (xstarto(2)==1) then
            call ApplyBC(phio, bc_phi_y_fix_lo, bc_phi_y_val_lo, 'y', 'L', 2, dyo)
        end if
        if (xstarto(3)==1) then
            call ApplyBC(phio, bc_phi_z_fix_lo, bc_phi_z_val_lo, 'z', 'L', 2, dzo)
        end if
        if (xendo(2)==nymro) then
            call ApplyBC(phio, bc_phi_y_fix_up, bc_phi_y_val_up, 'y', 'U', 2, dyo)
        end if
        if (xendo(3)==nzmro) then
            call ApplyBC(phio, bc_phi_z_fix_up, bc_phi_z_val_up, 'z', 'U', 2, dzo)
        end if
    end if

    call interpolate_xyz_old_to_new_ref(phio, phi(1:nxmr,:,:))

    call DestroyReal3DArray(phio)

    return

end subroutine InterpInputPhi