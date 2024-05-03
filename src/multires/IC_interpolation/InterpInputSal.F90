subroutine InterpInputSal
    
    use param
    use input_grids
    use afid_salinity
    use mpih
    use decomp_2d
    use AuxiliaryRoutines
    use HermiteInterpolations, only: interpolate_xyz_old_to_new_ref
    use afid_sides
    implicit none

    integer  :: ic,jc

    real :: Sup, Slo
    real :: dyo, dzo

    real, allocatable, dimension(:,:,:) :: salo
    
    sal(:,:,:) = 0.d0

!=========================================================
!     Interpolation of S

    ! Allocate and read in old S
    call AllocateReal3DArray(salo, -1, nxro+1, xstarto(2)-lvlhalo, xendo(2)+lvlhalo, xstarto(3)-lvlhalo, xendo(3)+lvlhalo)

    salo(:,:,:) = 0.d0

    call HdfReadContinua(nzro, nyro, nxro, xstarto(2), xendo(2), xstarto(3), xendo(3), 5, &
            salo(1:nxro, xstarto(2)-lvlhalo:xendo(2)+lvlhalo, xstarto(3)-lvlhalo:xendo(3)+lvlhalo))

    ! Salinity BCs
    if (RAYS>0) then
        Sup = 0.5
        Slo = -0.5
    else
        Sup = -0.5
        Slo = +0.5
    end if
    if (phasefield) then
        Sup = 0.0
        Slo = 1.0
    end if
        
    !-- Boundary points
    do ic=xstarto(3),xendo(3)
        do jc=xstarto(2),xendo(2)
            if (SfixS==1) then
                salo(0,jc,ic) = 2.0*Slo - salo(1,jc,ic)
            else
                salo(0,jc,ic) = salo(1,jc,ic)
            end if
            if (SfixN==1) then
                salo(nxro,jc,ic) = 2.0*Sup - salo(nxmro,jc,ic)
            else
                salo(nxro,jc,ic) = salo(nxmro,jc,ic)
            end if
        end do
    end do

    call update_halo(salo, lvlhalo)

    !! Add side boundary conditions if using
    if (sidewall) then
        dyo = ycro(2) - ycro(1)
        dzo = zcro(2) - zcro(1)
        if (xstarto(2)==1) then
            call ApplyBC(salo, bc_sal_y_fix_lo, bc_sal_y_val_lo, 'y', 'L', 2, dyo)
        end if
        if (xstarto(3)==1) then
            call ApplyBC(salo, bc_sal_z_fix_lo, bc_sal_z_val_lo, 'z', 'L', 2, dzo)
        end if
        if (xendo(2)==nymro) then
            call ApplyBC(salo, bc_sal_y_fix_up, bc_sal_y_val_up, 'y', 'U', 2, dyo)
        end if
        if (xendo(3)==nzmro) then
            call ApplyBC(salo, bc_sal_z_fix_up, bc_sal_z_val_up, 'z', 'U', 2, dzo)
        end if
    end if

    call interpolate_xyz_old_to_new_ref(salo, sal(1:nxmr,:,:))

    call DestroyReal3DArray(salo)

    return

end subroutine InterpInputSal