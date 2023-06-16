!> Impose the boundary conditions at the edge of the domain in y
!! if we are using the discrete cosine transform for the pressure
!! solution. By default, we are imposing zero velocity (incl. no slip)
!! and zero flux for scalar fields. In z the default is to
!! impose free-slip walls so that 2D simulations are possible
!! N.B. There is a lot of repeated code here that could probably benefit
!! from an ApplyBC function with arbitrary input argument...
subroutine SetSidewallBCs
    use local_arrays, only: vx, vy, vz, temp
    use decomp_2d, only: xstart, xend
    use param
    use afid_salinity, only: sal
    use afid_phasefield, only: phi
    use afid_moisture, only: humid
    implicit none
    real :: dyy, dzz    !! Grid spacing

    dyy = yc(2)
    dzz = zc(2)

    !! Left wall
    if (xstart(2)==1) then
        if (bc_vx_y_fix_lo) then
            vx(:,0,:) = 2.0*bc_vx_y_val_lo - vx(:,1,:)
        else
            vx(:,0,:) = vx(:,1,:) - dyy*bc_vx_y_val_lo
        end if
        !! Always impose zero normal velocity at walls
        vy(:,1,:) = 0.d0
        if (bc_vz_y_fix_lo) then
            vz(:,0,:) = 2.0*bc_vz_y_val_lo - vz(:,1,:)
        else
            vz(:,0,:) = vz(:,1,:) - dyy*bc_vz_y_val_lo
        end if
        if (bc_temp_y_fix_lo) then
            temp(:,0,:) = 2.0*bc_temp_y_val_lo - temp(:,1,:)
        else
            temp(:,0,:) = temp(:,1,:) - dyy*bc_temp_y_val_lo
        end if
        if (salinity) then
            if (bc_sal_y_fix_lo) then
                sal(:,0,:) = 2.0*bc_sal_y_val_lo - sal(:,1,:)
            else
                sal(:,0,:) = sal(:,1,:) - dyy*bc_sal_y_val_lo
            end if
        end if
        if (phasefield) then
            if (bc_phi_y_fix_lo) then
                phi(:,0,:) = 2.0*bc_phi_y_val_lo - phi(:,1,:)
            else
                phi(:,0,:) = phi(:,1,:) - dyy*bc_phi_y_val_lo
            end if
        end if
        if (moist) then
            if (bc_humid_y_fix_lo) then
                humid(:,0,:) = 2.0*bc_humid_y_val_lo - humid(:,1,:)
            else
                humid(:,0,:) = humid(:,1,:) - dyy*bc_humid_y_val_lo
            end if
        end if
    end if
    if (xstart(3)==1) then
        if (bc_vx_z_fix_lo) then
            vx(:,:,0) = 2.0*bc_vx_z_val_lo - vx(:,:,1)
        else
            vx(:,:,0) = vx(:,:,1) - dzz*bc_vx_z_val_lo
        end if
        if (bc_vy_z_fix_lo) then
            vy(:,:,0) = 2.0*bc_vy_z_val_lo - vy(:,:,1)
        else
            vy(:,:,0) = vy(:,:,1) - dzz*bc_vy_z_val_lo
        end if
        !! Always impose zero normal velocity at walls
        vz(:,:,1) = 0.d0
        if (bc_temp_z_fix_lo) then
            temp(:,:,0) = 2.0*bc_temp_z_val_lo - temp(:,:,1)
        else
            temp(:,:,0) = temp(:,:,1) - dzz*bc_temp_z_val_lo
        end if
        if (salinity) then
            if (bc_sal_z_fix_lo) then
                sal(:,:,0) = 2.0*bc_sal_z_val_lo - sal(:,:,1)
            else
                sal(:,:,0) = sal(:,:,1) - dzz*bc_sal_z_val_lo
            end if
        end if
        if (phasefield) then
            if (bc_phi_z_fix_lo) then
                phi(:,:,0) = 2.0*bc_phi_z_val_lo - phi(:,:,1)
            else
                phi(:,:,0) = phi(:,:,1) - dzz*bc_phi_z_val_lo
            end if
        end if
        if (moist) then
            if (bc_humid_z_fix_lo) then
                humid(:,:,0) = 2.0*bc_humid_z_val_lo - humid(:,:,1)
            else
                humid(:,:,0) = humid(:,:,1) - dzz*bc_humid_z_val_lo
            end if
        end if
    end if

    !! Right wall
    if (xend(2)==nym) then
        if (bc_vx_y_fix_up) then
            vx(:,ny,:) = 2.0*bc_vx_y_val_up - vx(:,nym,:)
        else
            vx(:,ny,:) = vx(:,nym,:) - dyy*bc_vx_y_val_up
        end if
        !! Always impose zero normal velocity at walls
        vy(:,ny,:) = 0.d0
        if (bc_vz_y_fix_up) then
            vz(:,ny,:) = 2.0*bc_vz_y_val_up - vz(:,nym,:)
        else
            vz(:,ny,:) = vz(:,nym,:) - dyy*bc_vz_y_val_up
        end if
        if (bc_temp_y_fix_up) then
            temp(:,ny,:) = 2.0*bc_temp_y_val_up - temp(:,nym,:)
        else
            temp(:,ny,:) = temp(:,nym,:) - dyy*bc_temp_y_val_up
        end if
        if (salinity) then
            if (bc_sal_y_fix_up) then
                sal(:,nyr,:) = 2.0*bc_sal_y_val_up - sal(:,nymr,:)
            else
                sal(:,nyr,:) = sal(:,nymr,:) - dyy*bc_sal_y_val_up
            end if
        end if
        if (phasefield) then
            if (bc_phi_y_fix_up) then
                phi(:,nyr,:) = 2.0*bc_phi_y_val_up - phi(:,nymr,:)
            else
                phi(:,nyr,:) = phi(:,nymr,:) - dyy*bc_phi_y_val_up
            end if
        end if
        if (moist) then
            if (bc_humid_y_fix_up) then
                humid(:,ny,:) = 2.0*bc_humid_y_val_up - humid(:,nym,:)
            else
                humid(:,ny,:) = humid(:,nym,:) - dyy*bc_humid_y_val_up
            end if
        end if
    end if
    if (xend(3)==nzm) then
        if (bc_vx_z_fix_up) then
            vx(:,:,nz) = 2.0*bc_vx_z_val_up - vx(:,:,nzm)
        else
            vx(:,:,nz) = vx(:,:,nzm) - dzz*bc_vx_z_val_up
        end if
        if (bc_vy_z_fix_up) then
            vy(:,:,nz) = 2.0*bc_vy_z_val_up - vy(:,:,nzm)
        else
            vy(:,:,nz) = vy(:,:,nzm) - dzz*bc_vy_z_val_up
        end if
        !! Always impose zero normal velocity at walls
        vz(:,:,nz) = 0.d0
        if (bc_temp_z_fix_up) then
            temp(:,:,nz) = 2.0*bc_temp_z_val_up - temp(:,:,nzm)
        else
            temp(:,:,nz) = temp(:,:,nzm) - dzz*bc_temp_z_val_up
        end if
        if (salinity) then
            if (bc_sal_z_fix_up) then
                sal(:,:,nzr) = 2.0*bc_sal_z_val_up - sal(:,:,nzmr)
            else
                sal(:,:,nzr) = sal(:,:,nzmr) - dzz*bc_sal_z_val_up
            end if
        end if
        if (phasefield) then
            if (bc_phi_z_fix_up) then
                phi(:,:,nzr) = 2.0*bc_phi_z_val_up - phi(:,:,nzmr)
            else
                phi(:,:,nzr) = phi(:,:,nzmr) - dzz*bc_phi_z_val_up
            end if
        end if
        if (moist) then
            if (bc_humid_z_fix_up) then
                humid(:,:,nz) = 2.0*bc_humid_z_val_up - humid(:,:,nzm)
            else
                humid(:,:,nz) = humid(:,:,nzm) - dzz*bc_humid_z_val_up
            end if
        end if
    end if
end subroutine SetSidewallBCs