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
    real :: dyyr, dzzr    !! Refined grid spacing
    integer :: n

    dyy = yc(2)
    dzz = zc(2)

    if (multires) then
        dyyr = ycr(2)
        dzzr = zcr(2)
    end if

    !! Left wall
    if (xstart(2)==1) then
        if (bc_vx_y_fix_lo) then
            do n=1,lvlhalo
                vx(:,1-n,:) = 2.0*bc_vx_y_val_lo - vx(:,n,:)
            end do
        else
            do n=1,lvlhalo
                vx(:,1-n,:) = vx(:,n,:) - dyy*bc_vx_y_val_lo
            end do
        end if
        !! Always impose zero normal velocity at walls
        vy(:,1,:) = 0.d0
        vy(:,0,:) = vy(:,2,:)
        if (bc_vz_y_fix_lo) then
            do n=1,lvlhalo
                vz(:,1-n,:) = 2.0*bc_vz_y_val_lo - vz(:,n,:)
            end do
        else
            do n=1,lvlhalo
                vz(:,1-n,:) = vz(:,n,:) - dyy*bc_vz_y_val_lo
            end do
        end if
        if (bc_temp_y_fix_lo) then
            do n=1,lvlhalo
                temp(:,1-n,:) = 2.0*bc_temp_y_val_lo - temp(:,n,:)
            end do
        else
            do n=1,lvlhalo
                temp(:,1-n,:) = temp(:,n,:) - dyy*bc_temp_y_val_lo
            end do
        end if
        if (salinity) then
            if (bc_sal_y_fix_lo) then
                do n=1,lvlhalo
                    sal(:,1-n,:) = 2.0*bc_sal_y_val_lo - sal(:,n,:)
                end do
            else
                do n=1,lvlhalo
                    sal(:,1-n,:) = sal(:,n,:) - dyyr*bc_sal_y_val_lo
                end do
            end if
        end if
        if (phasefield) then
            if (bc_phi_y_fix_lo) then
                do n=1,lvlhalo
                    phi(:,1-n,:) = 2.0*bc_phi_y_val_lo - phi(:,n,:)
                end do
            else
                do n=1,lvlhalo
                    phi(:,1-n,:) = phi(:,n,:) - dyyr*bc_phi_y_val_lo
                end do
            end if
        end if
        if (moist) then
            if (bc_humid_y_fix_lo) then
                do n=1,lvlhalo
                    humid(:,1-n,:) = 2.0*bc_humid_y_val_lo - humid(:,n,:)
                end do
            else
                do n=1,lvlhalo
                    humid(:,1-n,:) = humid(:,n,:) - dyy*bc_humid_y_val_lo
                end do
            end if
        end if
    end if
    if (xstart(3)==1) then
        if (  .not. periodic_bc(3)) then
        if (bc_vx_z_fix_lo) then
            do n=1,lvlhalo
                vx(:,:,1-n) = 2.0*bc_vx_z_val_lo - vx(:,:,n)
            end do
        else
            do n=1,lvlhalo
                vx(:,:,1-n) = vx(:,:,n) - dzz*bc_vx_z_val_lo
            end do
        end if
        if (bc_vy_z_fix_lo) then
            do n=1,lvlhalo
                vy(:,:,1-n) = 2.0*bc_vy_z_val_lo - vy(:,:,n)
            end do
        else
            do n=1,lvlhalo
                vy(:,:,1-n) = vy(:,:,n) - dzz*bc_vy_z_val_lo
            end do
        end if
        !! Always impose zero normal velocity at walls
        vz(:,:,1) = 0.d0
        vz(:,:,0) = vz(:,:,2)
        if (bc_temp_z_fix_lo) then
            do n=1,lvlhalo
                temp(:,:,1-n) = 2.0*bc_temp_z_val_lo - temp(:,:,n)
            end do
        else
            do n=1,lvlhalo
                temp(:,:,1-n) = temp(:,:,n) - dzz*bc_temp_z_val_lo
            end do
        end if
        if (salinity) then
            if (bc_sal_z_fix_lo) then
                do n=1,lvlhalo
                    sal(:,:,1-n) = 2.0*bc_sal_z_val_lo - sal(:,:,n)
                end do
            else
                do n=1,lvlhalo
                    sal(:,:,1-n) = sal(:,:,n) - dzzr*bc_sal_z_val_lo
                end do
            end if
        end if
        if (phasefield) then
            if (bc_phi_z_fix_lo) then
                do n=1,lvlhalo
                    phi(:,:,1-n) = 2.0*bc_phi_z_val_lo - phi(:,:,n)
                end do
            else
                do n=1,lvlhalo
                    phi(:,:,1-n) = phi(:,:,n) - dzzr*bc_phi_z_val_lo
                end do
            end if
        end if
        if (moist) then
            if (bc_humid_z_fix_lo) then
                do n=1,lvlhalo
                    humid(:,:,1-n) = 2.0*bc_humid_z_val_lo - humid(:,:,n)
                end do
            else
                do n=1,lvlhalo
                    humid(:,:,1-n) = humid(:,:,n) - dzz*bc_humid_z_val_lo
                end do
            end if
        end if
    end if
    end if
    !! Right wall
    if (xend(2)==nym) then
        if (bc_vx_y_fix_up) then
            do n=1,lvlhalo
                vx(:,nym+n,:) = 2.0*bc_vx_y_val_up - vx(:,nym+1-n,:)
            end do
        else
            do n=1,lvlhalo
                vx(:,nym+n,:) = vx(:,nym+1-n,:) - dyy*bc_vx_y_val_up
            end do
        end if
        !! Always impose zero normal velocity at walls
        vy(:,ny,:) = 0.d0
        if (bc_vz_y_fix_up) then
            do n=1,lvlhalo
                vz(:,nym+n,:) = 2.0*bc_vz_y_val_up - vz(:,nym+1-n,:)
            end do
        else
            do n=1,lvlhalo
                vz(:,nym+n,:) = vz(:,nym+1-n,:) - dyy*bc_vz_y_val_up
            end do
        end if
        if (bc_temp_y_fix_up) then
            do n=1,lvlhalo
                temp(:,nym+n,:) = 2.0*bc_temp_y_val_up - temp(:,nym+1-n,:)
            end do
        else
            do n=1,lvlhalo
                temp(:,nym+n,:) = temp(:,nym+1-n,:) - dyy*bc_temp_y_val_up
            end do
        end if
        if (salinity) then
            if (bc_sal_y_fix_up) then
                do n=1,lvlhalo
                    sal(:,nymr+n,:) = 2.0*bc_sal_y_val_up - sal(:,nymr+1-n,:)
                end do
            else
                do n=1,lvlhalo
                    sal(:,nymr+n,:) = sal(:,nymr+1-n,:) - dyyr*bc_sal_y_val_up
                end do
            end if
        end if
        if (phasefield) then
            if (bc_phi_y_fix_up) then
                do n=1,lvlhalo
                    phi(:,nymr+n,:) = 2.0*bc_phi_y_val_up - phi(:,nymr+1-n,:)
                end do
            else
                do n=1,lvlhalo
                    phi(:,nymr+n,:) = phi(:,nymr+1-n,:) - dyyr*bc_phi_y_val_up
                end do
            end if
        end if
        if (moist) then
            if (bc_humid_y_fix_up) then
                do n=1,lvlhalo
                    humid(:,nym+n,:) = 2.0*bc_humid_y_val_up - humid(:,nym+1-n,:)
                end do
            else
                do n=1,lvlhalo
                    humid(:,nym+n,:) = humid(:,nym+1-n,:) - dyy*bc_humid_y_val_up
                end do
            end if
        end if
    end if
    if (xend(3)==nzm) then
        if ( .not. periodic_bc(3)) then
        if (bc_vx_z_fix_up) then
            do n=1,lvlhalo
                vx(:,:,nzm+n) = 2.0*bc_vx_z_val_up - vx(:,:,nzm+1-n)
            end do
        else
            do n=1,lvlhalo
                vx(:,:,nzm+n) = vx(:,:,nzm+1-n) - dzz*bc_vx_z_val_up
            end do
        end if
        if (bc_vy_z_fix_up) then
            do n=1,lvlhalo
                vy(:,:,nzm+n) = 2.0*bc_vy_z_val_up - vy(:,:,nzm+1-n)
            end do
        else
            do n=1,lvlhalo
                vy(:,:,nzm+n) = vy(:,:,nzm+1-n) - dzz*bc_vy_z_val_up
            end do
        end if
        !! Always impose zero normal velocity at walls
        vz(:,:,nz) = 0.d0
        if (bc_temp_z_fix_up) then
            do n=1,lvlhalo
                temp(:,:,nzm+n) = 2.0*bc_temp_z_val_up - temp(:,:,nzm+1-n)
            end do
        else
            do n=1,lvlhalo
                temp(:,:,nzm+n) = temp(:,:,nzm+1-n) - dzz*bc_temp_z_val_up
            end do
        end if
        if (salinity) then
            if (bc_sal_z_fix_up) then
                do n=1,lvlhalo
                    sal(:,:,nzmr+n) = 2.0*bc_sal_z_val_up - sal(:,:,nzmr+1-n)
                end do
            else
                do n=1,lvlhalo
                    sal(:,:,nzmr+n) = sal(:,:,nzmr+1-n) - dzzr*bc_sal_z_val_up
                end do
            end if
        end if
        if (phasefield) then
            if (bc_phi_z_fix_up) then
                do n=1,lvlhalo
                    phi(:,:,nzmr+n) = 2.0*bc_phi_z_val_up - phi(:,:,nzmr+1-n)
                end do
            else
                do n=1,lvlhalo
                    phi(:,:,nzmr+n) = phi(:,:,nzmr+1-n) - dzzr*bc_phi_z_val_up
                end do
            end if
        end if
        if (moist) then
            if (bc_humid_z_fix_up) then
                do n=1,lvlhalo
                    humid(:,:,nzm+n) = 2.0*bc_humid_z_val_up - humid(:,:,nzm+1-n)
                end do
            else
                do n=1,lvlhalo
                    humid(:,:,nzm+n) = humid(:,:,nzm+1-n) - dzz*bc_humid_z_val_up
                end do
            end if
        end if
    end if
end if
end subroutine SetSidewallBCs