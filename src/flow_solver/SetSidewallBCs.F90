!> Impose the boundary conditions at the edge of the domain in y
!! if we are using the discrete cosine transform for the pressure
!! solution. By default, we are imposing zero velocity (incl. no slip)
!! and zero flux for temperature. Maybe later we can add a flag 
!! to make this more flexible.
!! We impose free-slip walls in z so that 2D simulations are possible
subroutine SetSidewallBCs
    use local_arrays, only: vx, vy, vz, temp
    use decomp_2d, only: xstart, xend
    use param
    use afid_salinity, only: sal
    use afid_phasefield, only: phi
    use afid_moisture, only: humid
    implicit none

    !! Left wall
    if (xstart(2)==1) then
        vx(:,0,:) = -vx(:,1,:)
        vy(:,1,:) = 0.d0
        vz(:,0,:) = -vz(:,1,:)
        temp(:,0,:) = temp(:,1,:)
        if (salinity) sal(:,0,:) = sal(:,1,:)
        if (phasefield) phi(:,0,:) = phi(:,1,:)
        if (moist) humid(:,0,:) = humid(:,1,:)
    end if
    if (xstart(3)==1) then
        vx(:,:,0) = vx(:,:,1)
        vy(:,:,0) = vy(:,:,1)
        vz(:,:,1) = 0.d0
        temp(:,:,0) = temp(:,:,1)
        if (salinity) sal(:,:,0) = sal(:,:,1)
        if (phasefield) phi(:,:,0) = phi(:,:,1)
        if (moist) humid(:,:,0) = humid(:,:,1)
    end if

    !! Right wall
    if (xend(2)==nym) then
        vx(:,ny,:) = -vx(:,nym,:)
        vy(:,ny,:) = 0.d0
        vz(:,ny,:) = -vz(:,nym,:)
        temp(:,ny,:) = temp(:,nym,:)
        if (salinity) sal(:,nyr,:) = sal(:,nymr,:)
        if (phasefield) phi(:,nyr,:) = phi(:,nymr,:)
        if (moist) humid(:,ny,:) = humid(:,nym,:)
    end if
    if (xend(3)==nzm) then
        vx(:,:,nz) = vx(:,:,nzm)
        vy(:,:,nz) = vy(:,:,nzm)
        vz(:,:,nz) = 0.d0
        temp(:,:,nz) = temp(:,:,nzm)
        if (salinity) sal(:,:,nzr) = sal(:,:,nzmr)
        if (phasefield) phi(:,:,nzr) = phi(:,:,nzmr)
        if (moist) humid(:,:,nz) = humid(:,:,nzm)
    end if
end subroutine SetSidewallBCs