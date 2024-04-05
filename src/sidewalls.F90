!> Module containing the flags, parameters and subroutines needed to simulate a box with sidewalls
!! enable this feature using the flag `sidewall` in the param module
module afid_sides

    !! Sidewall boundary conditions
    logical :: bc_vx_y_fix_lo = .true.      !! Dirichlet/Neumann flag for lower y BC for vx
    logical :: bc_vx_y_fix_up = .true.      !! Dirichlet/Neumann flag for upper y BC for vx
    logical :: bc_vx_z_fix_lo = .false.     !! Dirichlet/Neumann flag for lower z BC for vx
    logical :: bc_vx_z_fix_up = .false.     !! Dirichlet/Neumann flag for upper z BC for vx

    logical :: bc_vy_z_fix_lo = .false.     !! Dirichlet/Neumann flag for lower z BC for vy
    logical :: bc_vy_z_fix_up = .false.     !! Dirichlet/Neumann flag for upper z BC for vy

    logical :: bc_vz_y_fix_lo = .true.      !! Dirichlet/Neumann flag for lower y BC for vz
    logical :: bc_vz_y_fix_up = .true.      !! Dirichlet/Neumann flag for upper y BC for vz

    logical :: bc_temp_y_fix_lo = .false.   !! Dirichlet/Neumann flag for lower y BC for temperature
    logical :: bc_temp_y_fix_up = .false.   !! Dirichlet/Neumann flag for upper y BC for temperature
    logical :: bc_temp_z_fix_lo = .false.   !! Dirichlet/Neumann flag for lower z BC for temperature
    logical :: bc_temp_z_fix_up = .false.   !! Dirichlet/Neumann flag for upper z BC for temperature

    logical :: bc_sal_y_fix_lo = .false.    !! Dirichlet/Neumann flag for lower y BC for salinity
    logical :: bc_sal_y_fix_up = .false.    !! Dirichlet/Neumann flag for upper y BC for salinity
    logical :: bc_sal_z_fix_lo = .false.    !! Dirichlet/Neumann flag for lower z BC for salinity
    logical :: bc_sal_z_fix_up = .false.    !! Dirichlet/Neumann flag for upper z BC for salinity

    logical :: bc_phi_y_fix_lo = .false.    !! Dirichlet/Neumann flag for lower y BC for phase-field
    logical :: bc_phi_y_fix_up = .false.    !! Dirichlet/Neumann flag for upper y BC for phase-field
    logical :: bc_phi_z_fix_lo = .false.    !! Dirichlet/Neumann flag for lower z BC for phase-field
    logical :: bc_phi_z_fix_up = .false.    !! Dirichlet/Neumann flag for upper z BC for phase-field

    logical :: bc_humid_y_fix_lo = .false.  !! Dirichlet/Neumann flag for lower y BC for humidity
    logical :: bc_humid_y_fix_up = .false.  !! Dirichlet/Neumann flag for upper y BC for humidity
    logical :: bc_humid_z_fix_lo = .false.  !! Dirichlet/Neumann flag for lower z BC for humidity
    logical :: bc_humid_z_fix_up = .false.  !! Dirichlet/Neumann flag for upper z BC for humidity

    real :: bc_vx_y_val_lo = 0.0        !! Boundary (flux?) value for lower y BC for vx
    real :: bc_vx_y_val_up = 0.0        !! Boundary (flux?) value for upper y BC for vx
    real :: bc_vx_z_val_lo = 0.0        !! Boundary (flux?) value for lower z BC for vx
    real :: bc_vx_z_val_up = 0.0        !! Boundary (flux?) value for upper z BC for vx

    real :: bc_vy_z_val_lo = 0.0        !! Boundary (flux?) value for lower z BC for vy
    real :: bc_vy_z_val_up = 0.0        !! Boundary (flux?) value for upper z BC for vy

    real :: bc_vz_y_val_lo = 0.0        !! Boundary (flux?) value for lower y BC for vz
    real :: bc_vz_y_val_up = 0.0        !! Boundary (flux?) value for upper y BC for vz

    real :: bc_temp_y_val_lo = 0.0      !! Boundary (flux?) value for lower y BC for temperature
    real :: bc_temp_y_val_up = 0.0      !! Boundary (flux?) value for upper y BC for temperature
    real :: bc_temp_z_val_lo = 0.0      !! Boundary (flux?) value for lower z BC for temperature
    real :: bc_temp_z_val_up = 0.0      !! Boundary (flux?) value for upper z BC for temperature

    real :: bc_sal_y_val_lo = 0.0       !! Boundary (flux?) value for lower y BC for salinity
    real :: bc_sal_y_val_up = 0.0       !! Boundary (flux?) value for upper y BC for salinity
    real :: bc_sal_z_val_lo = 0.0       !! Boundary (flux?) value for lower z BC for salinity
    real :: bc_sal_z_val_up = 0.0       !! Boundary (flux?) value for upper z BC for salinity

    real :: bc_phi_y_val_lo = 0.0       !! Boundary (flux?) value for lower y BC for phase-field
    real :: bc_phi_y_val_up = 0.0       !! Boundary (flux?) value for upper y BC for phase-field
    real :: bc_phi_z_val_lo = 0.0       !! Boundary (flux?) value for lower z BC for phase-field
    real :: bc_phi_z_val_up = 0.0       !! Boundary (flux?) value for upper z BC for phase-field

    real :: bc_humid_y_val_lo = 0.0     !! Boundary (flux?) value for lower y BC for humidity
    real :: bc_humid_y_val_up = 0.0     !! Boundary (flux?) value for upper y BC for humidity
    real :: bc_humid_z_val_lo = 0.0     !! Boundary (flux?) value for lower z BC for humidity
    real :: bc_humid_z_val_up = 0.0     !! Boundary (flux?) value for upper z BC for humidity

    contains

!> Apply a sidewall boundary condition to the variable array `var`
!! ONLY call this routine from an MPI process that contains the relevant boundary
subroutine ApplyBC(var, bc_fix, bc_val, axis, edge, halolvl, gspace)
    real, intent(inout) :: var(:,:,:)   !< Variable array (3D field)
    logical :: bc_fix       !< Flag saying whether the BC is fixed value (.True.) or fixed gradient (.False.)
    real :: bc_val          !< Value to specify at the boundary
    character :: axis       !< specifies the axis ('y' or 'z') of the boundary
    character :: edge       !< specifies which boundary to apply the condition to ('L' or 'U')
    integer :: halolvl      !< Number of halo points beyond boundary
    real :: gspace          !< Grid spacing (needed for gradient conditions)

    integer :: n, nb

    if (edge=='L') then
        !!!
        !   o  o || o  o
        !   1  2    3  4
        !     hlvl
        if (axis=='y') then
            if (bc_fix) then
                do n=1,halolvl
                    var(:,halolvl+1-n,:) = 2.0*bc_val - var(:,halolvl+n,:)
                end do
            else
                do n=1,halolvl
                    var(:,halolvl+1-n,:) = var(:,halolvl+n,:) - gspace*bc_val
                end do
            end if
        elseif (axis=='z') then
            if (bc_fix) then
                do n=1,halolvl
                    var(:,:,halolvl+1-n) = 2.0*bc_val - var(:,:,halolvl+n)
                end do
            else
                do n=1,halolvl
                    var(:,:,halolvl+1-n) = var(:,:,halolvl+n) - gspace*bc_val
                end do
            end if
        end if
    elseif (edge=='U') then
        !!!
        !   o  o || o  o
        !     nb      nb+hlvl
        if (axis=='y') then
            nb = size(var, dim=2) - halolvl
            if (bc_fix) then
                do n=1,halolvl
                    var(:,nb+n,:) = 2.0*bc_val - var(:,nb+1-n,:)
                end do
            else
                do n=1,halolvl
                    var(:,nb+n,:) = var(:,nb+1-n,:) - gspace*bc_val
                end do
            end if
        elseif (axis=='z') then
            nb = size(var, dim=3) - halolvl
            if (bc_fix) then
                do n=1,halolvl
                    var(:,:,nb+n) = 2.0*bc_val - var(:,:,nb+1-n)
                end do
            else
                do n=1,halolvl
                    var(:,:,nb+n) = var(:,:,nb+1-n) - gspace*bc_val
                end do
            end if
        end if
    end if

end subroutine

!> Impose the boundary conditions at the edge of the domain in y
!! if we are using the discrete cosine transform for the pressure
!! solution. By default, we are imposing zero velocity (incl. no slip)
!! and zero flux for scalar fields. In z the default is to
!! impose free-slip walls so that 2D simulations are possible

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
        call ApplyBC(vx, bc_vx_y_fix_lo, bc_vx_y_val_lo, 'y', 'L', lvlhalo, dyy)
        !! Always impose zero normal velocity at walls
        vy(:,1,:) = 0.d0
        vy(:,0,:) = vy(:,2,:)
        call ApplyBC(vz, bc_vz_y_fix_lo, bc_vz_y_val_lo, 'y', 'L', lvlhalo, dyy)
        call ApplyBC(temp, bc_temp_y_fix_lo, bc_temp_y_val_lo, 'y', 'L', lvlhalo, dyy)
        if (salinity) then
            call ApplyBC(sal, bc_sal_y_fix_lo, bc_sal_y_val_lo, 'y', 'L', lvlhalo, dyyr)
        end if
        if (phasefield) then
            call ApplyBC(phi, bc_phi_y_fix_lo, bc_phi_y_val_lo, 'y', 'L', lvlhalo, dyyr)
        end if
        if (moist) then
            call ApplyBC(humid, bc_humid_y_fix_lo, bc_humid_y_val_lo, 'y', 'L', lvlhalo, dyy)
        end if
    end if
    if (xstart(3)==1) then
        call ApplyBC(vx, bc_vx_z_fix_lo, bc_vx_z_val_lo, 'z', 'L', lvlhalo, dzz)
        call ApplyBC(vy, bc_vy_z_fix_lo, bc_vy_z_val_lo, 'z', 'L', lvlhalo, dzz)
        !! Always impose zero normal velocity at walls
        vz(:,:,1) = 0.d0
        vz(:,:,0) = vz(:,:,2)
        call ApplyBC(temp, bc_temp_z_fix_lo, bc_temp_z_val_lo, 'z', 'L', lvlhalo, dzz)
        if (salinity) then
            call ApplyBC(sal, bc_sal_z_fix_lo, bc_sal_z_val_lo, 'z', 'L', lvlhalo, dzzr)
        end if
        if (phasefield) then
            call ApplyBC(phi, bc_phi_z_fix_lo, bc_phi_z_val_lo, 'z', 'L', lvlhalo, dzzr)
        end if
        if (moist) then
            call ApplyBC(humid, bc_humid_z_fix_lo, bc_humid_z_val_lo, 'z', 'L', lvlhalo, dzz)
        end if
    end if

    !! Right wall
    if (xend(2)==nym) then
        call ApplyBC(vx, bc_vx_y_fix_up, bc_vx_y_val_up, 'y', 'U', lvlhalo, dyy)
        !! Always impose zero normal velocity at walls
        vy(:,ny,:) = 0.d0
        ! vy(:,ny+1,:) = vy(:,nym,:)
        call ApplyBC(vz, bc_vz_y_fix_up, bc_vz_y_val_up, 'y', 'U', lvlhalo, dyy)
        call ApplyBC(temp, bc_temp_y_fix_up, bc_temp_y_val_up, 'y', 'U', lvlhalo, dyy)
        if (salinity) then
            call ApplyBC(sal, bc_sal_y_fix_up, bc_sal_y_val_up, 'y', 'U', lvlhalo, dyyr)
        end if
        if (phasefield) then
            call ApplyBC(phi, bc_phi_y_fix_up, bc_phi_y_val_up, 'y', 'U', lvlhalo, dyyr)
        end if
        if (moist) then
            call ApplyBC(humid, bc_humid_y_fix_up, bc_humid_y_val_up, 'y', 'U', lvlhalo, dyy)
        end if
    end if
    if (xend(3)==nzm) then
        call ApplyBC(vx, bc_vx_z_fix_up, bc_vx_z_val_up, 'z', 'U', lvlhalo, dzz)
        call ApplyBC(vy, bc_vy_z_fix_up, bc_vy_z_val_up, 'z', 'U', lvlhalo, dzz)
        !! Always impose zero normal velocity at walls
        vz(:,:,nz) = 0.d0
        ! vz(:,:,nz+1) = vz(:,:,nzm)
        call ApplyBC(temp, bc_temp_z_fix_up, bc_temp_z_val_up, 'z', 'U', lvlhalo, dzz)
        if (salinity) then
            call ApplyBC(sal, bc_sal_z_fix_up, bc_sal_z_val_up, 'z', 'U', lvlhalo, dzzr)
        end if
        if (phasefield) then
            call ApplyBC(phi, bc_phi_z_fix_up, bc_phi_z_val_up, 'z', 'U', lvlhalo, dzzr)
        end if
        if (moist) then
            call ApplyBC(humid, bc_humid_z_fix_up, bc_humid_z_val_up, 'z', 'U', lvlhalo, dzz)
        end if
    end if
end subroutine SetSidewallBCs


!> Read the sidewall.in input file to specify the boundary conditions in y and z
subroutine ReadSidewallInput
    use param
    implicit none
    integer :: i, io
    character(len=4) :: dummy

    open(newunit=io, file='sidewall.in', status='old', action='read')
        do i=1,7
            read(io,301) dummy
        end do
        read(io,*) bc_vx_y_fix_lo, bc_vx_y_val_lo, bc_vx_y_fix_up,  bc_vx_y_val_up
        read(io,301) dummy
        read(io,*) bc_vx_z_fix_lo, bc_vx_z_val_lo, bc_vx_z_fix_up,  bc_vx_z_val_up
        read(io,301) dummy
        read(io,301) dummy
        read(io,301) dummy
        read(io,*) bc_vy_z_fix_lo, bc_vy_z_val_lo, bc_vy_z_fix_up,  bc_vy_z_val_up
        read(io,301) dummy
        read(io,301) dummy
        read(io,301) dummy
        read(io,*) bc_vz_y_fix_lo, bc_vz_y_val_lo, bc_vz_y_fix_up,  bc_vz_y_val_up
        read(io,301) dummy
        read(io,301) dummy
        read(io,301) dummy
        read(io,*) bc_temp_y_fix_lo, bc_temp_y_val_lo, bc_temp_y_fix_up,  bc_temp_y_val_up
        read(io,301) dummy
        read(io,*) bc_temp_z_fix_lo, bc_temp_z_val_lo, bc_temp_z_fix_up,  bc_temp_z_val_up
        read(io,301) dummy
        read(io,301) dummy
        read(io,301) dummy
        read(io,*) bc_sal_y_fix_lo, bc_sal_y_val_lo, bc_sal_y_fix_up,  bc_sal_y_val_up
        read(io,301) dummy
        read(io,*) bc_sal_z_fix_lo, bc_sal_z_val_lo, bc_sal_z_fix_up,  bc_sal_z_val_up
        read(io,301) dummy
        read(io,301) dummy
        read(io,301) dummy
        read(io,*) bc_phi_y_fix_lo, bc_phi_y_val_lo, bc_phi_y_fix_up,  bc_phi_y_val_up
        read(io,301) dummy
        read(io,*) bc_phi_z_fix_lo, bc_phi_z_val_lo, bc_phi_z_fix_up,  bc_phi_z_val_up
        read(io,301) dummy
        read(io,301) dummy
        read(io,301) dummy
        read(io,*) bc_humid_y_fix_lo, bc_humid_y_val_lo, bc_humid_y_fix_up,  bc_humid_y_val_up
        read(io,301) dummy
        read(io,*) bc_humid_z_fix_lo, bc_humid_z_val_lo, bc_humid_z_fix_up,  bc_humid_z_val_up
    301     format(a4)
    close(io)
end subroutine ReadSidewallInput

end module afid_sides