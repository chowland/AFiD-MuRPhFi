!! Subroutine to update location of the immersed boundary from the phase field variable.
!! Assumes that the solid phase lies beneath the
subroutine UpdateIBMLocation
    use mgrd_arrays, only: phi
    use ibm_param
    use IBMTools
    
    ! Integrate phi vertically to calculate height profile for solid
    call calc_interface_height(phi,solid_height)
    ! if (ismaster) write(*,*) 'height at first index:', solid_height(1,1)

    ! Interpolate solid height profile from refined grid to each velocity grid
    call interp_height_to_vel_grid(solid_height, height_vx, height_vy, height_vz)

    ! Build `ibmaskX` variables from height profile
    ! call mask_below_height(height_vx, ibmaskx, 'x')
    ! call mask_below_height(height_vy, ibmasky, 'y')
    ! call mask_below_height(height_vz, ibmaskz, 'z')
    call mask_above_height(height_vx, ibmaskx, 'x')
    call mask_above_height(height_vy, ibmasky, 'y')
    call mask_above_height(height_vz, ibmaskz, 'z')

    ! For each velocity grid, store the interpolation values for the boundary points
    call calc_IBM_interpolation(height_vx, ibmaskx, distx, 'x')
    call calc_IBM_interpolation(height_vy, ibmasky, disty, 'y')
    call calc_IBM_interpolation(height_vz, ibmaskz, distz, 'z')
    
end subroutine UpdateIBMLocation