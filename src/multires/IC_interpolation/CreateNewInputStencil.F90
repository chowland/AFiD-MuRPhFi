subroutine CreateInputStencil
    use param
    use mgrd_arrays
    use input_grids
    use stencil_mod
    implicit none

    ! Construct stencils in x direction
    ! cell edge grid (vx)
    call interpolation_indices(irangc, xco, nxo, xc, nx)
    call construct_stencil(cxvx, xco, nxo, xc, nx, irangc, "x")

    ! cell centre grid (vy, vz, temp)
    call interpolation_indices(irangs, xmo, nxmo, xm, nxm)
    call construct_stencil(cxvy, xmo, nxmo, xm, nxm, irangs, "x")
    cxvz(:,:) = cxvy(:,:)
    cxrs(:,:) = cxvy(:,:)

    ! Construct stencils in y direction
    ! cell edge grid (vy)
    call interpolation_indices(jrangc, yco, nymo, yc, nym)
    call construct_stencil(cyvy, yco, nymo, yc, nym, jrangc, "y")

    ! cell centre grid (vx, vz, temp)
    call interpolation_indices(jrangs, ymo, nymo, ym, nym)
    call construct_stencil(cyvx, ymo, nymo, ym, nym, jrangs, "y")
    cyvz(:,:) = cyvx(:,:)
    cyrs(:,:) = cyvx(:,:)

    ! Construct stencils in z direction
    ! cell edge grid (vz)
    call interpolation_indices(krangc, zco, nzmo, zc, nzm)
    call construct_stencil(czvz, zco, nzmo, zc, nzm, krangc, "z")

    ! cell centre grid (vx, vy, temp)
    call interpolation_indices(krangs, zmo, nzmo, zm, nzm)
    call construct_stencil(czvx, zmo, nzmo, zm, nzm, krangs, "z")
    czvy(:,:) = czvx(:,:)
    czrs(:,:) = czvx(:,:)

    return
end subroutine CreateInputStencil