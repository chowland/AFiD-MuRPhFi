subroutine CreateInputStencil
    use param
    use mgrd_arrays
    use input_grids
    use stencil_mod
    implicit none

    ! Construct stencils in x direction
    ! cell edge grid (vx)
    call interpolation_indices(irangc, xco(1:nxo), xc(1:nx), alx3)
    call construct_stencil(cxvx, xco(1:nxo), xc(1:nx), alx3, irangc, "x")

    ! cell centre grid (vy, vz, temp)
    call interpolation_indices(irangs, xmo(1:nxmo), xm(1:nxm), alx3)
    call construct_stencil(cxvy, xmo(1:nxmo), xm(1:nxm), alx3, irangs, "x")
    cxvz(:,:) = cxvy(:,:)
    cxrs(:,:) = cxvy(:,:)

    ! Construct stencils in y direction
    ! cell edge grid (vy)
    call interpolation_indices(jrangc, yco(1:nymo), yc(1:nym), ylen)
    call construct_stencil(cyvy, yco(1:nymo), yc(1:nym), ylen, jrangc, "y")

    ! cell centre grid (vx, vz, temp)
    call interpolation_indices(jrangs, ymo(1:nymo), ym(1:nym), ylen)
    call construct_stencil(cyvx, ymo(1:nymo), ym(1:nym), ylen, jrangs, "y")
    cyvz(:,:) = cyvx(:,:)
    cyrs(:,:) = cyvx(:,:)

    ! Construct stencils in z direction
    ! cell edge grid (vz)
    call interpolation_indices(krangc, zco(1:nzmo), zc(1:nzm), zlen)
    call construct_stencil(czvz, zco(1:nzmo), zc(1:nzm), zlen, krangc, "z")

    ! cell centre grid (vx, vy, temp)
    call interpolation_indices(krangs, zmo(1:nzmo), zm(1:nzm), zlen)
    call construct_stencil(czvx, zmo(1:nzmo), zm(1:nzm), zlen, krangs, "z")
    czvy(:,:) = czvx(:,:)
    czrs(:,:) = czvx(:,:)

    return
end subroutine CreateInputStencil