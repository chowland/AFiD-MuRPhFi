subroutine CreateMgrdStencil
    use param
    use mgrd_arrays
    use stencil_mod
    implicit none

    ! Construct stencils in x direction
    ! cell edge grid (vx)
    call interpolation_indices(irangc, xc, nx, xcr, nxr)
    call construct_stencil(cxvx, xc, nx, xcr, nxr, irangc, "x")

    ! cell centre grid (vy, vz, temp)
    call interpolation_indices(irangs, xm, nxm, xmr, nxmr)
    call construct_stencil(cxvy, xm, nxm, xmr, nxmr, irangs, "x")
    cxvz(:,:) = cxvy(:,:)
    cxrs(:,:) = cxvy(:,:)

    ! refined cell centre grid (sal, phi)
    call interpolation_indices(irangr, xmr, nxmr, xm, nxm)
    call construct_stencil(cxsalc, xmr, nxmr, xm, nxm, irangr, "x")

    ! Construct stencils in y direction
    ! cell edge grid (vy)
    call interpolation_indices(jrangc, yc, nym, ycr, nymr)
    call construct_stencil(cyvy, yc, nym, ycr, nymr, jrangc, "y")

    ! cell centre grid (vx, vz, temp)
    call interpolation_indices(jrangs, ym, nym, ymr, nymr)
    call construct_stencil(cyvx, ym, nym, ymr, nymr, jrangs, "y")
    cyvz(:,:) = cyvx(:,:)
    cyrs(:,:) = cyvx(:,:)

    ! refined cell centre grid (sal, phi)
    call interpolation_indices(jrangr, ymr, nymr, ym, nym)
    call construct_stencil(cysalc, ymr, nymr, ym, nym, jrangr, "y")

    ! Construct stencils in z direction
    ! cell edge grid (vz)
    call interpolation_indices(krangc, zc, nzm, zcr, nzmr)
    call construct_stencil(czvz, zc, nzm, zcr, nzmr, krangc, "z")

    ! cell centre grid (vx, vy, temp)
    call interpolation_indices(krangs, zm, nzm, zmr, nzmr)
    call construct_stencil(czvx, zm, nzm, zmr, nzmr, krangs, "z")
    czvy(:,:) = czvx(:,:)
    czrs(:,:) = czvx(:,:)

    ! refined cell centre grid (sal, phi)
    call interpolation_indices(krangr, zmr, nzmr, zm, nzm)
    call construct_stencil(czsalc, zmr, nzmr, zm, nzm, krangr, "z")

    return
end subroutine CreateMgrdStencil