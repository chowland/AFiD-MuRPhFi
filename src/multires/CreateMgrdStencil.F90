subroutine CreateMgrdStencil
    use param
    use mgrd_arrays
    use HermiteInterpolations
    implicit none

    integer, dimension(0:nxr) :: irange !< temporary index counter in x
    integer, dimension(0:nyr) :: jrange !< temporary index counter in y
    integer, dimension(0:nzr) :: krange !< temporary index counter in z

    ! Construct stencils in x direction
    ! cell edge grid (vx)
    call interpolation_indices(irangc, xc(1:nx), xcr(1:nxr), alx3)
    call construct_stencil(cxvx, xc(1:nx), xcr(1:nxr), alx3, irangc, "x")

    ! cell centre grid (vy, vz, temp)
    call interpolation_indices(irangs, xm(1:nxm), xmr(1:nxmr), alx3)
    call construct_stencil(cxvy, xm(1:nxm), xmr(1:nxmr), alx3, irangs, "x")
    cxvz(:,:) = cxvy(:,:)
    cxrs(:,:) = cxvy(:,:)

    ! refined cell centre grid (sal, phi)
    call interpolation_indices(irangr, xmr(1:nxmr), xm(1:nxm), alx3)
    call construct_stencil(cxphic, xmr(1:nxmr), xm(1:nxm), alx3, irangr, "x")
    if (gAxis==1) then
        call interpolation_indices(irange, xmr(1:nxmr), xc(1:nx), alx3)
        call construct_stencil(cxsalc, xmr(1:nxmr), xc(1:nx), alx3, irange, "x")
        call interpolation_indices(irangb, xc(1:nx), xmr(1:nxmr), alx3)
    else
        irangb = irangs
        cxsalc(:,1:nxm) = cxphic(:,:)
    end if

    ! Construct stencils in y direction
    ! cell edge grid (vy)
    call interpolation_indices(jrangc, yc(1:nym), ycr(1:nymr), ylen)
    call construct_stencil(cyvy, yc(1:nym), ycr(1:nymr), ylen, jrangc, "y")

    ! cell centre grid (vx, vz, temp)
    call interpolation_indices(jrangs, ym(1:nym), ymr(1:nymr), ylen)
    call construct_stencil(cyvx, ym(1:nym), ymr(1:nymr), ylen, jrangs, "y")
    cyvz(:,:) = cyvx(:,:)
    cyrs(:,:) = cyvx(:,:)

    ! refined cell centre grid (sal, phi)
    call interpolation_indices(jrangr, ymr(1:nymr), ym(1:nym), ylen)
    call construct_stencil(cyphic, ymr(1:nymr), ym(1:nym), ylen, jrangr, "y")
    if (gAxis==2) then
        call interpolation_indices(jrange, ymr(1:nymr), yc(1:nym), ylen)
        call construct_stencil(cysalc, ymr(1:nymr), yc(1:nym), ylen, jrange, "y")
        call interpolation_indices(jrangb, yc(1:nym), ymr(1:nymr), ylen)
    else
        jrangb = jrangs
        cysalc(:,:) = cyphic(:,:)
    end if

    ! Construct stencils in z direction
    ! cell edge grid (vz)
    call interpolation_indices(krangc, zc(1:nzm), zcr(1:nzmr), zlen)
    call construct_stencil(czvz, zc(1:nzm), zcr(1:nzmr), zlen, krangc, "z")

    ! cell centre grid (vx, vy, temp)
    call interpolation_indices(krangs, zm(1:nzm), zmr(1:nzmr), zlen)
    call construct_stencil(czvx, zm(1:nzm), zmr(1:nzmr), zlen, krangs, "z")
    czvy(:,:) = czvx(:,:)
    czrs(:,:) = czvx(:,:)

    ! refined cell centre grid (sal, phi)
    call interpolation_indices(krangr, zmr(1:nzmr), zm(1:nzm), zlen)
    call construct_stencil(czphic, zmr(1:nzmr), zm(1:nzm), zlen, krangr, "z")
    if (gAxis==3) then
        call interpolation_indices(krange, zmr(1:nzmr), zc(1:nzm), zlen)
        call construct_stencil(czsalc, zmr(1:nzmr), zc(1:nzm), zlen, krange, "z")
        call interpolation_indices(krangb, zc(1:nzm), zmr(1:nzmr), zlen)
    else
        krangb = krangs
        czsalc(:,:) = czphic(:,:)
    end if

    return
end subroutine CreateMgrdStencil


subroutine CreatePFStencil
    use param
    use mgrd_arrays
    use HermiteInterpolations

    implicit none

    integer, dimension(0:nyr) :: jrange !< temporary index counter in y
    integer, dimension(0:nzr) :: krange !< temporary index counter in z

    ! Save indices for fast interpolation
    call interpolation_indices(yc_to_ymr, yc(1:nym), ymr(1:nymr), ylen)
    call interpolation_indices(zc_to_zmr, zc(1:nzm), zmr(1:nzmr), zlen)

    ! Create temporary index counters for stencil creation
    call interpolation_indices(jrange, ymr(1:nymr), yc(1:nym), ylen)
    call interpolation_indices(krange, zmr(1:nzmr), zc(1:nzm), zlen)

    ! Create stencils
    call construct_stencil(cych, ymr(1:nymr), yc(1:nym), ylen, jrange, "y")
    call construct_stencil(czch, zmr(1:nzmr), zc(1:nzm), zlen, krange, "z")

end subroutine CreatePFStencil