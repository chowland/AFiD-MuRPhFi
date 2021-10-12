subroutine CreateSalStencil
    use param
    use mgrd_arrays
    use input_grids
    use HermiteInterpolations
    implicit none

    ! x direction
    call interpolation_indices(irangsr, xmro(1:nxmro), xmr(1:nxmr), alx3)
    call construct_stencil(cxrs, xmro(1:nxmro), xmr(1:nxmr), alx3, irangsr, "x")

    ! y direction
    call interpolation_indices(jrangsr, ymro(1:nymro), ymr(1:nymr), ylen)
    call construct_stencil(cyrs, ymro(1:nymro), ymr(1:nymr), ylen, jrangsr, "y")

    ! z direction
    call interpolation_indices(krangsr, zmro(1:nzmro), zmr(1:nzmr), zlen)
    call construct_stencil(czrs, zmro(1:nzmro), zmr(1:nzmr), zlen, krangsr, "z")

    return
end subroutine CreateSalStencil