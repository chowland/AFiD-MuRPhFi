subroutine CreateSalStencil
    use param
    use mgrd_arrays
    use input_grids
    use stencil_mod
    implicit none

    ! x direction
    call interpolation_indices(irangsr, xmro, nxmro, xmr, nxmr)
    call construct_stencil(cxrs, xmro, nxmro, xmr, nxmr, irangsr, "x")

    ! y direction
    call interpolation_indices(jrangsr, ymro, nymro, ymr, nymr)
    call construct_stencil(cyrs, ymro, nymro, ymr, nymr, jrangsr, "y")

    ! z direction
    call interpolation_indices(krangsr, zmro, nzmro, zmr, nzmr)
    call construct_stencil(czrs, zmro, nzmro, zmr, nzmr, krangsr, "z")

    return
end subroutine CreateSalStencil