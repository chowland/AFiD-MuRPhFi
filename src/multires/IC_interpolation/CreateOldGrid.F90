!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
!    FILE: CreateOldGrid.F90                              !
!    CONTAINS: subroutine CreateOldGrid                   !
!                                                         !
!    PURPOSE: Compute the grids for an old input          !
!               file, ready for interpolation             !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine CreateOldGrid
    use param
    use input_grids
    use AuxiliaryRoutines
    use GridModule
    implicit none

    nxmo = nxo - 1
    nymo = nyo - 1
    nzmo = nzo - 1
    nxmro = nxro - 1
    nymro = nyro - 1
    nzmro = nzro - 1

    ! Uniform grids (horizontal directions)
    
    call uniform_grid(zco(1:nzo), zmo(1:nzmo), nzmo, zlen)

    if (multires) then
        call uniform_grid(zcro(1:nzro), zmro(1:nzmro), nzmro, zlen)
        zmro(0) = 2.d0*zmro(1) - zmro(2)
        zmro(nzro) = 2.d0*zmro(nzmro) - zmro(nzmro - 1)
    end if

    call uniform_grid(yco(1:nyo), ymo(1:nymo), nymo, ylen)

    if (multires) then
        call uniform_grid(ycro(1:nyro), ymro(1:nymro), nymro, ylen)
        ymro(0) = 2.d0*ymro(1) - ymro(2)
        ymro(nyro) = 2.d0*ymro(nymro) - ymro(nymro - 1)
    end if

    ! Vertical coordinate definition

    ! Option 0: Uniform clustering

    if (istr3o==0) call uniform_grid(xco(1:nxo),xmo(1:nxmo),nxmo,alx3)

    if (multires) then
        if (istr3ro==0) call uniform_grid(xcro(1:nxro),xmro(1:nxmro),nxmro,alx3)
    end if

    ! Option 1: Centre-focused clustering

    if (istr3o==1) call centre_focus_grid(xco(1:nxo),xmo(1:nxmo),nxmo,alx3,str3o)

    if (multires) then
        if (istr3ro==1) call centre_focus_grid(xcro(1:nxro),xmro(1:nxmro),nxmro,alx3,str3o)
    end if

    ! Option 4: Hyperbolic tangent-type clustering

    if (istr3o==4) call tanh_grid(xco(1:nxo),xmo(1:nxmo),nxmo,alx3,str3o)

    if (multires) then
        if (istr3ro==4) call tanh_grid(xcro(1:nxro),xmro(1:nxmro),nxmro,alx3,str3o)
    end if

    ! Option 6: Clipped Chebychev-type clustering

    if (istr3o==6) call cheb_grid(xco(1:nxo),xmo(1:nxmo),nxmo,alx3,str3o)

    if (multires) then
        if (istr3ro==6) call cheb_grid(xcro(1:nxro),xmro(1:nxmro),nxmro,alx3,str3o)
    end if

    ! Option 7: One-sided clipped Chebychev

    if (istr3o==7) call asym_cheb_grid(xco(1:nxo),xmo(1:nxmo),nxmo,alx3,str3o)

    if (multires) then
        if (istr3ro==7) call asym_cheb_grid(xcro(1:nxro),xmro(1:nxmro),nxmro,alx3,str3o)
    end if

    xmo(nxo) = 2*xco(nxo) - xmo(nxmo)
    if (multires) then
        xmro(nxro) = 2*xcro(nxro) - xmro(nxmro)
    end if

    return

end subroutine CreateOldGrid