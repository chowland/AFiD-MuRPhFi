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
    implicit none

    real :: x1, x2, x3
    real :: delet, etain, tstr3, z2dp

    real, allocatable, dimension(:) :: etaz, etazm

    integer :: i, j, k, nxp, nclip

    nxmo = nxo - 1
    nymo = nyo - 1
    nzmo = nzo - 1
    nxmro = nxro - 1
    nymro = nyro - 1
    nzmro = nzro - 1

    ! Uniform grids (horizontal directions)
    
    do i=1,nzo
        x1 = real(i-1)/real(nzmo)
        zco(i) = zlen*x1
    end do
    do i=1,nzmo
        zmo(i) = 0.5d0*(zco(i) + zco(i+1))
    end do

    if (multires) then
        do i=1,nzro
            x1 = real(i-1)/real(nzmro)
            zcro(i) = zlen*x1
        end do
        do i=1,nzmro
            zmro(i) = 0.5d0*(zcro(i) + zcro(i+1))
        end do
        zmro(0) = 2.d0*zmro(1) - zmro(2)
        zmro(nzro) = 2.d0*zmro(nzmro) - zmro(nzmro - 1)
    end if

    do j=1,nyo
        x2 = real(j-1)/real(nymo)
        yco(j) = ylen*x2
    end do
    do j=1,nymo
        ymo(j) = 0.5d0*(yco(j) + yco(j+1))
    end do

    if (multires) then
        do j=1,nyro
            x2 = real(j-1)/real(nymro)
            ycro(j) = ylen*x2
        end do
        do j=1,nymro
            ymro(j) = 0.5d0*(ycro(j) + ycro(j+1))
        end do
        ymro(0) = 2.d0*ymro(1) - ymro(2)
        ymro(nyro) = 2.d0*ymro(nymro) - ymro(nymro - 1)
    end if

    ! Vertical coordinate definition

    if (multires) then
        call AllocateReal1DArray(etaz,1,nxro+500)
        call AllocateReal1DArray(etazm,1,nxro+500)
    else
        call AllocateReal1DArray(etaz,1,nxo+500)
        call AllocateReal1DArray(etazm,1,nxo+500)
    end if

    ! Option 0: Uniform clustering
    if (istr3o.eq.0) then
        do k=1,nxo
            x3 = real(k-1)/real(nxmo)
            etaz(k) = alx3*x3
            xco(k) = etaz(k)
        end do
    end if
    if (multires) then
        if (istr3ro.eq.0) then
            do k=1,nxro
                x3 = real(k-1)/real(nxmro)
                etaz(k) = alx3*x3
                xcro(k) = etaz(k)
            end do
        end if
    end if

    ! Option 4: Hyperbolic tangent-type clustering
    tstr3=tanh(str3o)
    if (istr3o.eq.4) then
        xco(1) = 0.0d0
        do k=2,nxo
            z2dp = float(2*k - nxo - 1)/float(nxmo)
            xco(k) = (1 + tanh(str3*z2dp)/tstr3)*0.5*alx3
        end do
    end if
    if (multires) then
        if (istr3ro.eq.4) then
            xcro(1) = 0.0d0
            do k=2,nxro
                z2dp = float(2*k - nxro - 1)/float(nxmro)
                xcro(k) = (1 + tanh(str3*z2dp)/tstr3)*0.5*alx3
            end do
        end if
    end if

    ! Option 6: Clipped Chebychev-type clustering
    if (istr3o.eq.6) then
        nclip = int(str3o)
        nxp = nxo + nclip + nclip
        do k=1,nxp
            etazm(k) = cos(pi*(float(k) - 0.5)/float(nxp))
        end do
        do k=1,nxo
            etaz(k) = etazm(k + nclip)
        end do
        delet = etaz(1) - etaz(nxo)
        etain = etaz(1)
        do k=1,nxo
            etaz(k) = etaz(k)/(0.5*delet)
        end do
        xco(1) = 0.d0
        do k=2,nxmo
            xco(k) = alx3*(1.d0 -  etaz(k))*0.5
        end do
        xco(nxo) = alx3
    end if
    if (multires) then
        if (istr3ro.eq.6) then
            nclip = int(str3o)
            nxp = nxro + nclip + nclip
            do k=1,nxp
                etazm(k) = cos(pi*(float(k) - 0.5)/float(nxp))
            end do
            do k=1,nxro
                etaz(k) = etazm(k + nclip)
            end do
            delet = etaz(1) - etaz(nxro)
            etain = etaz(1)
            do k=1,nxro
                etaz(k) = etaz(k)/(0.5*delet)
            end do
            xcro(1) = 0.d0
            do k=2,nxmro
                xcro(k) = alx3*(1.d0 -  etaz(k))*0.5
            end do
            xcro(nxro) = alx3
        end if
    end if

    ! Option 7: One-sided clipped Chebychev
    if (istr3o.eq.7) then
        nclip = int(str3o)
        nxp = nxo + nclip
        do k=1,nxp
            etazm(k) = cos(pi*float(k)/float(nxp))
        end do
        do k=1,nxo
            etaz(k) = etazm(k + nclip)
        end do
        delet = etaz(1)
        do k=1,nxo
            etaz(k) = etaz(k)/delet
        end do
        xco(1) = 0.d0
        do k=2,nxmo
            xco(k) = alx3*(1.d0 -  etaz(k))
        end do
        xco(nxo) = alx3
    end if
    if (multires) then
        if (istr3ro.eq.7) then
            nclip = int(str3o)
            nxp = nxro + nclip
            do k=1,nxp
                etazm(k) = cos(pi*float(k)/float(nxp))
            end do
            do k=1,nxro
                etaz(k) = etazm(k + nclip)
            end do
            delet = etaz(1)
            do k=1,nxro
                etaz(k) = etaz(k)/delet
            end do
            xcro(1) = 0.d0
            do k=2,nxmro
                xcro(k) = alx3*(1.d0 -  etaz(k))
            end do
            xcro(nxro) = alx3
        end if
    end if

    call DestroyReal1DArray(etaz)
    call DestroyReal1DArray(etazm)

    do k=1,nxmo
        xmo(k)=(xco(k) + xco(k+1))*0.5d0
    end do
    xmo(nxo) = 2*xco(nxo) - xmo(nxmo)
    if (multires) then
        do k=1,nxmro
            xmro(k)=(xcro(k) + xcro(k+1))*0.5d0
        end do
        xmro(nxro) = 2*xcro(nxro) - xmro(nxmro)
    end if

    return

end subroutine CreateOldGrid