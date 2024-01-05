!> Module containing the routines needed to interpolate
!! between two grids
module HermiteInterpolations
    use decomp_2d, only: xstart, xstartr, xend, xendr
    use mgrd_arrays, only: irangs, jrangs, krangs, &
                            irangr, jrangr, krangr, &
                            irangb, jrangb, krangb, &
                            cxrs, cyrs, czrs, &
                            cxsalc, cysalc, czsalc, &
                            cxphic, cyphic, czphic
    use param, only: nx, nxm, nxmr, lvlhalo, ny, nym, nz, nzm
    use input_grids, only: nxmo, nxmro, xstarto, xendo, &
                            irangsr, jrangsr, krangsr
    implicit none

    private
    public interpolation_indices, construct_stencil, &
            interpolate_xyz_to_refined, &
            interpolate_xyz_to_coarse, &
            interpolate_xyz_to_coarse_fast, &
            interpolate_xyz_old_to_new, &
            interpolate_xyz_old_to_new_ref

contains

!> Construct index array `idx` for the interpolation from grid
!! `x_old` to grid `x_new`
subroutine interpolation_indices(idx, x_old, x_new, x_len)
    integer, intent(out) :: idx(0:)     !! Index array output
    real, intent(in) :: x_old(:)        !! Grid to interpolate from
    real, intent(in) :: x_new(:)        !! Grid being interpolated to
    real, intent(in) :: x_len           !! Upper boundary position
    real, allocatable :: xn(:)
    integer :: n_old, n_new, io, in

    n_old = size(x_old)
    n_new = size(x_new)

    allocate(xn(0:n_new+1))
    xn(1:n_new) = x_new(1:n_new)
    xn(0) = -xn(1)
    if (xn(0)==xn(1)) xn(0) = -xn(2)
    xn(n_new+1) = 2.0*x_len - xn(n_new)
    if (xn(n_new+1)==xn(n_new)) xn(n_new+1) = 2.0*x_len - xn(n_new-1)
    
    idx(0) = 1
    do io=1,n_old
        do in=0,n_new
            if (xn(in) < x_old(io) .and. xn(in+1) >= x_old(io)) then
                idx(io) = in + 1
            end if
        end do
    end do
    idx(n_old+1) = n_new + 1

    deallocate(xn)

end subroutine interpolation_indices

!> Subroutine for extending the grid `x_og` at the lower (x=0) and
!! upper (x=`x_len`) boundaries by `n_ext` points in each direction
subroutine extend_grid(x_ext, x_og, x_len, n_ext)
    real, allocatable, intent(out) :: x_ext(:) !! Extended grid
    real, intent(in) :: x_og(:)    !! Original grid
    real, intent(in) :: x_len   !! Upper boundary position
    integer, intent(in) :: n_ext   !! Number of points to extend by

    integer :: n_og, i

    n_og = size(x_og)
    allocate(x_ext(1-n_ext:n_og+n_ext))

    x_ext(1:n_og) = x_og(1:n_og)

    ! Check if first grid point is on the boundary
    if (abs(x_og(1)) < 1e-10) then
        ! index location and boundary:
        !   -1  0  |1|  2  3
        do i=1,n_ext
            x_ext(1-i) = -x_og(1+i)
        end do
    else
        ! index location and boundary:
        !   -1  0  ||  1  2
        do i=1,n_ext
            x_ext(1-i) = -x_og(i)
        end do
    end if

    ! Check if last grid point is on the boundary
    if (abs(x_og(n_og) - x_len) < 1e-10) then
        ! index location and boundary:
        !   n-2  n-1  |n|  n+1  n+2
        do i=1,n_ext
            x_ext(n_og+i) = 2*x_len - x_og(n_og-i)
        end do
    else
        do i=1,n_ext
            x_ext(n_og+i) = 2*x_len - x_og(n_og+1-i)
        end do
    end if

end subroutine extend_grid

!> Compute stencil coefficients for the interpolation between `x_old`
!! and `x_new` along a given axis. 
subroutine construct_stencil(cx, x_old, x_new, x_len, idx, axis)
    real, intent(out) :: cx(:,:)    !! Stencil coefficient array
    real, intent(in) :: x_old(:)    !! Starting grid
    real, intent(in) :: x_new(:)    !! Grid interpolated to
    real, intent(in) :: x_len       !! Upper boundary position
    integer, intent(in) :: idx(0:)  !! Index array from `interpolation_indices`
    character, intent(in) :: axis   !! Axis of grid ('x', 'y', or 'z')

    integer :: n_old, n_new, io, in, n_ext
    real, allocatable :: xo(:), xn(:)
    real :: t, dlc, dlm, dlp, h00, h01, h10, h11

    if (scan("xyz",axis)==0) then
        write(*,*) 'WARNING: invalid value for axis used in construct_stencil!'
        write(*,*) '         Please set axis as one of "x", "y", or "z".'
    end if

    n_old = size(x_old)
    n_new = size(x_new)

    cx(:,:) = 0.0

    ! Extend grid by 1 point over boundary if x, 2 points if not
    if (axis=="x") then
        n_ext = 1
    else
        n_ext = 2
    end if
    call extend_grid(xo, x_old, x_len, n_ext)
    call extend_grid(xn, x_new, x_len, n_ext)
    
    do io=0,n_old
        ! Use linear interpolation if by solid boundary
        if (axis=="x" .and. (io==0 .or. io==n_old)) then
            dlc = xo(io+1) - xo(io)
            do in=max(idx(io),1),min(idx(io+1)-1,n_new)
                t = (xn(in) - xo(io))/dlc
                cx(1,in) = 0.0
                cx(2,in) = 1 - t
                cx(3,in) = t
                cx(4,in) = 0.0
            end do
        ! Otherwise, use second order Hermite interpolation
        else
            dlm = xo(io) - xo(io-1)
            dlc = xo(io+1) - xo(io)
            dlp = xo(io+2) - xo(io+1)
            do in=max(idx(io),1),min(idx(io+1)-1,n_new)
                t = (xn(in) - xo(io))/dlc
                if (t<1e-12) then
                    cx(:,in) = 0.0
                    cx(2,in) = 1.0
                else
                    h00 = (1.0 + 2.0*t)*(1.0 - t)**2
                    h10 = t*(1.0 - t)**2
                    h01 = (1.0 + 2.0*(1.0 - t))*t**2
                    h11 = -(1.0 - t)*t**2
                    cx(1,in) = -h10*dlc**2/dlm/(dlc + dlm)
                    cx(2,in) = h00 - h11*dlp/(dlp + dlc) &
                                    + h10*(dlc - dlm)/dlm
                    cx(3,in) = h01 + h10*dlm/(dlm + dlc) &
                                    + h11*(dlp - dlc)/dlp
                    cx(4,in) = h11*dlc**2/dlp/(dlp + dlc)
                end if
            end do
        end if
    end do

    deallocate(xo)
    deallocate(xn)

end subroutine construct_stencil

!> Interpolate the 3D array `cvar` from the coarse grid to the refined
!! grid and store the result in `rvar`
subroutine interpolate_xyz_to_refined(cvar, rvar)
    real, dimension(-1:,xstart(2)-lvlhalo:,xstart(3)-lvlhalo:), intent(in) :: cvar
    real, dimension(:,xstartr(2)-lvlhalo:,xstartr(3)-lvlhalo:), intent(out) :: rvar

    real, dimension(4,4,4) :: qv3
    real, dimension(4,4) :: qv2
    real, dimension(4) :: qv1

    integer :: ic, jc, kc, icr, jcr, kcr

    do ic=xstart(3)-1,xend(3)
        do jc=xstart(2)-1,xend(2)
            do kc=0,nxm

                qv3 = cvar(kc-1:kc+2,jc-1:jc+2,ic-1:ic+2)

                do icr=max(krangs(ic),xstartr(3)),min(krangs(ic+1)-1,xendr(3))
                    qv2(:,:) = qv3(:,:,1)*czrs(1,icr) + qv3(:,:,2)*czrs(2,icr) &
                                + qv3(:,:,3)*czrs(3,icr) + qv3(:,:,4)*czrs(4,icr)
                    do jcr=max(jrangs(jc),xstartr(2)),min(jrangs(jc+1)-1,xendr(2))
                        qv1(:) = qv2(:,1)*cyrs(1,jcr) + qv2(:,2)*cyrs(2,jcr) &
                                + qv2(:,3)*cyrs(3,jcr) + qv2(:,4)*cyrs(4,jcr)
                        do kcr=max(irangs(kc),1),min(irangs(kc+1)-1,nxmr)
                            rvar(kcr,jcr,icr) = sum(qv1(1:4)*cxrs(1:4,kcr))
                        end do
                    end do
                end do

            end do
        end do
    end do

end subroutine interpolate_xyz_to_refined

!> Interpolate the 3D array `rvar` from the refined grid to the
!! coarse grid and store the result in `cvar`
subroutine interpolate_xyz_to_coarse(rvar, cvar)
    real, dimension(-1:,xstartr(2)-lvlhalo:,xstartr(3)-lvlhalo:), intent(in) :: rvar
    real, dimension(:,xstart(2)-lvlhalo:,xstart(3)-lvlhalo:), intent(out) :: cvar

    real, dimension(4,4,4) :: qv3
    real, dimension(4,4) :: qv2
    real, dimension(4) :: qv1

    integer :: ic, jc, kc, icr, jcr, kcr

    do icr=xstartr(3)-1,xendr(3)
        do jcr=xstartr(2)-1,xendr(2)
            do kcr=0,nxmr

                qv3 = rvar(kcr-1:kcr+2,jcr-1:jcr+2,icr-1:icr+2)

                do ic=max(krangr(icr),xstart(3)),min(krangr(icr+1)-1,xend(3))
                    qv2(:,:) = qv3(:,:,1)*czphic(1,ic) + qv3(:,:,2)*czphic(2,ic) &
                                + qv3(:,:,3)*czphic(3,ic) + qv3(:,:,4)*czphic(4,ic)
                    do jc=max(jrangr(jcr),xstart(2)),min(jrangr(jcr+1)-1,xend(2))
                        qv1(:) = qv2(:,1)*cyphic(1,jc) + qv2(:,2)*cyphic(2,jc) &
                                + qv2(:,3)*cyphic(3,jc) + qv2(:,4)*cyphic(4,jc)
                        do kc=max(irangr(kcr),1),min(irangr(kcr+1)-1,nxm)
                            cvar(kc,jc,ic) = sum(qv1(1:4)*cxphic(1:4,kc))
                        end do
                    end do
                end do

            end do
        end do
    end do

end subroutine interpolate_xyz_to_coarse

!> A fast implementation of interpolating the 3D array `rvar`
!! from the refined grid to the coarse grid, storing the result
!! in `cvar`. The faster implementation comes from looping over the
!! coarse grid indices rather than refined grid indices.
subroutine interpolate_xyz_to_coarse_fast(rvar, cvar, vname)
    real, dimension(-1:,xstartr(2)-lvlhalo:,xstartr(3)-lvlhalo:), intent(in) :: rvar
    real, dimension(:,xstart(2)-lvlhalo:,xstart(3)-lvlhalo:), intent(out) :: cvar
    character(len=3), intent(in) :: vname

    real, dimension(4,4,4) :: qv3
    real, dimension(4,4) :: qv2
    real, dimension(4) :: qv1

    real, dimension(4,nxm) :: cx
    real, dimension(4,nym) :: cy
    real, dimension(4,nzm) :: cz

    integer, dimension(0:nx) :: irang
    integer, dimension(0:ny) :: jrang
    integer, dimension(0:nz) :: krang

    integer :: ic, jc, kc, icr, jcr, kcr

    if (vname=="sal") then
        cx(:,:) = cxsalc(:,1:nxm)
        cy(:,:) = cysalc(:,:)
        cz(:,:) = czsalc(:,:)
        irang(0:nx) = irangb(0:nx)
        jrang = jrangb
        krang = krangb
    else
        cx(:,:) = cxphic(:,:)
        cy(:,:) = cyphic(:,:)
        cz(:,:) = czphic(:,:)
        irang(0:nx) = irangs(0:nx)
        jrang = jrangs
        krang = krangs
    end if

    do ic=xstart(3),xend(3)
        icr = krang(ic)
        do jc=xstart(2),xend(2)
            jcr = jrang(jc)
            do kc=1,nxm
                kcr = irang(kc)

                qv3 = rvar(kcr-2:kcr+1,jcr-2:jcr+1,icr-2:icr+1)

                qv2(:,:) = qv3(:,:,1)*cz(1,ic) + qv3(:,:,2)*cz(2,ic) &
                        + qv3(:,:,3)*cz(3,ic) + qv3(:,:,4)*cz(4,ic)

                qv1(:) = qv2(:,1)*cy(1,jc) + qv2(:,2)*cy(2,jc) &
                        + qv2(:,3)*cy(3,jc) + qv2(:,4)*cy(4,jc)
                
                cvar(kc,jc,ic) = sum(qv1(1:4)*cx(1:4,kc))

            end do
        end do
    end do

end subroutine interpolate_xyz_to_coarse_fast

    subroutine interpolate_xyz_old_to_new(ovar, nvar)
        real, dimension(-1:,xstarto(2)-lvlhalo:,xstarto(3)-lvlhalo:), intent(in) :: ovar
        real, dimension(:,xstart(2)-lvlhalo:,xstart(3)-lvlhalo:), intent(out) :: nvar

        real, dimension(4,4,4) :: qv3
        real, dimension(4,4) :: qv2
        real, dimension(4) :: qv1

        integer :: ic, jc, kc, icr, jcr, kcr

        do ic=xstarto(3)-1,xendo(3)
            do jc=xstarto(2)-1,xendo(2)
                do kc=0,nxmo
    
                    qv3 = ovar(kc-1:kc+2,jc-1:jc+2,ic-1:ic+2)

                    do icr=max(krangs(ic),xstart(3)),min(krangs(ic+1)-1,xend(3))
                        qv2(:,:) = qv3(:,:,1)*czrs(1,icr) + qv3(:,:,2)*czrs(2,icr) &
                                 + qv3(:,:,3)*czrs(3,icr) + qv3(:,:,4)*czrs(4,icr)
                        do jcr=max(jrangs(jc),xstart(2)),min(jrangs(jc+1)-1,xend(2))
                            qv1(:) = qv2(:,1)*cyrs(1,jcr) + qv2(:,2)*cyrs(2,jcr) &
                                   + qv2(:,3)*cyrs(3,jcr) + qv2(:,4)*cyrs(4,jcr)
                            do kcr=max(irangs(kc),1),min(irangs(kc+1)-1,nxm)
                                nvar(kcr,jcr,icr) = sum(qv1(1:4)*cxrs(1:4,kcr))
                            end do
                        end do
                    end do

                end do
            end do
        end do

    end subroutine interpolate_xyz_old_to_new

    subroutine interpolate_xyz_old_to_new_ref(ovar, nvar)
        real, dimension(-1:,xstarto(2)-lvlhalo:,xstarto(3)-lvlhalo:), intent(in) :: ovar
        real, dimension(:,xstartr(2)-lvlhalo:,xstartr(3)-lvlhalo:), intent(out) :: nvar

        real, dimension(4,4,4) :: qv3
        real, dimension(4,4) :: qv2
        real, dimension(4) :: qv1

        integer :: ic, jc, kc, icr, jcr, kcr

        do ic=xstarto(3)-1,xendo(3)
            do jc=xstarto(2)-1,xendo(2)
                do kc=0,nxmro
    
                    qv3 = ovar(kc-1:kc+2,jc-1:jc+2,ic-1:ic+2)

                    do icr=max(krangsr(ic),xstartr(3)),min(krangsr(ic+1)-1,xendr(3))
                        qv2(:,:) = qv3(:,:,1)*czrs(1,icr) + qv3(:,:,2)*czrs(2,icr) &
                                 + qv3(:,:,3)*czrs(3,icr) + qv3(:,:,4)*czrs(4,icr)
                        do jcr=max(jrangsr(jc),xstartr(2)),min(jrangsr(jc+1)-1,xendr(2))
                            qv1(:) = qv2(:,1)*cyrs(1,jcr) + qv2(:,2)*cyrs(2,jcr) &
                                   + qv2(:,3)*cyrs(3,jcr) + qv2(:,4)*cyrs(4,jcr)
                            do kcr=max(irangsr(kc),1),min(irangsr(kc+1)-1,nxmr)
                                nvar(kcr,jcr,icr) = sum(qv1(1:4)*cxrs(1:4,kcr))
                            end do
                        end do
                    end do

                end do
            end do
        end do

    end subroutine interpolate_xyz_old_to_new_ref


end module HermiteInterpolations