module IBMTools
    use decomp_2d
    use param
    use mgrd_arrays

contains

!! May have to make h(1,:,:) so that we can update halo
!! Or alternatively, include halo in loop, ensuring that phi has an up-to-date halo
subroutine calc_interface_height(ph, h)
    real, intent(in) :: ph(1:,xstartr(2)-lvlhalo:,xstartr(3)-lvlhalo:)
    !< input variable 
    real, intent(out) :: h(xstartr(2)-lvlhalo:,xstartr(3)-lvlhalo:)
    real :: dxx(nxmr)
    integer :: i,j,k,kp

    h(:,:) = 0.0

    dxx = d3xcr(1:nxmr)/dxr
    do i=xstartr(3)-lvlhalo,xendr(3)+lvlhalo
        do j=xstartr(2)-lvlhalo,xendr(2)+lvlhalo
            do k=1,nxmr
                h(j,i) = h(j,i) + ph(k,j,i)*dxx(k)
                ! kp = k + 1
                ! if (ph(k,j,i) <= 0.5 .and. ph(kp,j,i) > 0.5) then
                !     h(j,i) = xmr(k) + (xmr(kp) - xmr(k))*(0.5 - ph(k,j,i))/(ph(kp,j,i) - ph(k,j,i))
                ! end if
            end do
        end do
    end do

end subroutine calc_interface_height

subroutine interp_height_to_vel_grid(hr, hcx, hcy, hcz)
    real, intent(in) :: hr(xstartr(2)-lvlhalo:,xstartr(3)-lvlhalo:)
    real, intent(out) :: hcx(xstart(2):,xstart(3):)
    real, intent(out) :: hcy(xstart(2):,xstart(3):)
    real, intent(out) :: hcz(xstart(2):,xstart(3):)

    real :: hv2(4,4)
    real :: hv1(4)
    integer :: ic, jc, icr, jcr, icc, jcc

    do ic=xstart(3),xend(3)
        icr = krangs(ic)
        icc = zc_to_zmr(ic)
        do jc=xstart(2),xend(2)
            jcr = jrangs(jc)
            jcc = yc_to_ymr(jc)

            hv2 = hr(jcr-2:jcr+1,icr-2:icr+1)
            hv1 = czphic(1,ic)*hv2(:,1) + czphic(2,ic)*hv2(:,2) &
                + czphic(3,ic)*hv2(:,3) + czphic(4,ic)*hv2(:,4)
            hcx(jc,ic) = sum(cyphic(1:4,jc)*hv1)

            hcy(jc,ic) = sum(cych(1:4,jc)*hv1)

            hv1 = czch(1,ic)*hv2(:,1) + czch(2,ic)*hv2(:,2) &
                + czch(3,ic)*hv2(:,3) + czch(4,ic)*hv2(:,4)
            hcz(jc,ic) = sum(cyphic(1:4,jc)*hv1)
        end do
    end do

end subroutine


!> Takes an input height profile `h` (function of y and z)
!>  and outputs a 3D integer array `ibmask` equal to 2 in liquid
!> and equal to 0 in the solid phase
!> grd should be 'x', 'y', or 'z' to denote whether the profile
!> is on the grid for vx, vy, or vz
!> ASSUMES SOLID BELOW LIQUID
subroutine mask_below_height(h, ibmask, grd)
    real, intent(in) :: h(xstart(2):,xstart(3):)
    integer, intent(out) :: ibmask(1:,xstart(2):,xstart(3):)
    character, intent(in) :: grd

    real :: ze, ye, xe
    integer :: i,j,k

    ibmask(:,:,:) = 2

    do i=xstart(3),xend(3)
        ze = zm(i)
        if (grd=='z') ze = zc(i)
        do j=xstart(2),xend(2)
            ye = ym(j)
            if (grd=='y') ye = yc(j)
            do k=1,nxm
                xe = xm(k)
                if (grd=='x') xe = xc(k)
                if (xe < h(j,i)) then
                    ibmask(k,j,i) = 0
                end if
            end do
        end do
    end do

end subroutine mask_below_height

!> Takes an input height profile `h` (function of y and z)
!>  and outputs a 3D integer array `ibmask` equal to 2 in liquid
!> and equal to 0 in the solid phase
!> grd should be 'x', 'y', or 'z' to denote whether the profile
!> is on the grid for vx, vy, or vz
!> ASSUMES LIQUID BELOW SOLID
subroutine mask_above_height(h, ibmask, grd)
    real, intent(in) :: h(xstart(2):,xstart(3):)
    integer, intent(out) :: ibmask(1:,xstart(2):,xstart(3):)
    character, intent(in) :: grd

    real :: ze, ye, xe
    integer :: i,j,k

    ibmask(:,:,:) = 2

    do i=xstart(3),xend(3)
        ze = zm(i)
        if (grd=='z') ze = zc(i)
        do j=xstart(2),xend(2)
            ye = ym(j)
            if (grd=='y') ye = yc(j)
            do k=1,nxm
                xe = alx3 - xm(k)
                if (grd=='x') xe = alx3 - xc(k)
                if (xe < h(j,i)) then
                    ibmask(k,j,i) = 0
                end if
            end do
        end do
    end do

end subroutine mask_above_height

subroutine calc_IBM_interpolation(h, ibmask, dist, grd)
    real, intent(in) :: h(xstart(2):,xstart(3):)
    integer, intent(inout) :: ibmask(1:,xstart(2):,xstart(3):)
    real, intent(out) :: dist(:)
    character, intent(in) :: grd

    real :: xe,ye,ze
    integer :: i,j,k,km,kp,n

    n = 0

    do i=xstart(3),xend(3)
        ze = zm(i)
        if (grd=='z') ze = zc(i)
        do j=xstart(2),xend(2)
            ye = ym(j)
            if (grd=='y') ye = yc(j)
            do k=1,nxm
                km = kmv(k)
                kp = kpv(k)
                xe = xm(k)
                xem = xm(km)
                xep = xm(kp)
                if (grd=='x') then
                    xe = xc(k)
                    xem = xc(km)
                    xep = xc(kp)
                end if

                ! Liquid over solid
                if ((ibmask(k,j,i)==2) .and. (ibmask(km,j,i)==0)) then
                    n = n + 1
                    d1x = xep - xe
                    d2x = xe - h(j,i)
                    dist(n) = d2x/(d1x + d2x)
                    ibmask(k,j,i) = 1

                ! Liquid under solid
                elseif ((ibmask(k,j,i)==2) .and. (ibmask(kp,j,i)==0)) then
                    n = n + 1
                    d1x = xe - xem
                    d2x = alx3 - h(j,i) - xe
                    dist(n) = d2x/(d1x + d2x)
                    ibmask(k,j,i) = -1
                end if
            end do
        end do
    end do
end subroutine calc_IBM_interpolation

end module IBMTools