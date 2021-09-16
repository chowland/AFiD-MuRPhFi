module stencil_mod
    implicit none

    private
    public interpolation_indices, construct_stencil

contains

    subroutine interpolation_indices(idx, x_old, x_new, x_len)
        integer, intent(out) :: idx(0:)
        real, intent(in) :: x_old(:), x_new(:), x_len
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

    end subroutine

    subroutine construct_stencil(cx, x_old, x_new, x_len, idx, axis)
        real, intent(out) :: cx(:,:)
        real, intent(in) :: x_old(:), x_new(:), x_len
        integer, intent(in) :: idx(0:)
        character, intent(in) :: axis

        integer :: n_old, n_new, io, in
        real, allocatable :: xo(:), xn(:)
        real :: t, dlc, dlm, dlp, h00, h01, h10, h11

        if (scan("xyz",axis)==0) then
            write(*,*) 'WARNING: invalid value for axis used in construct_stencil!'
            write(*,*) '         Please set axis as one of "x", "y", or "z".'
        end if

        n_old = size(x_old)
        n_new = size(x_new)

        if (axis=="x") then
            allocate(xo(0:n_old+1))
            allocate(xn(0:n_new+1))
        else
            allocate(xo(-1:n_old+2))
            allocate(xn(-1:n_new+2))
        end if
        xo(1:n_old) = x_old(1:n_old)
        xn(1:n_new) = x_new(1:n_new)
        if (axis=="x") then
            xo(0) = -xo(2)
            xo(n_old+1) = 2.0*x_len - xo(n_old-1)
            xn(0) = -xn(2)
            xn(n_new+1) = 2.0*x_len - xn(n_new-1)
        else
            xo(0) = -xo(1)
            xo(-1) = -xo(2)
            xo(n_old+1) = 2.0*x_len - xo(n_old)
            xo(n_old+2) = 2.0*x_len - xo(n_old-1)
            xn(0) = -xn(1)
            xn(-1) = -xn(2)
            xn(n_new+1) = 2.0*x_len - xn(n_new)
            xn(n_new+2) = 2.0*x_len - xn(n_new-1)
        end if
        
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
                end do
            end if
        end do

        deallocate(xo)
        deallocate(xn)

    end subroutine

end module stencil_mod