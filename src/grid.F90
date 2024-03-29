module GridModule
    use param, only: pi
    implicit none

    private
    public uniform_grid, tanh_grid, cheb_grid, asym_cheb_grid, &
            second_derivative_coeff, centre_focus_grid, &
            natural_BL_grid, sym_natural_BL_grid, scallop_grid

contains

    subroutine uniform_grid(c_grd, m_grd, Nm, grd_len)
        real, intent(out) :: c_grd(:), m_grd(:)
        integer, intent(in) :: Nm
        real, intent(in) :: grd_len

        integer :: i

        do i=1,Nm+1
            c_grd(i) = grd_len*real(i-1)/real(Nm)
        end do

        do i=1,Nm
            m_grd(i) = 0.5*(c_grd(i) + c_grd(i+1))
        end do

    end subroutine

    subroutine tanh_grid(c_grd, m_grd, Nm, grd_len, str)
        real, intent(out) :: c_grd(:), m_grd(:)
        integer, intent(in) :: Nm
        real, intent(in) :: grd_len, str

        integer :: i
        real :: tstr, z2

        tstr = tanh(str)

        c_grd(1) = 0.0
        do i=2,Nm+1
            z2 = real(2*(i-1) - Nm)/real(Nm)
            c_grd(i) = 0.5*(1 + tanh(str*z2)/tstr)*grd_len
            if (c_grd(i) < 0.0 .or. c_grd(i) > grd_len) then
                write(*,*) 'Refined grid is too streched: ','c_grd(',i,')=',c_grd(i)
                stop
            end if
        end do

        do i=1,Nm
            m_grd(i) = 0.5*(c_grd(i) + c_grd(i+1))
        end do

    end subroutine

    subroutine cheb_grid(c_grd, m_grd, Nm, grd_len, str)
        real, intent(out) :: c_grd(:), m_grd(:)
        integer, intent(in) :: Nm
        real, intent(in) :: grd_len, str

        integer :: nclip, N_ext, i
        real, allocatable :: etaz(:), etazm(:)
        real :: delet

        nclip = int(str)
        N_ext = Nm + 1 + 2*nclip
        
        allocate(etazm(1:N_ext))
        allocate(etaz(1:Nm+1))

        do i=1,N_ext
            etazm(i) = cos(pi*(real(i) - 0.5)/real(N_ext))
        end do
        do i=1,Nm+1
            etaz(i) = etazm(i+nclip)
        end do
        delet = etaz(1) - etaz(Nm+1)
        do i=1,Nm+1
            etaz(i) = 2.0*etaz(i)/delet
        end do

        c_grd(1) = 0.0
        do i=2,Nm
            c_grd(i) = 0.5*grd_len*(1.0 - etaz(i))
        end do
        c_grd(Nm+1) = grd_len

        do i=1,Nm
            m_grd(i) = 0.5*(c_grd(i) + c_grd(i+1))
        end do

        deallocate(etazm)
        deallocate(etaz)

    end subroutine

    subroutine asym_cheb_grid(c_grd, m_grd, Nm, grd_len, str)
        real, intent(out) :: c_grd(:), m_grd(:)
        integer, intent(in) :: Nm
        real, intent(in) :: grd_len, str

        integer :: nclip, N_ext, i
        real, allocatable :: etaz(:), etazm(:)
        real :: delet

        nclip = int(str)
        N_ext = Nm + 1 + nclip
        
        allocate(etazm(1:N_ext))
        allocate(etaz(1:Nm+1))

        do i=1,N_ext
            etazm(i) = cos(pi*real(i)/real(N_ext)/2.0)
        end do
        do i=1,Nm+1
            etaz(i) = etazm(i+nclip)
        end do
        delet = etaz(1)
        do i=1,Nm+1
            etaz(i) = etaz(i)/delet
        end do

        c_grd(1) = 0.0
        do i=2,Nm
            c_grd(i) = grd_len*(1.0 - etaz(i))
        end do
        c_grd(Nm+1) = grd_len

        do i=1,Nm
            m_grd(i) = 0.5*(c_grd(i) + c_grd(i+1))
        end do

        deallocate(etazm)
        deallocate(etaz)

    end subroutine

    subroutine centre_focus_grid(c_grd, m_grd, Nm, grd_len, str)
        real, intent(out) :: c_grd(:), m_grd(:)
        integer, intent(in) :: Nm
        real, intent(in) :: grd_len, str

        integer :: i

        c_grd(1) = 0.0
        do i=2,Nm
            c_grd(i) = grd_len/2.0*( &
                    str*(2.0*(i-1)/Nm - 1)**7 + 2.0*(i-1)/Nm - 1 &
                )/(str + 1) + grd_len/2.0
        end do
        c_grd(Nm+1) = grd_len

        do i=1,Nm
            m_grd(i) = 0.5*(c_grd(i) + c_grd(i+1))
        end do

    end subroutine centre_focus_grid

    subroutine natural_BL_grid(c_grd, m_grd, Nm, grd_len)
        !! Natural grid space scaling for turbulent boundary layers
        !! following Pirozzoli & Orlandi (2021) J. Comp. Phys.
        real, intent(out) :: c_grd(:), m_grd(:)
        integer, intent(in) :: Nm
        real, intent(in) :: grd_len

        integer :: j, jm
        real :: c_eta, alpha, jb, dyw

        jb = 16.0
        c_eta = 0.8
        alpha = 1.5
        dyw = 0.1

        do j=1,Nm+1
            jm = j - 1
            c_grd(j) = 1/(1.0 + (jm/jb)**2)*(dyw*jm + &
                    (0.75*alpha*c_eta*jm)**(4.0/3.0)*(jm/jb)**2)
        end do
        ! Rescale grid to match domain size
        c_grd(:) = c_grd(:)/c_grd(Nm+1)*grd_len

        do j=1,Nm
            m_grd(j) = 0.5*(c_grd(j) + c_grd(j+1))
        end do

    end subroutine natural_BL_grid

    subroutine sym_natural_BL_grid(c_grd, m_grd, Nm, grd_len, Schmidt)
        !! Natural grid space scaling for turbulent boundary layers
        !! following Pirozzoli & Orlandi (2021) J. Comp. Phys.
        real, intent(out) :: c_grd(:), m_grd(:)
        integer, intent(in) :: Nm
        real, intent(in) :: grd_len
        real, intent(in) :: Schmidt

        integer :: j, jm, nh
        real :: c_eta, alpha, jb, dyw, Lh

        jb = 16.0
        c_eta = 0.8
        alpha = 1.5
        dyw = 0.1

        nh = Nm/2
        do j=1,nh+1
            jm = j - 1
            c_grd(j) = 1/(1.0 + (jm/jb)**2)*(dyw*jm + &
                    (0.75*alpha*c_eta*jm/Schmidt**0.5)**(4.0/3.0)*(jm/jb)**2)
        end do
        Lh = 2.0*c_grd(nh+1)
        do j=1,nh
            c_grd(Nm + 2 - j) = Lh - c_grd(j)
        end do
        ! Rescale grid to match domain size
        c_grd(:) = c_grd(:)/c_grd(Nm+1)*grd_len

        do j=1,Nm
            m_grd(j) = 0.5*(c_grd(j) + c_grd(j+1))
        end do

    end subroutine sym_natural_BL_grid

    subroutine scallop_grid(c_grd, m_grd, Nm, grd_len, Retau, dw)
        !! Natural grid space scaling for turbulent boundary layers
        !! following Pirozzoli & Orlandi (2021) J. Comp. Phys.
        real, intent(out) :: c_grd(:), m_grd(:)
        integer, intent(in) :: Nm
        real, intent(in) :: grd_len, Retau, dw

        integer :: k, ks
        real :: alpha, kb, sig, dxlo, dxup, dxsmooth

        ! Index of roughness height
        kb = 0.2*Retau/dw
        ! Scale to Kolmogorov of upper grid spacing
        alpha = 8.0/(Nm + 1 - kb) * Retau**0.5*(1.0 - 5**(-0.5))
        if (alpha > 2) then
            write(*,*) "WARNING: upper grid spacing predicted >2x Kolmogorov"
        elseif (alpha < 0) then
            write(*,*) "ERROR: wall spacing too small"
        end if

        ! Smoothing region in k
        ks = 20

        do k=1,Nm
            ! Uniform grid spacing in scallop region
            dxlo = dw
            ! Match LK+ == 0.25 x+^0.5
            dxup = alpha/4.0*(alpha/8.0*(k - kb) + (Retau/5.0)**0.5)
            ! Sigmoid function for smooth transition
            sig = 0.5*(1.0 + tanh((k - kb)/ks))
            dxsmooth = (1.0 - sig)*dxlo + sig*dxup
            c_grd(k+1) = c_grd(k) + dxsmooth
        end do
        ! Rescale grid to match domain size
        c_grd(:) = c_grd(:)/c_grd(Nm+1)*grd_len

        do k=1,Nm
            m_grd(k) = 0.5*(c_grd(k) + c_grd(k+1))
        end do

    end subroutine scallop_grid

    subroutine second_derivative_coeff(ap3 ,ac3, am3, x, xlen, fix_up, fix_low)
        ! Returns second derivative coefficients ap3, ac3, am3
        ! Inputs:
        ! -- grid vector x
        ! -- domain size xlen (e.g. alx3, ylen, zlen)
        ! -- fix_up/low (=1 if upper/lower boundary condition is fixed value)
        !               (=0 if upper/lower boundary condition is zero gradient)
        real, intent(out) :: ap3(:), ac3(:), am3(:)
        real, intent(in) :: x(:), xlen
        integer, intent(in) :: fix_up, fix_low

        integer :: k, nx
        real :: a33

        nx = size(x)

        k = 1
        a33 = 2.0/(x(k) + x(k+1))
        ap3(k) = a33/(x(k+1) - x(k))
        am3(k) = 0.0
        ac3(k) = -ap3(k) - fix_low*a33/x(k)
        if (x(k)==0.0) then
            ap3(k) = 0.0
            ac3(k) = 1.0
        end if

        do k=2,nx-1
            a33 = 2.0/(x(k+1) - x(k-1))
            ap3(k) = a33/(x(k+1) - x(k))
            am3(k) = a33/(x(k) - x(k-1))
            ac3(k) = -ap3(k) - am3(k)
        end do

        k = nx
        a33 = 2.0/(2.0*xlen - x(k) - x(k-1))
        ap3(k) = 0.0
        am3(k) = a33/(x(k) - x(k-1))
        ac3(k) = -am3(k) - fix_up*a33/(xlen - x(k))
        if (x(k)==xlen) then
            am3(k) = 0.0
            ac3(k) = 1.0
        end if

    end subroutine second_derivative_coeff

end module GridModule