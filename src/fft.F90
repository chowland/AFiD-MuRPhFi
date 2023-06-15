!> This module provides specific parameters and temporary arrays for use with
!! the fast Fourier transform
!! N.B. The bulk of the FFT framework is actually in the modified 2decomp library 2decomp_fft
module afid_fft
    use iso_c_binding
    implicit none

    type, bind(C) :: fftw_iodim
        integer(c_int) :: n, is, os
    end type fftw_iodim

    integer, parameter :: FFTW_PATIENT = 32     !! FFTW patient planner flag
    integer, parameter :: FFTW_ESTIMATE = 64    !! FFTW estimate planner flag
    integer, parameter :: FFTW_FORWARD = -1     !! FFTW sign flag (forward transform)
    integer, parameter :: FFTW_BACKWARD = 1     !! FFTW sign flag (backward transform)
    integer, parameter :: FFTW_REDFT01 = 4      !! FFTW cosine transform method (forward)
    integer, parameter :: FFTW_REDFT10 = 5      !! FFTW cosine transform method (backward)

    integer, parameter :: C_FFTW_R2R_KIND = C_INT32_T   !! 

    type(c_ptr) :: fwd_guruplan_y       !! Plan storage for forward y transform
    type(c_ptr) :: bwd_guruplan_y       !! Plan storage for backward y transform
    type(c_ptr) :: fwd_guruplan_z       !! Plan storage for forward z transform
    type(c_ptr) :: bwd_guruplan_z       !! Plan storage for backward z transform

    logical :: planned = .false.        !! Flag determining whether or not the FFT has been planned yet

    real, allocatable, dimension(:,:,:) :: ry1  !! Temporary real y-pencil array
    real, allocatable, dimension(:,:,:) :: ry2  !! Temporary real y-pencil array (for DCT)
    real, allocatable, dimension(:,:,:) :: rz2  !! Temporary real z-pencil array (for DCT)
    real, allocatable, dimension(:,:,:) :: dphr !! Temporary real x-pencil array (for DCT tridiagonal solve)

    complex, allocatable, dimension(:,:,:) :: cy1       !! Temporary complex y-pencil array (for DFT)
    complex, allocatable, dimension(:,:,:) :: cz1       !! Temporary complex z-pencil array (for DFT)
    complex, allocatable, dimension(:,:,:) :: dphc      !! Temporary complex x-pencil array (for DFT tridiagonal solve)
    complex, allocatable, dimension(:,:,:) :: fouvar1   !! Temporary array for spectra routines
    complex, allocatable, dimension(:,:,:) :: fouvar2   !! Temporary array for spectra routines

    interface
        type(C_PTR) function fftw_plan_guru_dft(rank,dims,howmany_rank,howmany_dims,in,out,sign,flags) &
                bind(C, name='fftw_plan_guru_dft')
            import
            integer(C_INT), value :: rank
            type(fftw_iodim), dimension(*), intent(in) :: dims
            integer(C_INT), value :: howmany_rank
            type(fftw_iodim), dimension(*), intent(in) :: howmany_dims
            complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: in
            complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: out
            integer(C_INT), value :: sign
            integer(C_INT), value :: flags
        end function fftw_plan_guru_dft

        type(C_PTR) function fftw_plan_guru_dft_r2c(rank,dims,howmany_rank,howmany_dims,in,out,flags) &
                bind(C, name='fftw_plan_guru_dft_r2c')
            import
            integer(C_INT), value :: rank
            type(fftw_iodim), dimension(*), intent(in) :: dims
            integer(C_INT), value :: howmany_rank
            type(fftw_iodim), dimension(*), intent(in) :: howmany_dims
            real(C_DOUBLE), dimension(*), intent(out) :: in
            complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: out
            integer(C_INT), value :: flags
        end function fftw_plan_guru_dft_r2c

        type(C_PTR) function fftw_plan_guru_dft_c2r(rank,dims,howmany_rank,howmany_dims,in,out,flags)  &
                bind(C, name='fftw_plan_guru_dft_c2r')
            import
            integer(C_INT), value :: rank
            type(fftw_iodim), dimension(*), intent(in) :: dims
            integer(C_INT), value :: howmany_rank
            type(fftw_iodim), dimension(*), intent(in) :: howmany_dims
            complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: in
            real(C_DOUBLE), dimension(*), intent(out) :: out
            integer(C_INT), value :: flags
        end function fftw_plan_guru_dft_c2r

        type(C_PTR) function fftw_plan_guru_r2r(rank,dims,howmany_rank,howmany_dims,in,out,kind,flags) &
                bind(C, name='fftw_plan_guru_r2r')
            import
            integer(C_INT), value :: rank
            type(fftw_iodim), dimension(*), intent(in) :: dims
            integer(C_INT), value :: howmany_rank
            type(fftw_iodim), dimension(*), intent(in) :: howmany_dims
            real(C_DOUBLE), dimension(*), intent(out) :: in
            real(C_DOUBLE), dimension(*), intent(out) :: out
            integer(C_FFTW_R2R_KIND), dimension(*), intent(in) :: kind
            integer(C_INT), value :: flags
        end function fftw_plan_guru_r2r

    end interface

end module