!> This module provides specific parameters and temporary arrays for use with
!! the fast Fourier transform
!! N.B. The bulk of the FFT framework is actually in the modified 2decomp library 2decomp_fft
module afid_fft
    use param, only: sidewall, nym, nzm, ismaster
    use mpih
    use decomp_2d_fft
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

    end Interface

contains

!> Create FFTW plan using DFTs if no sidewalls
!! and using DCTs if `sidewall` is set to true
subroutine PlanFourierTransform
    type(fftw_iodim),dimension(1) :: iodim
    type(fftw_iodim),dimension(2) :: iodim_howmany
    integer :: info
    integer(C_FFTW_R2R_KIND), dimension(1) :: kind_forw, kind_back

    kind_forw(1) = FFTW_REDFT10
    kind_back(1) = FFTW_REDFT01

    iodim(1) % n = nzm
    iodim(1) % is = (sp%zen(1)-sp%zst(1)+1)*(sp%zen(2)-sp%zst(2)+1)
    iodim(1) % os = (sp%zen(1)-sp%zst(1)+1)*(sp%zen(2)-sp%zst(2)+1)

    iodim_howmany(1) % n = (sp%zen(1)-sp%zst(1)+1)
    iodim_howmany(1) % is = 1
    iodim_howmany(1) % os = 1
    iodim_howmany(2) % n = (sp%zen(2)-sp%zst(2)+1)
    iodim_howmany(2) % is = (sp%zen(1)-sp%zst(1)+1)
    iodim_howmany(2) % os = (sp%zen(1)-sp%zst(1)+1)

    ! Construct forward plan for z transform
    if (sidewall) then
        fwd_guruplan_z = fftw_plan_guru_r2r(1, iodim, &
            2, iodim_howmany, rz2, rz2, &
            kind_forw, FFTW_ESTIMATE)
    else
        fwd_guruplan_z = fftw_plan_guru_dft(1, iodim, &
            2, iodim_howmany, cz1, cz1, &
            FFTW_FORWARD, FFTW_ESTIMATE)
    end if
    
    iodim(1) % n = nzm
    ! Construct backward plan for z transform
    if (sidewall) then
        bwd_guruplan_z = fftw_plan_guru_r2r(1, iodim, &
            2,iodim_howmany,rz2,rz2, &
            kind_back,FFTW_ESTIMATE)
    else
        bwd_guruplan_z = fftw_plan_guru_dft(1, iodim, &
            2,iodim_howmany,cz1,cz1, &
            FFTW_BACKWARD,FFTW_ESTIMATE)
    end if
    
    if (.not.c_associated(bwd_guruplan_z)) then
        if (ismaster) print*,'Failed to create guru plan. You should'
        if (ismaster) print*,'link with FFTW3 before MKL'
        if (ismaster) print*,'Please check linking order.'
        call MPI_Abort(MPI_COMM_WORLD,1,info)
    end if
    
    iodim(1) % n = nym
    iodim(1) % is = ph%yen(1)-ph%yst(1)+1
    iodim(1) % os = sp%yen(1)-sp%yst(1)+1
    
    iodim_howmany(1) % n = (ph%yen(1)-ph%yst(1)+1)
    iodim_howmany(1) % is = 1
    iodim_howmany(1) % os = 1
    iodim_howmany(2) % n = (ph%yen(3)-ph%yst(3)+1)
    iodim_howmany(2) % is = (ph%yen(1)-ph%yst(1)+1) &
                          * (ph%yen(2)-ph%yst(2)+1)
    iodim_howmany(2) % os = (sp%yen(1)-sp%yst(1)+1) &
                          * (sp%yen(2)-sp%yst(2)+1)

    ! Construct forward plan for y transform
    if (sidewall) then
        fwd_guruplan_y = fftw_plan_guru_r2r(1, iodim, &
            2, iodim_howmany, ry1, ry2, &
            kind_forw, FFTW_ESTIMATE)
    else
        fwd_guruplan_y = fftw_plan_guru_dft_r2c(1, iodim, &
            2, iodim_howmany, ry1, cy1, &
            FFTW_ESTIMATE)
    end if
    
    iodim(1) % n = nym
    iodim(1) % is = sp%yen(1)-sp%yst(1)+1
    iodim(1) % os = ph%yen(1)-ph%yst(1)+1

    iodim_howmany(1) % n=(sp%yen(1)-sp%yst(1)+1)
    iodim_howmany(1) % is=1
    iodim_howmany(1) % os=1
    iodim_howmany(2) % n=(sp%yen(3)-sp%yst(3)+1)
    iodim_howmany(2) % is=(sp%yen(1)-sp%yst(1)+1) &
                        * (sp%yen(2)-sp%yst(2)+1)
    iodim_howmany(2) % os=(ph%yen(1)-ph%yst(1)+1) &
                        * (ph%yen(2)-ph%yst(2)+1)
    
    ! Construct backward plan for y transform
    if (sidewall) then
        bwd_guruplan_y = fftw_plan_guru_r2r(1,iodim, &
            2,iodim_howmany,ry2,ry1, &
            kind_back, FFTW_ESTIMATE)
    else
        bwd_guruplan_y = fftw_plan_guru_dft_c2r(1,iodim, &
            2,iodim_howmany,cy1,ry1, &
            FFTW_ESTIMATE)
    end if

    ! Save that we have planned the FFT
    planned=.true.

end subroutine PlanFourierTransform

end module