!> Module for calculating energy/power spectra of various quantities
!! Spectra are stored as time averages and written out at the end of a simulation
module afid_spectra
    use afid_fft
    use afid_pressure
    use decomp_2d_fft

    real, allocatable, dimension(:,:,:) :: vx_spec !! Power spectrum of vx
    real, allocatable, dimension(:,:,:) :: vy_spec !! Power spectrum of vy
    real, allocatable, dimension(:,:,:) :: vz_spec !! Power spectrum of vz
    real, allocatable, dimension(:,:,:) :: te_spec !! Power spectrum of temp

    real, allocatable, dimension(:,:,:) :: wSr_spec !! Cospectrum of wall-normal S-flux (real part)
    real, allocatable, dimension(:,:,:) :: wSi_spec !! Cospectrum of wall-normal S-flux (imaginary part)

contains

!> Allocate memory for 3D arrays of energy spectra
subroutine InitSpectra
    use AuxiliaryRoutines, only: AllocateReal3DArray

    allocate(fouvar1(sp%xst(1):sp%xen(1), &
            sp%xst(2):sp%xen(2),sp%xst(3):sp%xen(3)))
    allocate(fouvar2(sp%xst(1):sp%xen(1), &
            sp%xst(2):sp%xen(2), sp%xst(3):sp%xen(3)))

    call AllocateReal3DArray(vx_spec,1,nxm,sp%xst(2),sp%xen(2),sp%xst(3),sp%xen(3))
    call AllocateReal3DArray(vy_spec,1,nxm,sp%xst(2),sp%xen(2),sp%xst(3),sp%xen(3))
    call AllocateReal3DArray(vz_spec,1,nxm,sp%xst(2),sp%xen(2),sp%xst(3),sp%xen(3))
    call AllocateReal3DArray(te_spec,1,nxm,sp%xst(2),sp%xen(2),sp%xst(3),sp%xen(3))

    if (salinity) then
        call AllocateReal3DArray(wSr_spec,1,nxm,sp%xst(2),sp%xen(2),sp%xst(3),sp%xen(3))
        call AllocateReal3DArray(wSi_spec,1,nxm,sp%xst(2),sp%xen(2),sp%xst(3),sp%xen(3))

        wSr_spec = 0.0
        wSi_spec = 0.0
    end if

    vx_spec = 0.0
    vy_spec = 0.0
    vz_spec = 0.0
    te_spec = 0.0
    
end subroutine InitSpectra

!> Free memory from arrays used to store energy spectra
!! (at end of simulation)
subroutine DeallocateSpectra
    use AuxiliaryRoutines, only: DestroyReal3DArray

    call DestroyReal3DArray(vx_spec)
    call DestroyReal3DArray(vy_spec)
    call DestroyReal3DArray(vz_spec)
    call DestroyReal3DArray(te_spec)
    call DestroyReal3DArray(wSr_spec)
    call DestroyReal3DArray(wSi_spec)

    deallocate(fouvar1, fouvar2)

end subroutine DeallocateSpectra

!> Perform a Fourier transform of the variable `var` in y and z,
!! returning the transformed complex array `fouvar`
subroutine CalcFourierCoef(var, fouvar)
    real, dimension(1:nxm, xstart(2):xend(2), xstart(3):xend(3)), intent(in) :: var
    complex, dimension(sp%xst(1):sp%xen(1), &
                sp%xst(2):sp%xen(2), sp%xst(3):sp%xen(3)), intent(out) :: fouvar

    real :: inyzm

    inyzm = 1.0/(nzm*nym)

    ! Transpose to y-pencil to perform FFT in y locally
    call transpose_x_to_y(var, ry1, ph)
    ! Plan the FFTs if this is not yet done (should be done by init pressure)
    if (.not. planned) call PlanFourierTransform
    ! Perform the FFT in y
    call dfftw_execute_dft_r2c(fwd_guruplan_y, ry1, cy1)

    ! Transpose to z-pencil to perform FFT in z locally
    call transpose_y_to_z(cy1, cz1, sp)
    ! Perform the FFT in z
    call dfftw_execute_dft(fwd_guruplan_z, cz1, cz1)
    ! Normalize the result (FFTW does not do this automatically)
    cz1 = cz1*inyzm

    !CJH: Is this really necessary? Could we work with z-pencils?
    call transpose_z_to_x(cz1, fouvar, sp)

end subroutine CalcFourierCoef

!> Calculate the power spectrum of the variable `var` as a function of
!! height and horizontal wavenumbers
subroutine AddRealSpectrum(var, spectrum)
    use param, only: dt
    real, dimension(1:nxm,xstart(2):xend(2),xstart(3):xend(3)), intent(in) :: var
    real, dimension(1:nxm,sp%xst(2):sp%xen(2),sp%xst(3):sp%xen(3)), intent(inout) :: spectrum
    integer :: i, j, k

    ! Compute FFT on the variable
    call CalcFourierCoef(var, fouvar1)

    do i=sp%xst(3),sp%xen(3)
        do j=sp%xst(2),sp%xen(2)
            do k=1,nxm
                spectrum(k,j,i) = spectrum(k,j,i) + dt*abs(fouvar1(k,j,i))**2
            end do
        end do
    end do

end subroutine AddRealSpectrum

!> Calculate the cospectrum of the variables `var1` and `var2` as a function of
!! height and horizontal wavenumbers and add the real part to `rspec` and the
!! imaginary part to ispec (scaled by dt)
subroutine AddCospectrum(var1, var2, rspec, ispec)
    use param, only: dt
    real, dimension(1:nxm,xstart(2):xend(2),xstart(3):xend(3)), intent(in) :: var1, var2
    real, dimension(1:nxm,sp%xst(2):sp%xen(2),sp%xst(3):sp%xen(3)), intent(inout) :: rspec, ispec
    complex :: cspec
    integer :: i, j, k

    ! Compute FFT on the variable
    call CalcFourierCoef(var1, fouvar1)
    call CalcFourierCoef(var2, fouvar2)

    do i=sp%xst(3),sp%xen(3)
        do j=sp%xst(2),sp%xen(2)
            do k=1,nxm
                cspec = fouvar1(k,j,i)*conjg(fouvar2(k,j,i))
                rspec(k,j,i) = rspec(k,j,i) + dt*real(cspec)
                ispec(k,j,i) = ispec(k,j,i) + dt*aimag(cspec)
            end do
        end do
    end do

end subroutine AddCospectrum

!> Update the spectra
subroutine UpdateSpectra
    use local_arrays, only: vx, vy, vz, temp, dph, dq
    use afid_salinity, only: salc

    call AddRealSpectrum(vx(1:nxm,xstart(2):xend(2),xstart(3):xend(3)), vx_spec)
    call AddRealSpectrum(vy(1:nxm,xstart(2):xend(2),xstart(3):xend(3)), vy_spec)
    call AddRealSpectrum(vz(1:nxm,xstart(2):xend(2),xstart(3):xend(3)), vz_spec)
    call AddRealSpectrum(temp(1:nxm,xstart(2):xend(2),xstart(3):xend(3)), te_spec)

    ! Interpolate vx and salc to the cell centre, temporarily storing in dph
    do i=xstart(3),xend(3)
        do j=xstart(2),xend(2)
            do k=1,nxm
                dph(k,j,i) = 0.5*(vx(k,j,i) + vx(k+1,j,i))
                dq(k,j,i) = 0.5*(salc(k,j,i) + salc(k,j+1,i))
            end do
        end do
    end do

    if (salinity) then
        call AddCospectrum(dq(1:nxm,xstart(2):xend(2),xstart(3):xend(3)), &
                    dph(1:nxm,xstart(2):xend(2),xstart(3):xend(3)), &
                    wSr_spec, wSi_spec)
    end if

end subroutine UpdateSpectra

!> Normalize the time-averaged spectra and write them out to files
subroutine WriteSpectra
    use h5_tools, only: write_3D_spectrum
    use afid_averaging, only: tav_start
    real :: tav_length  !! Duration of time averaging interval
    real :: i_tav
    character(len=30) :: filename

    tav_length = time - tav_start
    i_tav = 1.0/tav_length

    vx_spec = vx_spec*i_tav
    vy_spec = vy_spec*i_tav
    vz_spec = vz_spec*i_tav
    te_spec = te_spec*i_tav
    if (salinity) then
        wSr_spec = wSr_spec*i_tav
        wSi_spec = wSi_spec*i_tav
    end if

    ! Save the arrays
    filename = 'outputdir/vx_spectrum.h5'
    call write_3D_spectrum(filename, vx_spec)
    filename = 'outputdir/vy_spectrum.h5'
    call write_3D_spectrum(filename, vy_spec)
    filename = 'outputdir/vz_spectrum.h5'
    call write_3D_spectrum(filename, vz_spec)
    filename = 'outputdir/te_spectrum.h5'
    call write_3D_spectrum(filename, te_spec)

    if (salinity) then
        filename = 'outputdir/wSr_spectrum.h5'
        call write_3D_spectrum(filename, wSr_spec)
        filename = 'outputdir/wSi_spectrum.h5'
        call write_3D_spectrum(filename, wSi_spec)
    end if

end subroutine WriteSpectra
    
end module afid_spectra