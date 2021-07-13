subroutine CalcFourierCoef(var,fouvar)
use, intrinsic :: iso_c_binding
use param
use fftw_params
use decomp_2d
use decomp_2d_fft
use mpih
implicit none
integer :: i,j,k,info
integer :: phpiv(nxm)
integer :: nymh
real,intent(in),dimension(1:nxm,xstart(2):xend(2),xstart(3):xend(3)) :: var
complex,intent(inout),dimension(sp%xst(1):sp%xen(1), &
                                sp%xst(2):sp%xen(2), &
                                sp%xst(3):sp%xen(3)) :: fouvar

type(fftw_iodim),dimension(1) :: iodim
type(fftw_iodim),dimension(2) :: iodim_howmany

!RO   Allocate variables for FFT transform

call CreateFFTTmpArrays

nymh=nym/2+1

call transpose_x_to_y(var,ry1,ph)

!RO   Plan FFT transforms if not planned previously

if (.not.planned) then
    iodim(1)%n=nzm
    iodim(1)%is=(sp%zen(1)-sp%zst(1)+1)*(sp%zen(2)-sp%zst(2)+1)
    iodim(1)%os=(sp%zen(1)-sp%zst(1)+1)*(sp%zen(2)-sp%zst(2)+1)
    iodim_howmany(1)%n=(sp%zen(1)-sp%zst(1)+1)
    iodim_howmany(1)%is=1
    iodim_howmany(1)%os=1
    iodim_howmany(2)%n=(sp%zen(2)-sp%zst(2)+1)
    iodim_howmany(2)%is=(sp%zen(1)-sp%zst(1)+1)
    iodim_howmany(2)%os=(sp%zen(1)-sp%zst(1)+1)
    fwd_guruplan_z=fftw_plan_guru_dft(1,iodim, &
            2,iodim_howmany,cz1,cz1, &
            FFTW_FORWARD,FFTW_ESTIMATE)
    iodim(1)%n=nzm
    bwd_guruplan_z=fftw_plan_guru_dft(1,iodim, &
            2,iodim_howmany,cz1,cz1, &
            FFTW_BACKWARD,FFTW_ESTIMATE)

    if (.not.c_associated(bwd_guruplan_z)) then
        if (ismaster) print*,'Failed to create guru plan. You should'
        if (ismaster) print*,'link with FFTW3 before MKL'
        if (ismaster) print*,'Please check linking order.'
        call MPI_Abort(MPI_COMM_WORLD,1,info)
    endif

    iodim(1)%n=nym
    iodim(1)%is=ph%yen(1)-ph%yst(1)+1
    iodim(1)%os=sp%yen(1)-sp%yst(1)+1
    iodim_howmany(1)%n=(ph%yen(1)-ph%yst(1)+1)
    iodim_howmany(1)%is=1
    iodim_howmany(1)%os=1
    iodim_howmany(2)%n=(ph%yen(3)-ph%yst(3)+1)
    iodim_howmany(2)%is=(ph%yen(1)-ph%yst(1)+1) &
            *(ph%yen(2)-ph%yst(2)+1)
    iodim_howmany(2)%os=(sp%yen(1)-sp%yst(1)+1) &
            *(sp%yen(2)-sp%yst(2)+1)
    fwd_guruplan_y=fftw_plan_guru_dft_r2c(1,iodim, &
            2,iodim_howmany,ry1,cy1, &
            FFTW_ESTIMATE)

    iodim(1)%n=nym
    iodim(1)%is=sp%yen(1)-sp%yst(1)+1
    iodim(1)%os=ph%yen(1)-ph%yst(1)+1
    iodim_howmany(1)%n=(sp%yen(1)-sp%yst(1)+1)
    iodim_howmany(1)%is=1
    iodim_howmany(1)%os=1
    iodim_howmany(2)%n=(sp%yen(3)-sp%yst(3)+1)
    iodim_howmany(2)%is=(sp%yen(1)-sp%yst(1)+1) &
            *(sp%yen(2)-sp%yst(2)+1)
    iodim_howmany(2)%os=(ph%yen(1)-ph%yst(1)+1) &
            *(ph%yen(2)-ph%yst(2)+1)
    bwd_guruplan_y=fftw_plan_guru_dft_c2r(1,iodim, &
            2,iodim_howmany,cy1,ry1, &
            FFTW_ESTIMATE)
    planned=.true.
endif

call dfftw_execute_dft_r2c(fwd_guruplan_y,ry1,cy1)

call transpose_y_to_z(cy1,cz1,sp)

call dfftw_execute_dft(fwd_guruplan_z,cz1,cz1)

!EP   Normalize. FFT does not do this
cz1 = cz1 / (nzm*nym)

call transpose_z_to_x(cz1,fouvar,sp)

call DestroyFFTTmpArrays

return
end subroutine CalcFourierCoef
!
!***********************************************************************
subroutine CalcSpec(var1,var2,specy,specz)
use, intrinsic :: iso_c_binding
use param
use fftw_params
use decomp_2d
use decomp_2d_fft
use mpih
implicit none

real,intent(in),dimension(1:nxm,xstart(2):xend(2),xstart(3):xend(3)) :: var1,var2
real,intent(inout),dimension(2,1:int((nym/2+1)*nxm)) :: specy
real,intent(inout),dimension(2,1:int((nzm/2+1)*nxm)) :: specz
real :: re_uhat,re_vhat,im_uhat,im_vhat
real :: re_spec,im_spec
integer :: i,j,k,izk,jyk,niz,nizk

!-- output is specy and specz
specy(:,:)=0.d0
specz(:,:)=0.d0

allocate(fouvar1(sp%xst(1):sp%xen(1), &
                    sp%xst(2):sp%xen(2), &
                    sp%xst(3):sp%xen(3)))
allocate(fouvar2(sp%xst(1):sp%xen(1), &
                    sp%xst(2):sp%xen(2), &
                    sp%xst(3):sp%xen(3)))

!-- Calculate DFT for var1 and var2
call CalcFourierCoef(var1,fouvar1)
call CalcFourierCoef(var2,fouvar2)

!-- Calculate the spectra in y and z (periodic dirs)
do i=sp%xst(3),sp%xen(3)
    do j=sp%xst(2),sp%xen(2)
        do k=1,nxm

            re_uhat= real(fouvar1(k,j,i))
            im_uhat=aimag(fouvar1(k,j,i))
            re_vhat= real(fouvar2(k,j,i))
            im_vhat=aimag(fouvar2(k,j,i))

            re_spec = re_uhat*re_vhat + im_uhat*im_vhat
            im_spec = im_uhat*re_vhat - re_uhat*im_vhat

            izk = (i-1)*nxm+k !-- i counter for specz
            jyk = (j-1)*nxm+k !-- j counter for specy
            niz = mod(nzm-(i-1),nzm)+1
            nizk = (niz-1)*nxm+k

            !-- specy
            if ( i.eq.1 .and. j.eq.1 ) then
            else
                specy(1,jyk) = specy(1,jyk) + re_spec
                specy(2,jyk) = specy(2,jyk) + im_spec
            end if

            !-- specz
            if ( j.eq.1 ) then
                if ( i.gt.1 .and. i.le.(nzm/2+1) ) then
                    specz(1,izk) = specz(1,izk) + re_spec
                    specz(2,izk) = specz(2,izk) + im_spec
                end if
            else if( j.eq.(nym/2+1) ) then
                if ( i.ge.1 .and. i.le.(nzm/2+1) ) then
                    specz(1,izk) = specz(1,izk) + re_spec
                    specz(2,izk) = specz(2,izk) + im_spec
                end if
            else
                if ( i.ge.1 .and. i.le.(nzm/2+1) ) then
                    specz(1,izk) = specz(1,izk) + re_spec
                    specz(2,izk) = specz(2,izk) + im_spec
                end if
                if ( niz.ge.1 .and. niz.le.(nzm/2+1) ) then
                    specz(1,nizk) = specz(1,nizk) + re_spec
                    specz(2,nizk) = specz(2,nizk) - im_spec
                end if
            end if

        end do
    end do
end do

if(allocated(fouvar1)) deallocate(fouvar1)
if(allocated(fouvar2)) deallocate(fouvar2)

!-- Save the files
call MpiSumReal1D(specy(1,:),(nym/2+1)*nxm)
call MpiSumReal1D(specy(2,:),(nym/2+1)*nxm)
call MpiSumReal1D(specz(1,:),(nzm/2+1)*nxm)
call MpiSumReal1D(specz(2,:),(nzm/2+1)*nxm)

! if(ismaster) then
!    open(unit=88,file=filnam1,status='unknown')
!       do k=1,(nym/2+1)*nxm
!          write(88,'(2(1x,e23.15))') specy(1,k), specy(2,k)
!       end do
!    close(88)
!    open(unit=89,file=filnam2,status='unknown')
!       do k=1,(nzm/2+1)*nxm
!          write(89,'(2(1x,e23.15))') specz(1,k), specz(2,k)
!       end do
!    close(89)
! end if

return
end subroutine CalcSpec
!
!***********************************************************************
subroutine CalcPowerSpec(var, idx, specy, specz)
    use, intrinsic :: iso_c_binding
    use param
    use fftw_params
    use decomp_2d
    use decomp_2d_fft
    use mpih
    implicit none

    real,intent(in),dimension(1:nxm,xstart(2):xend(2),xstart(3):xend(3)) :: var
    integer,intent(in),dimension(5) :: idx
    real,intent(inout),dimension(5,1:int(nym/2+1)) :: specy
    real,intent(inout),dimension(5,1:int(nzm/2+1)) :: specz
    real :: re_uhat,re_vhat,im_uhat,im_vhat
    real :: re_spec,im_spec
    integer :: i,j,k,izk,jyk,niz,nizk,kk
    logical :: dsetexists

!-- output is specy and specz
    specy(:,:)=0.d0
    specz(:,:)=0.d0
    
    allocate(fouvar1(sp%xst(1):sp%xen(1), &
                    sp%xst(2):sp%xen(2), &
                    sp%xst(3):sp%xen(3)))
    
    !-- Calculate DFT for var
    call CalcFourierCoef(var,fouvar1)
        
    !-- Calculate the spectra in y and z (periodic dirs)
    do i=sp%xst(3),sp%xen(3)
        do j=sp%xst(2),sp%xen(2)
            do k=1,5
                kk = idx(k)

                re_uhat= real(fouvar1(kk,j,i))
                im_uhat=aimag(fouvar1(kk,j,i))

                re_spec = re_uhat**2 + im_uhat**2
    
                niz = mod(nzm-(i-1),nzm)+1

                !-- specy
                if ( i.eq.1 .and. j.eq.1 ) then
                else
                    specy(k,j) = specy(k,j) + re_spec
                end if

                !-- specz
                if ( j.eq.1 ) then
                    if ( i.gt.1 .and. i.le.(nzm/2+1) ) then
                        specz(k,i) = specz(k,i) + re_spec
                    end if
                else if( j.eq.(nym/2+1) ) then
                    if ( i.ge.1 .and. i.le.(nzm/2+1) ) then
                        specz(k,i) = specz(k,i) + re_spec
                    end if
                else
                    if ( i.ge.1 .and. i.le.(nzm/2+1) ) then
                        specz(k,i) = specz(k,i) + re_spec
                    end if
                    if ( niz.ge.1 .and. niz.le.(nzm/2+1) ) then
                        specz(k,niz) = specz(k,niz) + re_spec
                    end if
                end if
            end do
        end do
    end do

    if(allocated(fouvar1)) deallocate(fouvar1)

!-- Save the files
    do k=1,5
        call MpiSumReal1D(specy(k,:),(nym/2+1))
        call MpiSumReal1D(specz(k,:),(nzm/2+1))
    end do
    
    return
end subroutine CalcPowerSpec
!
!***********************************************************************
subroutine WritePowerSpec
    use param
    use local_arrays, only: vz,vy,vx,temp
    use mgrd_arrays, only: sal
    use decomp_2d, only: xstart,xend,xstartr,xendr
    use stat_arrays
    use mpih
    use hdf5
    implicit none

    character(70) :: filename
    character(70) :: dataname
    character( 5) :: frame
    real,dimension(5,1:int(nym/2+1)) :: vx_specy,vy_specy,vz_specy,te_specy
    real,dimension(5,1:int(nzm/2+1)) :: vx_specz,vy_specz,vz_specz,te_specz
    real,dimension(5,1:int(nymr/2+1)) :: sa_specy
    real,dimension(5,1:int(nzmr/2+1)) :: sa_specz

    real     :: tprfi
    integer  :: ndims
    integer        :: hdf_error
    integer(HID_T) :: file_id
    integer(HID_T) :: filespace
    integer(HID_T) :: dset_var
    integer(HSIZE_T), dimension(2) :: dims
    logical :: dsetexists

    ! calculate spectra
    vx_specy=0.d0; vx_specz=0.d0
    call CalcPowerSpec(vx(1:nxm,xstart(2):xend(2),xstart(3):xend(3)),&
            spec_idx, vx_specy,vx_specz)
    vy_specy=0.d0; vy_specz=0.d0
    call CalcPowerSpec(vy(1:nxm,xstart(2):xend(2),xstart(3):xend(3)),&
            spec_idx, vy_specy,vy_specz)
    vz_specy=0.d0; vz_specz=0.d0
    call CalcPowerSpec(vz(1:nxm,xstart(2):xend(2),xstart(3):xend(3)),&
            spec_idx, vz_specy,vz_specz)
    te_specy=0.d0; te_specz=0.d0
    call CalcPowerSpec(temp(1:nxm,xstart(2):xend(2),xstart(3):xend(3)),&
            spec_idx, te_specy,te_specz)
    ! sa_specy=0.d0; sa_specz=0.d0
    ! call CalcPowerSpecMgrd(sal(1:nxmr,xstartr(2):xendr(2),xstartr(3):xendr(3)),&
    !         spec_idxr, sa_specy,sa_specz)

    ! Record frame number and filename as strings
    write(frame,"(i5.5)")nint(time/tframe)
    filename = trim("outputdir/spectra.h5")

    if (ismaster) then

        !=========
        !  specy
        !=========
        call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, hdf_error)
    
        ndims = 2
        dims(1) = 5
        dims(2) = nym/2+1

        dataname = trim("specy/vx/"//frame)
        call h5lexists_f(file_id, dataname, dsetexists, hdf_error)
        if (dsetexists) call h5ldelete_f(file_id, dataname, hdf_error)
        call h5screate_simple_f(ndims, dims, filespace, hdf_error)
        call h5dcreate_f(file_id, trim(dataname), H5T_NATIVE_DOUBLE, filespace, dset_var, hdf_error)
        call h5dwrite_f(dset_var, H5T_NATIVE_DOUBLE, vx_specy(1:5,1:int(nym/2+1)), dims, hdf_error)
        call h5dclose_f(dset_var, hdf_error)
        call h5sclose_f(filespace, hdf_error)

        dataname = trim("specy/vy/"//frame)
        call h5lexists_f(file_id, dataname, dsetexists, hdf_error)
        if (dsetexists) call h5ldelete_f(file_id, dataname, hdf_error)
        call h5screate_simple_f(ndims, dims, filespace, hdf_error)
        call h5dcreate_f(file_id, trim(dataname), H5T_NATIVE_DOUBLE, filespace, dset_var, hdf_error)
        call h5dwrite_f(dset_var, H5T_NATIVE_DOUBLE, vy_specy(1:5,1:int(nym/2+1)), dims, hdf_error)
        call h5dclose_f(dset_var, hdf_error)
        call h5sclose_f(filespace, hdf_error)

        dataname = trim("specy/vz/"//frame)
        call h5lexists_f(file_id, dataname, dsetexists, hdf_error)
        if (dsetexists) call h5ldelete_f(file_id, dataname, hdf_error)
        call h5screate_simple_f(ndims, dims, filespace, hdf_error)
        call h5dcreate_f(file_id, trim(dataname), H5T_NATIVE_DOUBLE, filespace, dset_var, hdf_error)
        call h5dwrite_f(dset_var, H5T_NATIVE_DOUBLE, vz_specy(1:5,1:int(nym/2+1)), dims, hdf_error)
        call h5dclose_f(dset_var, hdf_error)
        call h5sclose_f(filespace, hdf_error)

        dataname = trim("specy/temp/"//frame)
        call h5lexists_f(file_id, dataname, dsetexists, hdf_error)
        if (dsetexists) call h5ldelete_f(file_id, dataname, hdf_error)
        call h5screate_simple_f(ndims, dims, filespace, hdf_error)
        call h5dcreate_f(file_id, trim(dataname), H5T_NATIVE_DOUBLE, filespace, dset_var, hdf_error)
        call h5dwrite_f(dset_var, H5T_NATIVE_DOUBLE, te_specy(1:5,1:int(nym/2+1)), dims, hdf_error)
        call h5dclose_f(dset_var, hdf_error)
        call h5sclose_f(filespace, hdf_error)

        ! dims(2) = nymr/2 + 1

        ! dataname = trim("specy/sal/"//frame)
        ! call h5lexists_f(file_id, dataname, dsetexists, hdf_error)
        ! if (dsetexists) call h5ldelete_f(file_id, dataname, hdf_error)
        ! call h5screate_simple_f(ndims, dims, filespace, hdf_error)
        ! call h5dcreate_f(file_id, trim(dataname), H5T_NATIVE_DOUBLE, filespace, dset_var, hdf_error)
        ! call h5dwrite_f(dset_var, H5T_NATIVE_DOUBLE, sa_specy(1:5,1:int(nymr/2+1)), dims, hdf_error)
        ! call h5dclose_f(dset_var, hdf_error)
        ! call h5sclose_f(filespace, hdf_error)

        dims(2) = nzm/2+1

        dataname = trim("specz/vx/"//frame)
        call h5lexists_f(file_id, dataname, dsetexists, hdf_error)
        if (dsetexists) call h5ldelete_f(file_id, dataname, hdf_error)
        call h5screate_simple_f(ndims, dims, filespace, hdf_error)
        call h5dcreate_f(file_id, trim(dataname), H5T_NATIVE_DOUBLE, filespace, dset_var, hdf_error)
        call h5dwrite_f(dset_var, H5T_NATIVE_DOUBLE, vx_specz(1:5,1:int(nzm/2+1)), dims, hdf_error)
        call h5dclose_f(dset_var, hdf_error)
        call h5sclose_f(filespace, hdf_error)

        dataname = trim("specz/vy/"//frame)
        call h5lexists_f(file_id, dataname, dsetexists, hdf_error)
        if (dsetexists) call h5ldelete_f(file_id, dataname, hdf_error)
        call h5screate_simple_f(ndims, dims, filespace, hdf_error)
        call h5dcreate_f(file_id, trim(dataname), H5T_NATIVE_DOUBLE, filespace, dset_var, hdf_error)
        call h5dwrite_f(dset_var, H5T_NATIVE_DOUBLE, vy_specz(1:5,1:int(nzm/2+1)), dims, hdf_error)
        call h5dclose_f(dset_var, hdf_error)
        call h5sclose_f(filespace, hdf_error)

        dataname = trim("specz/vz/"//frame)
        call h5lexists_f(file_id, dataname, dsetexists, hdf_error)
        if (dsetexists) call h5ldelete_f(file_id, dataname, hdf_error)
        call h5screate_simple_f(ndims, dims, filespace, hdf_error)
        call h5dcreate_f(file_id, trim(dataname), H5T_NATIVE_DOUBLE, filespace, dset_var, hdf_error)
        call h5dwrite_f(dset_var, H5T_NATIVE_DOUBLE, vz_specz(1:5,1:int(nzm/2+1)), dims, hdf_error)
        call h5dclose_f(dset_var, hdf_error)
        call h5sclose_f(filespace, hdf_error)

        dataname = trim("specz/temp/"//frame)
        call h5lexists_f(file_id, dataname, dsetexists, hdf_error)
        if (dsetexists) call h5ldelete_f(file_id, dataname, hdf_error)
        call h5screate_simple_f(ndims, dims, filespace, hdf_error)
        call h5dcreate_f(file_id, trim(dataname), H5T_NATIVE_DOUBLE, filespace, dset_var, hdf_error)
        call h5dwrite_f(dset_var, H5T_NATIVE_DOUBLE, te_specz(1:5,1:int(nzm/2+1)), dims, hdf_error)
        call h5dclose_f(dset_var, hdf_error)
        call h5sclose_f(filespace, hdf_error)

        ! dims(2) = nzmr/2 + 1

        ! dataname = trim("specz/sal/"//frame)
        ! call h5lexists_f(file_id, dataname, dsetexists, hdf_error)
        ! if (dsetexists) call h5ldelete_f(file_id, dataname, hdf_error)
        ! call h5screate_simple_f(ndims, dims, filespace, hdf_error)
        ! call h5dcreate_f(file_id, trim(dataname), H5T_NATIVE_DOUBLE, filespace, dset_var, hdf_error)
        ! call h5dwrite_f(dset_var, H5T_NATIVE_DOUBLE, sa_specz(1:5,1:int(nzmr/2+1)), dims, hdf_error)
        ! call h5dclose_f(dset_var, hdf_error)
        ! call h5sclose_f(filespace, hdf_error)

        call h5fclose_f(file_id, hdf_error)

    end if

    return

end subroutine WritePowerSpec

!
!***********************************************************************
subroutine WriteSpec
use param
use local_arrays, only: vz,vy,vx
use decomp_2d, only: xstart,xend
use stat_arrays
use mpih
use hdf5
implicit none

character*70 :: specfilename
character*70 :: dataname
real,dimension(2,1:int((nym/2+1)*nxm)) :: vxvx_specy,vyvy_specy
real,dimension(2,1:int((nzm/2+1)*nxm)) :: vxvx_specz,vyvy_specz

real     :: tprfi
integer  :: ndims
integer        :: hdf_error
integer(HID_T) :: file_id
integer(HID_T) :: filespace
integer(HID_T) :: dset_var
integer(HSIZE_T), dimension(2) :: dims

! calculate spectra
vxvx_specy=0.d0; vxvx_specz=0.d0
call CalcSpec(vx(1:nxm,xstart(2):xend(2),xstart(3):xend(3)),&
        vx(1:nxm,xstart(2):xend(2),xstart(3):xend(3)),vxvx_specy,vxvx_specz)
vyvy_specy=0.d0; vyvy_specz=0.d0
call CalcSpec(vy(1:nxm,xstart(2):xend(2),xstart(3):xend(3)),&
        vy(1:nxm,xstart(2):xend(2),xstart(3):xend(3)),vyvy_specy,vyvy_specz)

! write spectra
tprfi = 1.d0/tout ! tframe

if (ismaster) then

    !=========
    !  specy
    !=========
    write(specfilename,'(a)')'outputdir/specyfield_master.h5'
    call h5fopen_f(trim(specfilename), H5F_ACC_RDWR_F, file_id, hdf_error)

    ndims = 2
    dims(1) = 2
    dims(2) = (nym/2+1)*nxm

    write(dataname,'(a,i5.5)')'vxvx_specy/T',nint(time*tprfi)
    call h5screate_simple_f(ndims, dims, filespace, hdf_error)
    call h5dcreate_f(file_id, trim(dataname), H5T_NATIVE_DOUBLE, filespace, dset_var, hdf_error)
    call h5dwrite_f(dset_var, H5T_NATIVE_DOUBLE, vxvx_specy(1:2,1:int((nym/2+1)*nxm)), dims, hdf_error)
    call h5dclose_f(dset_var, hdf_error)
    call h5sclose_f(filespace, hdf_error)

    write(dataname,'(a,i5.5)')'vyvy_specy/T',nint(time*tprfi)
    call h5screate_simple_f(ndims, dims, filespace, hdf_error)
    call h5dcreate_f(file_id, trim(dataname), H5T_NATIVE_DOUBLE, filespace, dset_var, hdf_error)
    call h5dwrite_f(dset_var, H5T_NATIVE_DOUBLE, vyvy_specy(1:2,1:int((nym/2+1)*nxm)), dims, hdf_error)
    call h5dclose_f(dset_var, hdf_error)
    call h5sclose_f(filespace, hdf_error)

    call h5fclose_f(file_id, hdf_error)

    !=========
    !  specz
    !=========
    write(specfilename,'(a)')'outputdir/speczfield_master.h5'
    call h5fopen_f(trim(specfilename), H5F_ACC_RDWR_F, file_id, hdf_error)

    ndims = 2
    dims(1) = 2
    dims(2) = (nzm/2+1)*nxm

    write(dataname,'(a,i5.5)')'vxvx_specz/T',nint(time*tprfi)
    call h5screate_simple_f(ndims, dims, filespace, hdf_error)
    call h5dcreate_f(file_id, trim(dataname), H5T_NATIVE_DOUBLE, filespace, dset_var, hdf_error)
    call h5dwrite_f(dset_var, H5T_NATIVE_DOUBLE, vxvx_specz(1:2,1:int((nzm/2+1)*nxm)), dims, hdf_error)
    call h5dclose_f(dset_var, hdf_error)
    call h5sclose_f(filespace, hdf_error)

    write(dataname,'(a,i5.5)')'vyvy_specz/T',nint(time*tprfi)
    call h5screate_simple_f(ndims, dims, filespace, hdf_error)
    call h5dcreate_f(file_id, trim(dataname), H5T_NATIVE_DOUBLE, filespace, dset_var, hdf_error)
    call h5dwrite_f(dset_var, H5T_NATIVE_DOUBLE, vyvy_specz(1:2,1:int((nzm/2+1)*nxm)), dims, hdf_error)
    call h5dclose_f(dset_var, hdf_error)
    call h5sclose_f(filespace, hdf_error)

    call h5fclose_f(file_id, hdf_error)

end if

return
end subroutine WriteSpec
!
!***********************************************************************
subroutine InitPowerSpec
    use mpih
    use param
    use hdf5
    
    implicit none
    
    integer :: hdf_error, i, k
    integer(HID_T) :: file_id
    integer(HID_T) :: group_id
    integer(HID_T) :: filespace
    character(70) ::  filename
    real :: xloc(5)
    logical :: fexist
    
    open(unit=15, file='spectra.in', status='old')
        read(15,*) xloc(1), xloc(2), xloc(3), xloc(4), xloc(5)
    close(15)
    
    if (ismaster) then
        write(*,*) 'Read from spectra.in:'
        write(*,*) xloc(1), xloc(2), xloc(3), xloc(4), xloc(5)
    end if
    spec_idx(:) = 1
    spec_idxr(:) = 1
    do i=1,5
        k = 1
        do while (xm(k) <= xloc(i))
            spec_idx(i) = k
            k = k + 1
        end do
        do while (xmr(k) <= xloc(i))
            spec_idxr(i) = k
            k = k + 1
        end do
    end do
    if (ismaster) then
        write(*,*) 'Spectra indices, idx: ',spec_idx
        write(*,*) 'Spectra indices, idxr: ',spec_idxr
    end if

    !   create hdf5 file if reset or new simulation
    if (ismaster) then

        filename = trim("outputdir/spectra.h5")
        ! reset file
        inquire(file=filename,exist=fexist)
        if (fexist) then
        else
            call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, hdf_error)
    
            call h5gcreate_f(file_id, 'specy', group_id, hdf_error)
            call h5gclose_f(group_id, hdf_error)
            call h5gcreate_f(file_id, 'specy/vx', group_id, hdf_error)
            call h5gclose_f(group_id, hdf_error)
            call h5gcreate_f(file_id, 'specy/vy', group_id, hdf_error)
            call h5gclose_f(group_id, hdf_error)
            call h5gcreate_f(file_id, 'specy/vz', group_id, hdf_error)
            call h5gclose_f(group_id, hdf_error)
            call h5gcreate_f(file_id, 'specy/temp', group_id, hdf_error)
            call h5gclose_f(group_id, hdf_error)
            call h5gcreate_f(file_id, 'specy/sal', group_id, hdf_error)
            call h5gclose_f(group_id, hdf_error)
            call h5gcreate_f(file_id, 'specz', group_id, hdf_error)
            call h5gclose_f(group_id, hdf_error)
            call h5gcreate_f(file_id, 'specz/vx', group_id, hdf_error)
            call h5gclose_f(group_id, hdf_error)
            call h5gcreate_f(file_id, 'specz/vy', group_id, hdf_error)
            call h5gclose_f(group_id, hdf_error)
            call h5gcreate_f(file_id, 'specz/vz', group_id, hdf_error)
            call h5gclose_f(group_id, hdf_error)
            call h5gcreate_f(file_id, 'specz/temp', group_id, hdf_error)
            call h5gclose_f(group_id, hdf_error)
            call h5gcreate_f(file_id, 'specz/sal', group_id, hdf_error)
            call h5gclose_f(group_id, hdf_error)
    
            call h5fclose_f(file_id, hdf_error)
        end if
    end if
    return
end subroutine InitPowerSpec
!
!***********************************************************************
subroutine InitSpec
use mpih
use param
use hdf5

implicit none

integer :: hdf_error
integer(HID_T) :: file_id
integer(HID_T) :: group_id
integer(HID_T) :: filespace
character*70 ::  specfilename
logical :: tag_exist, grp_exist

!   create hdf5 file if reset or new simulation
if (ismaster) then

    ! reset file
    if(resetlogstime .or. .not.readflow) then

        !=========
        !  specy
        !=========
        write(specfilename,'(a)')'outputdir/specyfield_master.h5'
        call h5fcreate_f(trim(specfilename), H5F_ACC_TRUNC_F, file_id, hdf_error)

        call h5gcreate_f(file_id, 'vxvx_specy', group_id, hdf_error)
        call h5gclose_f(group_id, hdf_error)
        call h5gcreate_f(file_id, 'vyvy_specy', group_id, hdf_error)
        call h5gclose_f(group_id, hdf_error)

        call h5fclose_f(file_id, hdf_error)

        !=========
        !  specz
        !=========
        write(specfilename,'(a)')'outputdir/speczfield_master.h5'
        call h5fcreate_f(trim(specfilename), H5F_ACC_TRUNC_F, file_id, hdf_error)

        call h5gcreate_f(file_id, 'vxvx_specz', group_id, hdf_error)
        call h5gclose_f(group_id, hdf_error)
        call h5gcreate_f(file_id, 'vyvy_specz', group_id, hdf_error)
        call h5gclose_f(group_id, hdf_error)

        call h5fclose_f(file_id, hdf_error)

    endif

    ! continue and check if file exist
    if(.not.resetlogstime .and. readflow) then

        !=========
        !  specy
        !=========
        write(specfilename,'(a)')'outputdir/specyfield_master.h5'
        inquire(file=specfilename,exist=tag_exist)

        ! if not exist, create file
        if(.not.tag_exist) then

            call h5fcreate_f(trim(specfilename), H5F_ACC_TRUNC_F, file_id, hdf_error)

            call h5gcreate_f(file_id, 'vxvx_specy', group_id, hdf_error)
            call h5gclose_f(group_id, hdf_error)
            call h5gcreate_f(file_id, 'vyvy_specy', group_id, hdf_error)
            call h5gclose_f(group_id, hdf_error)

            call h5fclose_f(file_id, hdf_error)

        else if(tag_exist) then

            call h5fopen_f(trim(specfilename), H5F_ACC_RDWR_F, file_id, hdf_error)

            call h5lexists_f(file_id, 'vxvx_specy', grp_exist, hdf_error)
            if(.not.grp_exist)then
                call h5gcreate_f(file_id, 'vxvx_specy', group_id, hdf_error)
                call h5gclose_f(group_id, hdf_error)
            endif

            call h5lexists_f(file_id, 'vyvy_specy', grp_exist, hdf_error)
            if(.not.grp_exist)then
                call h5gcreate_f(file_id, 'vyvy_specy', group_id, hdf_error)
                call h5gclose_f(group_id, hdf_error)
            endif

            call h5fclose_f(file_id, hdf_error)

        endif

        !=========
        !  specz
        !=========
        write(specfilename,'(a)')'outputdir/speczfield_master.h5'
        inquire(file=specfilename,exist=tag_exist)

        ! if not exist, create file
        if(.not.tag_exist) then

            call h5fcreate_f(trim(specfilename), H5F_ACC_TRUNC_F, file_id, hdf_error)

            call h5gcreate_f(file_id, 'vxvx_specz', group_id, hdf_error)
            call h5gclose_f(group_id, hdf_error)
            call h5gcreate_f(file_id, 'vyvy_specz', group_id, hdf_error)
            call h5gclose_f(group_id, hdf_error)

            call h5fclose_f(file_id, hdf_error)

        else if(tag_exist) then

            call h5fopen_f(trim(specfilename), H5F_ACC_RDWR_F, file_id, hdf_error)

            call h5lexists_f(file_id, 'vxvx_specz', grp_exist, hdf_error)
            if(.not.grp_exist)then
                call h5gcreate_f(file_id, 'vxvx_specz', group_id, hdf_error)
                call h5gclose_f(group_id, hdf_error)
            endif

            call h5lexists_f(file_id, 'vyvy_specz', grp_exist, hdf_error)
            if(.not.grp_exist)then
                call h5gcreate_f(file_id, 'vyvy_specz', group_id, hdf_error)
                call h5gclose_f(group_id, hdf_error)
            endif

            call h5fclose_f(file_id, hdf_error)

        endif

    endif
endif

return
end subroutine InitSpec