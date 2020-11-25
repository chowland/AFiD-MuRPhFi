!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: CalcMeanProfiles.f90                           !
!    CONTAINS: subroutine CalcPlateNu                     !
!                                                         ! 
!    PURPOSE: Calculate profiles averaged in the          !
!     wall-parallel directions and output to a file.      !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine CalcMeanProfiles
    
    use param
    use local_arrays, only: vx,vy,vz,temp,sal
    use mgrd_arrays, only: vxr,vyr,vzr
    use mpih
    use decomp_2d, only: xstart,xend,xstartr,xendr
    
    implicit none

    integer :: i, j, k
    real :: inym, inzm, inymr, inzmr, tdx, tdxr
    real :: tii(2)
    real, dimension(nym) :: Tbar, Trms, chiT
    real, dimension(nym) :: vxT, vyT, vzT
    real, dimension(nym) :: vxrms, vyrms, vzrms
    real, dimension(nym) :: vybar, vzbar
    real, dimension(nym) :: epsilon
    real, dimension(nym) :: vxvy, vxvz
    real, dimension(nymr) :: Sbar, Srms, chiS
    real, dimension(nymr) :: vxS, vyS, vzS
    character(5) :: nstat
    character(30) :: dsetname,filename
    logical :: fexist

    tii(1) = MPI_WTIME()

    Tbar(:) =0.0;   vybar(:)=0.0;   vzbar(:)=0.0
    vxT(:)  =0.0;   vyT(:)  =0.0;   vzT(:)  =0.0
    vxrms(:)=0.0;   vyrms(:)=0.0;   vzrms(:)=0.0
    Trms(:) =0.0;   vxvy(:) =0.0;   vxvz(:) =0.0
    chiT(:) =0.0;   epsilon(:)=0.0

    Sbar(:)=0.0;    Srms(:)=0.0;    chiS(:)=0.0
    vxS(:) =0.0;    vyS(:) =0.0;    vzS(:) =0.0

    inym = 1.d0/nym
    inzm = 1.d0/nzm

    do i=xstart(3),xend(3)
        do j=xstart(2),xend(2)
            do k=1,nxm
                Tbar(k) = Tbar(k) + temp(k,j,i)*inym*inzm
                vybar(k) = vybar(k) + 0.5*(vy(k,j,i)+vy(k,j+1,i))*inym*inzm
                vzbar(k) = vzbar(k) + 0.5*(vz(k,j,i)+vz(k,j,i+1))*inym*inzm
                
                vxT(k) = vxT(k) + 0.5*(vx(k,j,i)+vx(k+1,j,i))*temp(k,j,i)*inym*inzm
                vyT(k) = vyT(k) + 0.5*(vy(k,j,i)+vy(k,j+1,i))*temp(k,j,i)*inym*inzm
                vzT(k) = vzT(k) + 0.5*(vz(k,j,i)+vz(k,j,i+1))*temp(k,j,i)*inym*inzm

                Trms(k) = Trms(k) + temp(k,j,i)**2*inym*inzm
                vxrms(k) = vxrms(k) + (0.5*(vx(k,j,i)+vx(k+1,j,i)))**2*inym*inzm
                vyrms(k) = vyrms(k) + (0.5*(vx(k,j,i)+vx(k+1,j,i)))**2*inym*inzm
                vzrms(k) = vzrms(k) + (0.5*(vx(k,j,i)+vx(k+1,j,i)))**2*inym*inzm

                vxvy(k) = vxvy(k) + 0.25*(vx(k,j,i)+vx(k+1,j,i))*(vy(k,j,i)+vy(k,j+1,i))*inym*inzm
                vxvz(k) = vxvz(k) + 0.25*(vx(k,j,i)+vx(k+1,j,i))*(vz(k,j,i)+vz(k,j,i+1))*inym*inzm
            end do
        end do
    end do

    call MpiAllSumReal1D(Tbar,nxm)
    call MpiAllSumReal1D(vybar,nxm)
    call MpiAllSumReal1D(vzbar,nxm)
    call MpiAllSumReal1D(vxT,nxm)
    call MpiAllSumReal1D(vyT,nxm)
    call MpiAllSumReal1D(vzT,nxm)
    call MpiAllSumReal1D(Trms,nxm)
    call MpiAllSumReal1D(vxrms,nxm)
    call MpiAllSumReal1D(vyrms,nxm)
    call MpiAllSumReal1D(vzrms,nxm)
    call MpiAllSumReal1D(vxvy,nxm)
    call MpiAllSumReal1D(vxvz,nxm)

    do i=xstart(3),xend(3)
        do j=xstart(2),xend(2)
            do k=1,nxm
                chiT(k) = chiT(k) + ((temp(k,j,i+1)-temp(k,j,i-1))*0.5*dz)**2*inym*inzm &
                                & + ((temp(k,j+1,i)-temp(k,j-1,i))*0.5*dy)**2*inym*inzm
            end do
            tdx = 0.5*dx/g3rm(1)
            chiT(1) = chiT(1) + ((temp(2,j,i)-temp(1,j,i)+2.0*TfixS*(temp(1,j,i)-tempbp(j,i)))&
                            & *tdx)**2*inym*inzm
            do k=2,nxm-1
                tdx = 0.5*dx/g3rm(k)
                chiT(k) = chiT(k) + ((temp(k+1,j,i)-temp(k-1,j,i)) * tdx)**2*inym*inzm
            end do
            tdx = 0.5*dx/g3rm(nxm)
            chiT(nxm) = chiT(nxm) + ((temp(nxm,j,i)-temp(nxm-1,j,i)+2.0*TfixN*(temptp(j,i)-temp(nxm,j,i)))&
            & *tdx)**2*inym*inzm
        end do
    end do

    call MpiAllSumReal1D(chiT,nxm)

    do i=xstart(3),xend(3)
        do j=xstart(2),xend(2)
            do k=1,nxm
                epsilon(k) = epsilon(k) + ((vz(k,j,i+1)-vz(k,j,i))*dz)**2*inym*inzm &   ! (dw/dz)^2
                                      & + (0.25*(vy(k,j+1,i+1)+vy(k,j,i+1)-vy(k,j+1,i-1)-vy(k,j,i-1))*dz)**2*inym*inzm & ! (dv/dz)^2
                                      & + (0.25*(vx(k+1,j,i+1)+vx(k,j,i+1)-vx(k+1,j,i-1)-vx(k,j,i-1))*dz)**2*inym*inzm & ! (du/dz)^2
                                      & + (0.25*(vz(k,j+1,i+1)+vz(k,j+1,i)-vz(k,j-1,i+1)-vz(k,j-1,i))*dy)**2*inym*inzm & ! (dw/dy)^2
                                      & + ((vy(k,j+1,i)-vy(k,j,i))*dy)**2*inym*inzm &   ! (dv/dy)^2
                                      & + (0.25*(vx(k+1,j+1,i)+vx(k,j+1,i)-vx(k+1,j-1,i)-vx(k,j-1,i))*dy)**2*inym*inzm & ! (du/dy)^2
                                      & + ((vx(k+1,j,i)-vx(k,j,i))*udx3m(k))**2*inym*inzm   ! (du/dx)^2
            end do
            tdx = 0.25*dx/g3rm(1)
            epsilon(1) = epsilon(1) + ((vz(2,j,i+1)+vz(2,j,i)+vz(1,j,i+1)+vz(1,j,i)          )*tdx)**2*inym*inzm &  
                                  & + ((vy(2,j+1,i)+vy(2,j,i)+vy(1,j+1,i)+vy(1,j,i)-2*xminusU)*tdx)**2*inym*inzm
            do k=2,nxm-1
                tdx = 0.25*dx/g3rm(k)
                epsilon(k) = epsilon(k) + ((vz(k+1,j,i+1)+vz(k+1,j,i)-vz(k-1,j,i+1)-vz(k-1,j,i))*tdx)**2*inym*inzm &    ! (dw/dx)^2
                                      & + ((vy(k+1,j+1,i)+vy(k+1,j,i)-vy(k-1,j+1,i)-vy(k-1,j,i))*tdx)**2*inym*inzm      ! (dv/dx)^2
            end do
            tdx = 0.25*dx/g3rm(nxm)
            epsilon(nxm) = epsilon(nxm) + ((        -vz(nxm-1,j,i+1)-vz(nxm-1,j,i)-vz(nxm,j,i+1)-vz(nxm,j,i))*tdx)**2*inym*inzm &
                                      & + ((2*xplusU-vy(nxm-1,j+1,i)-vy(nxm-1,j,i)-vy(nxm,j+1,i)-vy(nxm,j,i))*tdx)**2*inym*inzm
        end do
    end do

    call MpiAllSumReal1D(epsilon,nxm)

    inymr = 1.d0/nymr
    inzmr = 1.d0/nzmr

    do i=xstartr(3),xendr(3)
        do j=xstartr(2),xendr(2)
            do k=1,nxmr
                Sbar(k) = Sbar(k) + sal(k,j,i)*inymr*inzmr

                vxS(k) = vxS(k) + 0.5*(vxr(k,j,i)+vxr(k+1,j,i))*sal(k,j,i)*inymr*inzmr
                vyS(k) = vyS(k) + 0.5*(vyr(k,j,i)+vyr(k,j+1,i))*sal(k,j,i)*inymr*inzmr
                vzS(k) = vzS(k) + 0.5*(vzr(k,j,i)+vzr(k,j,i+1))*sal(k,j,i)*inymr*inzmr

                Srms(k) = Srms(k) + sal(k,j,i)**2*inymr*inzmr
            end do
        end do
    end do

    call MpiAllSumReal1D(Sbar,nxmr)
    call MpiAllSumReal1D(vxS,nxmr)
    call MpiAllSumReal1D(vyS,nxmr)
    call MpiAllSumReal1D(vzS,nxmr)
    call MpiAllSumReal1D(Srms,nxmr)

    do i=xstartr(3),xendr(3)
        do j=xstartr(2),xendr(2)
            do k=1,nxmr
                chiS(k) = chiS(k) + ((sal(k,j,i+1)-sal(k,j,i-1))*0.5*dzr)**2*inymr*inzmr
                chiS(k) = chiS(k) + ((sal(k,j+1,i)-sal(k,j-1,i))*0.5*dyr)**2*inymr*inzmr
            end do
            tdxr = 0.5*dxr/g3rmr(1)
            chiS(1) = chiS(1) + ((sal(2,j,i)-sal(1,j,i)+2.0*SfixS*(sal(1,j,i)-salbp(j,i)))&
                            & *tdx)**2*inymr*inzmr
            do k=2,nxmr-1
                tdx = 0.5*dxr/g3rmr(k)
                chiS(k) = chiS(k) + ((sal(k+1,j,i)-sal(k-1,j,i)) * tdx)**2*inymr*inzmr
            end do
            tdx = 0.5*dxr/g3rmr(nxmr)
            chiS(nxmr) = chiS(nxmr) + ((sal(nxmr,j,i)-sal(nxmr-1,j,i)+2.0*SfixN*(saltp(j,i)-sal(nxmr,j,i)))&
            & *tdx)**2*inymr*inzmr
        end do
    end do

    call MpiAllSumReal1D(chiS,nxmr)

    write(nstat,"(i5.5)")nint(time/tframe)
    filename = trim("outputdir/means.h5")
    
    inquire(file=filename,exist=fexist)
    if (.not.fexist) then
        if (ismaster) then
            call HdfCreateMeansFile(filename)
        end if
    end if

    if (ismaster) then
        dsetname = trim("Tbar/"//nstat)
        call HdfSerialWriteReal1D(dsetname,filename,Tbar,nxm)

        dsetname = trim("vybar/"//nstat)
        call HdfSerialWriteReal1D(dsetname,filename,vybar,nxm)

        dsetname = trim("vzbar/"//nstat)
        call HdfSerialWriteReal1D(dsetname,filename,vzbar,nxm)

        dsetname = trim("vxT/"//nstat)
        call HdfSerialWriteReal1D(dsetname,filename,vxT,nxm)

        dsetname = trim("vyT/"//nstat)
        call HdfSerialWriteReal1D(dsetname,filename,vyT,nxm)

        dsetname = trim("vzT/"//nstat)
        call HdfSerialWriteReal1D(dsetname,filename,vzT,nxm)

        dsetname = trim("Trms/"//nstat)
        call HdfSerialWriteReal1D(dsetname,filename,Trms,nxm)

        dsetname = trim("vxrms/"//nstat)
        call HdfSerialWriteReal1D(dsetname,filename,vxrms,nxm)
        
        dsetname = trim("vyrms/"//nstat)
        call HdfSerialWriteReal1D(dsetname,filename,vyrms,nxm)

        dsetname = trim("vzrms/"//nstat)
        call HdfSerialWriteReal1D(dsetname,filename,vzrms,nxm)

        dsetname = trim("chiT/"//nstat)
        call HdfSerialWriteReal1D(dsetname,filename,chiT,nxm)
        
        dsetname = trim("epsilon/"//nstat)
        call HdfSerialWriteReal1D(dsetname,filename,chiT,nxm)
        
        dsetname = trim("Sbar/"//nstat)
        call HdfSerialWriteReal1D(dsetname,filename,Sbar,nxmr)

        dsetname = trim("vxS/"//nstat)
        call HdfSerialWriteReal1D(dsetname,filename,vxS,nxmr)

        dsetname = trim("vyS/"//nstat)
        call HdfSerialWriteReal1D(dsetname,filename,vyS,nxmr)

        dsetname = trim("vzS/"//nstat)
        call HdfSerialWriteReal1D(dsetname,filename,vzS,nxmr)

        dsetname = trim("Srms/"//nstat)
        call HdfSerialWriteReal1D(dsetname,filename,Srms,nxmr)

        dsetname = trim("chiS/"//nstat)
        call HdfSerialWriteReal1D(dsetname,filename,chiS,nxmr)
    end if

    call MpiBarrier

    tii(2) = MPI_WTIME()
    if (ismaster) write(*,"(a,f8.3,a)") "Profile save duration: ",tii(2)-tii(1),"s"

end subroutine


! Subroutine for creating HDF5 file with structure to contain mean profiles

subroutine HdfCreateMeansFile(filename)
    use hdf5
    
    implicit none

    character(30),intent(in) :: filename
    integer(HID_T) :: file_id, group_id
    integer :: hdf_error

    call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, hdf_error)

    call h5gcreate_f(file_id,"Tbar",group_id,hdf_error)
    call h5gclose_f(group_id,hdf_error)
    call h5gcreate_f(file_id,"vybar",group_id,hdf_error)
    call h5gclose_f(group_id,hdf_error)
    call h5gcreate_f(file_id,"vzbar",group_id,hdf_error)
    call h5gclose_f(group_id,hdf_error)
    call h5gcreate_f(file_id,"vxT",group_id,hdf_error)
    call h5gclose_f(group_id,hdf_error)
    call h5gcreate_f(file_id,"vyT",group_id,hdf_error)
    call h5gclose_f(group_id,hdf_error)
    call h5gcreate_f(file_id,"vzT",group_id,hdf_error)
    call h5gclose_f(group_id,hdf_error)
    call h5gcreate_f(file_id,"Trms",group_id,hdf_error)
    call h5gclose_f(group_id,hdf_error)
    call h5gcreate_f(file_id,"vxrms",group_id,hdf_error)
    call h5gclose_f(group_id,hdf_error)
    call h5gcreate_f(file_id,"vyrms",group_id,hdf_error)
    call h5gclose_f(group_id,hdf_error)
    call h5gcreate_f(file_id,"vzrms",group_id,hdf_error)
    call h5gclose_f(group_id,hdf_error)
    call h5gcreate_f(file_id,"chiT",group_id,hdf_error)
    call h5gclose_f(group_id,hdf_error)
    call h5gcreate_f(file_id,"epsilon",group_id,hdf_error)
    call h5gclose_f(group_id,hdf_error)
    call h5gcreate_f(file_id,"Sbar",group_id,hdf_error)
    call h5gclose_f(group_id,hdf_error)
    call h5gcreate_f(file_id,"vxS",group_id,hdf_error)
    call h5gclose_f(group_id,hdf_error)
    call h5gcreate_f(file_id,"vyS",group_id,hdf_error)
    call h5gclose_f(group_id,hdf_error)
    call h5gcreate_f(file_id,"vzS",group_id,hdf_error)
    call h5gclose_f(group_id,hdf_error)
    call h5gcreate_f(file_id,"Srms",group_id,hdf_error)
    call h5gclose_f(group_id,hdf_error)
    call h5gcreate_f(file_id,"chiS",group_id,hdf_error)
    call h5gclose_f(group_id,hdf_error)
    
    call h5fclose_f(file_id, hdf_error)
end subroutine HdfCreateMeansFile

! subroutine HdfSerialWriteProfile1D(grpname,filename,var,sz)
    
!     use hdf5
!     use param

!     implicit none

!     character(30),intent(in) :: grpname, filename
!     integer, intent(in) :: sz
!     real, dimension(sz), intent(in) :: var
!     integer(HID_T) :: file_id, group_id
!     integer(HID_T) :: dset, filespace
!     integer :: hdf_error
!     integer(HSIZE_T) :: dims(1)
!     character(5) :: dsetname
!     logical :: linkexists

!     call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, hdf_error)
!     call h5gopen_f(file_id, grpname, group_id, hdf_error)

!     dims(1)=sz

!     write(dsetname,"(i5.5)")nint(time/tframe)

!     call h5screate_simple_f(1, dims, filespace, hdf_error)
!     call h5lexists_f(group_id, dsetname, linkexists, hdf_error)
!     if (linkexists) call h5ldelete_f(group_id, dsetname, hdf_error)

!     call h5dcreate_f(group_id, dsetname, H5T_NATIVE_DOUBLE, &
!                         & filespace, dset, hdf_error)
!     call h5dwrite_f(dset, H5T_NATIVE_DOUBLE, var(1:sz), dims, hdf_error)
!     call h5dclose_f(dset, hdf_error)

!     call h5sclose_f(filespace, hdf_error)
!     call h5gclose_f(group_id, hdf_error)
!     call h5fclose_f(file_id, hdf_error)

! end subroutine HdfSerialWriteProfile1D