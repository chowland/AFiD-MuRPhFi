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
    use local_arrays, only: vx,vy,vz,temp
    use mpih
    use decomp_2d, only: xstart,xend
    use afid_salinity, only: CalcSalStats, CreateSalinityH5Groups
    use afid_moisture, only: CalcMoistStats, CreateMoistH5Groups
    use afid_phasefield, only: CalcPhiStats, CreatePhaseH5Groups
    
    implicit none

    integer :: i, j, k, im, ip, jm, jp, km, kp
    real :: inym, inzm, tdx
    real :: tii(2)
    real, dimension(nxm) :: Tbar, Trms, chiT!, chiT2
    real, dimension(nxm) :: vxT, vyT, vzT
    real, dimension(nxm) :: vxrms, vyrms, vzrms
    real, dimension(nxm) :: vybar, vzbar
    real, dimension(nxm) :: epsilon
    real, dimension(nxm) :: vxvy, vxvz, vyvz

    character(5) :: nstat
    character(30) :: dsetname,filename
    logical :: fexist

    tii(1) = MPI_WTIME()

    Tbar(:) =0.0;   vybar(:)=0.0;   vzbar(:)=0.0
    vxT(:)  =0.0;   vyT(:)  =0.0;   vzT(:)  =0.0
    vxrms(:)=0.0;   vyrms(:)=0.0;   vzrms(:)=0.0
    Trms(:) =0.0;   vxvy(:) =0.0;   vxvz(:) =0.0
    chiT(:) =0.0;   epsilon(:)=0.0; vyvz(:) =0.0

    inym = 1.d0/nym
    inzm = 1.d0/nzm

    do i=xstart(3),xend(3)
        ip = i + 1
        do j=xstart(2),xend(2)
            jp = j + 1
            do k=1,nxm
                kp = k + 1
                Tbar(k) = Tbar(k) + temp(k,j,i)
                vybar(k) = vybar(k) + 0.5*(vy(k,j,i)+vy(k,jp,i))
                vzbar(k) = vzbar(k) + 0.5*(vz(k,j,i)+vz(k,j,ip))
                
                vxT(k) = vxT(k) + 0.5*(vx(k,j,i)+vx(kp,j,i))*temp(k,j,i)
                vyT(k) = vyT(k) + 0.5*(vy(k,j,i)+vy(k,jp,i))*temp(k,j,i)
                vzT(k) = vzT(k) + 0.5*(vz(k,j,i)+vz(k,j,ip))*temp(k,j,i)

                Trms(k) = Trms(k) + temp(k,j,i)**2
                vxrms(k) = vxrms(k) + 0.5*(vx(k,j,i)**2+vx(kp,j,i)**2)
                vyrms(k) = vyrms(k) + 0.5*(vy(k,j,i)**2+vy(k,jp,i)**2)
                vzrms(k) = vzrms(k) + 0.5*(vz(k,j,i)**2+vz(k,j,ip)**2)

                vxvy(k) = vxvy(k) + 0.25*(vx(k,j,i)+vx(kp,j,i))*(vy(k,j,i)+vy(k,jp,i))
                vxvz(k) = vxvz(k) + 0.25*(vx(k,j,i)+vx(kp,j,i))*(vz(k,j,i)+vz(k,j,ip))
                vyvz(k) = vyvz(k) + 0.25*(vy(k,j,i)+vy(k,jp,i))*(vz(k,j,i)+vz(k,j,ip))
            end do
        end do
    end do

    call MpiSumReal1D(Tbar,nxm)
    call MpiSumReal1D(vybar,nxm)
    call MpiSumReal1D(vzbar,nxm)
    call MpiSumReal1D(vxT,nxm)
    call MpiSumReal1D(vyT,nxm)
    call MpiSumReal1D(vzT,nxm)
    call MpiSumReal1D(Trms,nxm)
    call MpiSumReal1D(vxrms,nxm)
    call MpiSumReal1D(vyrms,nxm)
    call MpiSumReal1D(vzrms,nxm)
    call MpiSumReal1D(vxvy,nxm)
    call MpiSumReal1D(vxvz,nxm)
    call MpiSumReal1D(vyvz,nxm)

    do k=1,nxm
        Tbar(k) = Tbar(k)*inym*inzm
        vybar(k) = vybar(k)*inym*inzm
        vzbar(k) = vzbar(k)*inym*inzm
        vxT(k) = vxT(k)*inym*inzm
        vyT(k) = vyT(k)*inym*inzm
        vzT(k) = vzT(k)*inym*inzm
        Trms(k) = sqrt(Trms(k)*inym*inzm)
        vxrms(k) = sqrt(vxrms(k)*inym*inzm)
        vyrms(k) = sqrt(vyrms(k)*inym*inzm)
        vzrms(k) = sqrt(vzrms(k)*inym*inzm)
        vxvy(k) = vxvy(k)*inym*inzm
        vxvz(k) = vxvz(k)*inym*inzm
        vyvz(k) = vyvz(k)*inym*inzm
    end do

    do i=xstart(3),xend(3)
        ip = i + 1
        im = i - 1
        do j=xstart(2),xend(2)
            jp = j + 1
            jm = j - 1
            do k=1,nxm
                chiT(k) = chiT(k) + ((temp(k,j,ip)-temp(k,j,im))*0.5*dz)**2 &
                                & + ((temp(k,jp,i)-temp(k,jm,i))*0.5*dy)**2
            end do
            tdx = 0.5*dx/g3rm(1)
            chiT(1) = chiT(1) + ((temp(2,j,i)-temp(1,j,i)+2.0*TfixS*(temp(1,j,i)-tempbp(1,j,i)))&
                            & *tdx)**2
            do k=2,nxm-1
                km = k - 1
                kp = k + 1
                tdx = 0.5*dx/g3rm(k)
                chiT(k) = chiT(k) + ((temp(kp,j,i)-temp(km,j,i)) * tdx)**2
            end do
            tdx = 0.5*dx/g3rm(nxm)
            chiT(nxm) = chiT(nxm) + ((temp(nxm,j,i)-temp(nxm-1,j,i)+2.0*TfixN*(temptp(1,j,i)-temp(nxm,j,i)))&
                            & *tdx)**2
        end do
    end do

    call MpiSumReal1D(chiT,nxm)

    ! do i=xstart(3),xend(3)
    !     ip = i + 1
    !     im = i - 1
    !     do j=xstart(2),xend(2)
    !         jp = j + 1
    !         jm = j - 1
    !         do k=1,nxm
    !             chiT2(k) = chiT2(k) + &     ! (dT/dz)^2 + (dT/dy)^2
    !                     0.5*dz**2*((temp(k,j,ip) - temp(k,j,i))**2 + (temp(k,j,i) - temp(k,j,im))**2) + &
    !                     0.5*dy**2*((temp(k,jp,i) - temp(k,j,i))**2 + (temp(k,j,i) - temp(k,jm,i))**2)
    !         end do
    !         k = 1
    !         kp = 2
    !         chiT2(k) = chiT2(k) + &
    !                 0.5*((temp(kp,j,i) - temp(k,j,i))**2*udx3c(kp)**2 + &
    !                     TfixS*(temp(k,j,i) - tempbp(1,j,i))**2/xm(k)**2)
    !         do k=2,nxm-1
    !             kp = k + 1
    !             km = k - 1
    !             chiT2(k) = chiT2(k) + &     ! (dT/dx)^2
    !                     0.5*((temp(kp,j,i) - temp(k,j,i))**2*udx3c(kp)**2 + &
    !                          (temp(k,j,i) - temp(km,j,i))**2*udx3c(k)**2)
    !         end do
    !         k = nxm
    !         km = nxm - 1
    !         chiT2(k) = chiT2(k) + &
    !                 0.5*(TfixN*(temptp(1,j,i) - temp(k,j,i))**2/(alx3 - xm(k))**2 + &
    !                     (temp(k,j,i) - temp(km,j,i))**2*udx3c(k)**2)
    !     end do
    ! end do

    ! call MpiSumReal1D(chiT2,nxm)

    do i=xstart(3),xend(3)
        ip = i + 1
        im = i - 1
        do j=xstart(2),xend(2)
            jp = j + 1
            jm = j - 1
            do k=1,nxm
                kp = k + 1
                km = k - 1
                epsilon(k) = epsilon(k) + ((vz(k,j,ip)-vz(k,j,i))*dz)**2 &   ! (dw/dz)^2
                                & + 0.125*dz**2*((vy(k,jp,ip) - vy(k,jp,im))**2 + (vy(k,j,ip) - vy(k,j,im))**2) & ! (dv/dz)^2
                                & + 0.125*dz**2*((vx(kp,j,ip) - vx(kp,j,im))**2 + (vx(k,j,ip) - vx(k,j,im))**2) & ! (du/dz)^2
                                & + 0.125*dy**2*((vz(k,jp,ip) - vz(k,jm,ip))**2 + (vz(k,jp,i) - vz(k,jm,i))**2) & ! (dw/dy)^2
                                & + ((vy(k,j+1,i)-vy(k,j,i))*dy)**2 &   ! (dv/dy)^2
                                & + 0.125*dy**2*((vx(kp,jp,i) - vx(k,jm,i))**2 + (vx(k,jp,i) - vx(k,jm,i))**2) & ! (du/dy)^2
                                & + ((vx(k+1,j,i)-vx(k,j,i))*udx3m(k))**2   ! (du/dx)^2
            end do
            tdx = 0.125*(dx/g3rm(1))**2
            epsilon(1) = epsilon(1) + tdx*((vz(2,j,ip) - vz(1,j,ip) + inslws*2.0*(vz(1,j,ip) - 0.0))**2 &
                                   + (vz(2,j,i) - vz(1,j,i) + inslws*2.0*(vz(1,j,i) - 0.0))**2) &
                            & + tdx*((vy(2,jp,i) - vy(1,jp,i) + inslws*2.0*(vy(1,jp,i) - xminusU))**2 &
                                   + (vy(2,j,i) - vy(1,j,i) + inslws*2.0*(vy(1,j,i) - xminusU))**2)
            do k=2,nxm-1
                kp = k + 1
                km = k - 1
                tdx = 0.125*(dx/g3rm(k))**2
                epsilon(k) = epsilon(k) + tdx*((vz(kp,j,ip) - vz(km,j,ip))**2 + (vz(kp,j,i) - vz(km,j,i))**2) &    ! (dw/dx)^2
                                & + tdx*((vy(kp,jp,i) - vy(km,jp,i))**2 + (vy(kp,j,i) - vy(km,j,i))**2)      ! (dv/dx)^2
            end do
            tdx = 0.125*(dx/g3rm(nxm))**2
            epsilon(nxm) = epsilon(nxm) + tdx*((vz(nxm,j,ip) - vz(nxm-1,j,ip) + inslwn*2.0*(0.0 - vz(nxm,j,ip)))**2 &
                                       + (vz(nxm,j,i) - vz(nxm-1,j,i) + inslwn*2.0*(0.0 - vz(nxm,j,i)))**2) &
                                & + tdx*((vy(nxm,jp,i) - vy(nxm-1,jp,i) + inslwn*2.0*(xplusU - vy(nxm,jp,i)))**2 &
                                       + (vy(nxm,j,i) - vy(nxm-1,j,i) + inslwn*2.0*(xplusU - vy(nxm,j,i)))**2)
        end do
    end do

    call MpiSumReal1D(epsilon,nxm)

    do k=1,nxm
        chiT(k) = chiT(k)*inym*inzm/pect
        ! chiT2(k) = chiT2(k)*inym*inzm/pect
        epsilon(k) = epsilon(k)*inym*inzm/ren
    end do

    write(nstat,"(i5.5)")nint(time/tout)
    filename = trim("outputdir/means.h5")
    
    inquire(file=filename,exist=fexist)
    if (.not.fexist) then
        if (ismaster) then
            call HdfCreateMeansFile(filename)
            if (salinity) call CreateSalinityH5Groups(filename)
            if (phasefield) call CreatePhaseH5Groups(filename)
            if (moist) call CreateMoistH5Groups(filename)
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

        dsetname = trim("vxvy/"//nstat)
        call HdfSerialWriteReal1D(dsetname,filename,vxvy,nxm)

        dsetname = trim("vxvz/"//nstat)
        call HdfSerialWriteReal1D(dsetname,filename,vxvz,nxm)

        dsetname = trim("vyvz/"//nstat)
        call HdfSerialWriteReal1D(dsetname,filename,vyvz,nxm)

        dsetname = trim("chiT/"//nstat)
        call HdfSerialWriteReal1D(dsetname,filename,chiT,nxm)
        
        ! dsetname = trim("chiT2/"//nstat)
        ! call HdfSerialWriteReal1D(dsetname,filename,chiT2,nxm)
        
        dsetname = trim("epsilon/"//nstat)
        call HdfSerialWriteReal1D(dsetname,filename,epsilon,nxm)

    end if

    if (salinity) call CalcSalStats

    if (phasefield) call CalcPhiStats

    if (moist) call CalcMoistStats

    call MpiBarrier

    tii(2) = MPI_WTIME()
    if (ismaster) write(*,"(a,f8.3,a)") "Profile save duration: ",tii(2)-tii(1),"s"

end subroutine


! Subroutine for creating HDF5 file with structure to contain mean profiles

subroutine HdfCreateMeansFile(filename)
    use hdf5
    use param
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
    call h5gcreate_f(file_id,"vxvy",group_id,hdf_error)
    call h5gclose_f(group_id,hdf_error)
    call h5gcreate_f(file_id,"vxvz",group_id,hdf_error)
    call h5gclose_f(group_id,hdf_error)
    call h5gcreate_f(file_id,"vyvz",group_id,hdf_error)
    call h5gclose_f(group_id,hdf_error)
    call h5gcreate_f(file_id,"chiT",group_id,hdf_error)
    call h5gclose_f(group_id,hdf_error)
    ! call h5gcreate_f(file_id,"chiT2",group_id,hdf_error)
    ! call h5gclose_f(group_id,hdf_error)
    call h5gcreate_f(file_id,"epsilon",group_id,hdf_error)
    call h5gclose_f(group_id,hdf_error)

    call h5fclose_f(file_id, hdf_error)
end subroutine HdfCreateMeansFile