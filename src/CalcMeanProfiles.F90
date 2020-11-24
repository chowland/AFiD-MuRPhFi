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
                chiT(k) = chiT(k) + ((temp(k,j,i+1)-temp(k,j,i-1))*0.5*dz)**2*inym*inzm
                chiT(k) = chiT(k) + ((temp(k,j+1,i)-temp(k,j-1,i))*0.5*dy)**2*inym*inzm
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
            call HdfCreateBlankFile(filename)
        end if
    end if

    if (ismaster) then
        dsetname = trim("Tbar"//nstat)
        call HdfSerialWriteReal1D(dsetname,filename,Tbar,nxm)

        dsetname = trim("vybar"//nstat)
        call HdfSerialWriteReal1D(dsetname,filename,vybar,nxm)

        dsetname = trim("vzbar"//nstat)
        call HdfSerialWriteReal1D(dsetname,filename,vzbar,nxm)

        dsetname = trim("vxT"//nstat)
        call HdfSerialWriteReal1D(dsetname,filename,vxT,nxm)

        dsetname = trim("vyT"//nstat)
        call HdfSerialWriteReal1D(dsetname,filename,vyT,nxm)

        dsetname = trim("vzT"//nstat)
        call HdfSerialWriteReal1D(dsetname,filename,vzT,nxm)

        dsetname = trim("Trms"//nstat)
        call HdfSerialWriteReal1D(dsetname,filename,Trms,nxm)

        dsetname = trim("vxrms"//nstat)
        call HdfSerialWriteReal1D(dsetname,filename,vxrms,nxm)
        
        dsetname = trim("vyrms"//nstat)
        call HdfSerialWriteReal1D(dsetname,filename,vyrms,nxm)

        dsetname = trim("vzrms"//nstat)
        call HdfSerialWriteReal1D(dsetname,filename,vzrms,nxm)

        dsetname = trim("chiT"//nstat)
        call HdfSerialWriteReal1D(dsetname,filename,chiT,nxm)
        
        dsetname = trim("Sbar"//nstat)
        call HdfSerialWriteReal1D(dsetname,filename,Sbar,nxmr)

        dsetname = trim("vxS"//nstat)
        call HdfSerialWriteReal1D(dsetname,filename,vxS,nxmr)

        dsetname = trim("vyS"//nstat)
        call HdfSerialWriteReal1D(dsetname,filename,vyS,nxmr)

        dsetname = trim("vzS"//nstat)
        call HdfSerialWriteReal1D(dsetname,filename,vzS,nxmr)

        dsetname = trim("Srms"//nstat)
        call HdfSerialWriteReal1D(dsetname,filename,Srms,nxmr)

        dsetname = trim("chiS"//nstat)
        call HdfSerialWriteReal1D(dsetname,filename,chiS,nxmr)
    end if

    call MpiBarrier

    tii(2) = MPI_WTIME()
    if (ismaster) write(*,"(a,f8.3,a)") "Profile save duration: ",tii(2)-tii(1),"s"

end subroutine