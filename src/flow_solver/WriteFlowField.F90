!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: WriteFlowField.F90                             !
!    CONTAINS: subroutine WriteFlowField                  !
!                                                         ! 
!    PURPOSE: Write down the full flow snapshot for       !
!     restarting the simulation at a later date           !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine WriteFlowField(final)
    use param
    use local_arrays, only: vz,vy,vx,temp
    use mgrd_arrays, only: sal,phi
    implicit none
    logical, intent(in) :: final
    character(30) :: filnam1,dsetname,basename
    character(5)  :: frame
    
    basename='outputdir/continua'
    if (.not.final) then
        write(frame,"(i5.5)")nint(time/save_3D)
        basename='outputdir/fields/'//frame
    end if
    filnam1 = trim(basename)//'_temp.h5'
    call HdfWriteRealHalo3D(filnam1,temp)
    filnam1 = trim(basename)//'_vx.h5'
    call HdfWriteRealHalo3D(filnam1,vx)
    filnam1 = trim(basename)//'_vy.h5'
    call HdfWriteRealHalo3D(filnam1,vy)
    filnam1 = trim(basename)//'_vz.h5'
    call HdfWriteRealHalo3D(filnam1,vz)
    if (salinity) then
        filnam1 = trim(basename)//'_sal.h5'
        call HdfWriteRealHalo3DR(filnam1,sal)
    end if
    if (phasefield) then
        filnam1 = trim(basename)//'_phi.h5'
        call HdfWriteRealHalo3DR(filnam1,phi)
    end if
    
    if (ismaster .and. final) then !EP only write once
        filnam1 = trim('outputdir/continua_master.h5')
        call HdfCreateBlankFile(filnam1)
        
        dsetname = trim('nx')
        call HdfSerialWriteIntScalar(dsetname,filnam1,nx)
        dsetname = trim('ny')
        call HdfSerialWriteIntScalar(dsetname,filnam1,ny)
        dsetname = trim('nz')
        call HdfSerialWriteIntScalar(dsetname,filnam1,nz)
        dsetname = trim('ylen')
        call HdfSerialWriteRealScalar(dsetname,filnam1,ylen)
        dsetname = trim('zlen')
        call HdfSerialWriteRealScalar(dsetname,filnam1,zlen)
        dsetname = trim('time')
        call HdfSerialWriteRealScalar(dsetname,filnam1,time)
        dsetname = trim('istr3')
        call HdfSerialWriteIntScalar(dsetname,filnam1,istr3)
        dsetname = trim('str3')
        call HdfSerialWriteRealScalar(dsetname,filnam1,str3)
        !CS Mgrd
        if (multires) then
            dsetname = trim('nxr')
            call HdfSerialWriteIntScalar(dsetname,filnam1,nxr)
            dsetname = trim('nyr')
            call HdfSerialWriteIntScalar(dsetname,filnam1,nyr)
            dsetname = trim('nzr')
            call HdfSerialWriteIntScalar(dsetname,filnam1,nzr)
            dsetname = trim('istr3r')
            call HdfSerialWriteIntScalar(dsetname,filnam1,istr3r)
        end if
        
    endif
    
end subroutine WriteFlowField
