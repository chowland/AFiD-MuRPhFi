!> Interpolate the humidity field provided by `continua_humid.h5`
!! from a previous grid onto the current grid specified by `bou.in`
subroutine InterpInputHum
    use input_grids
    use afid_sides
    use afid_moisture
    use HermiteInterpolations, only: interpolate_xyz_old_to_new
    real, allocatable, dimension(:,:,:) :: qold !< Humidity on the old grid
    real :: qup     !< Upper boundary condition for humidity
    real :: qlo     !< Lower boundary condition for humidity
    integer :: ic, jc
    real :: dyo, dzo
    
    ! Read field from h5 file and store in qold
    call AllocateReal3DArray(qold, -1, nxo+1, &
            xstarto(2)-lvlhalo, xendo(2)+lvlhalo, &
            xstarto(3)-lvlhalo, xendo(3)+lvlhalo)
    qold(:,:,:) = 0.d0
    call HdfReadContinua(nzo, nyo, nxo, &
        xstarto(2), xendo(2), xstarto(3), xendo(3), 8, &
        qold(1:nxo,xstarto(2)-lvlhalo:xendo(2)+lvlhalo,xstarto(3)-lvlhalo:xendo(3)+lvlhalo))
    
    ! Boundary conditions (called after SetHumidityBCs)
    ! ASSUMES HUMTP, HUMBP ARE CONSTANT)
    qup = humtp(1,xstart(2),xstart(3))
    qlo = humbp(1,xstart(2),xstart(3))

    ! Apply boundary conditions in x-halo cells
    do ic=xstarto(3),xendo(3)
        do jc=xstarto(2),xendo(2)
            if (qfixS==1) then
                qold(0,jc,ic) = 2.0*qlo - qold(1,jc,ic)
            else
                qold(0,jc,ic) = qold(1,jc,ic)
            end if
            if (qfixN==1) then
                qold(nxo,jc,ic) = 2.0*qup - qold(nxmo,jc,ic)
            else
                qold(nxo,jc,ic) = qold(nxmo,jc,ic)
            end if
        end do
    end do
    call update_halo(qold, lvlhalo)

    !! Add side boundary conditions if using
    if (sidewall) then
        dyo = yco(2) - yco(1)
        dzo = zco(2) - zco(1)
        if (xstarto(2)==1) then
            call ApplyBC(qold, bc_humid_y_fix_lo, bc_humid_y_val_lo, 'y', 'L', 2, dyo)
        end if
        if (xstarto(3)==1) then
            call ApplyBC(qold, bc_humid_z_fix_lo, bc_humid_z_val_lo, 'z', 'L', 2, dzo)
        end if
        if (xendo(2)==nymo) then
            call ApplyBC(qold, bc_humid_y_fix_up, bc_humid_y_val_up, 'y', 'U', 2, dyo)
        end if
        if (xendo(3)==nzmo) then
            call ApplyBC(qold, bc_humid_z_fix_up, bc_humid_z_val_up, 'z', 'U', 2, dzo)
        end if
    end if

    call interpolate_xyz_old_to_new(qold, humid(1:nxm,:,:))
    
    call DestroyReal3DArray(qold)

end subroutine InterpInputHum