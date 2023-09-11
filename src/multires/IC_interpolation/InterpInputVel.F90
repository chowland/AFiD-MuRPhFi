subroutine InterpInputVel

    use param
    use input_grids
    use local_arrays, only: vx,vy,vz,temp
    use mgrd_arrays
    use mpih
    use decomp_2d
    use AuxiliaryRoutines
    use HermiteInterpolations, only: interpolate_xyz_old_to_new
    implicit none

    integer  :: ic,jc,kc, ip,jp,kp, icr,jcr,kcr
    integer  :: jc0,jcr0, ic0,icr0
    integer  :: comm_col,comm_row,ierror

    real,dimension(4,4) :: qv2
    real,dimension(4) :: qv1

    real,allocatable,dimension(:,:) :: vyxzc,vyxzr, vzxyc,vzxyr
    real, allocatable, dimension(:,:,:) :: tempo, vxo, vyo, vzo
    real, allocatable, dimension(:,:,:) :: tpdvo

    real  :: lxa, lya, lza, dyo, dzo, Tup, Tlo

      !-- Allocate temporary arrays for velocities and gradients
    real,allocatable,dimension(:,:) :: dvyloc,dvybot, dvzloc,dvzbot

    call AllocateReal3DArray(tpdvo,-1,nxo+1, &
            xstarto(2)-2,xendo(2)+2,xstarto(3)-2,xendo(3)+2)

    tpdvo(:,:,:) = 0.d0     ! Temporary gradient array - old
    tpdv(:,:,:) = 0.d0      ! Temporary gradient array - new
    vx(:,:,:) = 0.d0
    vy(:,:,:) = 0.d0
    vz(:,:,:) = 0.d0

    dyo = 1.d0/(yco(2) - yco(1))
    dzo = 1.d0/(zco(2) - zco(1))

!=========================================================
!     Interpolation of vx

    ! Allocate and read in old vx
    call AllocateReal3DArray(vxo, 0, nxo+1, xstarto(2)-lvlhalo,xendo(2)+lvlhalo,xstarto(3)-lvlhalo,xendo(3)+lvlhalo)

    call HdfReadContinua(nzo, nyo, nxo, xstarto(2), xendo(2), xstarto(3), xendo(3), 1, &
            vxo(1:nxo, xstarto(2)-lvlhalo:xendo(2)+lvlhalo,xstarto(3)-lvlhalo:xendo(3)+lvlhalo))

    call update_halo(vxo, lvlhalo)


    !-- Interpolate dvx/dx to cell center after transpose
    do ic=xstarto(3),xendo(3)
        ip=ic+1
        do jc=xstarto(2),xendo(2)
            jp=jc+1
    
            !-- Interior points
            do kc=1,nxmo
                kp=kc+1
                tpdvo(kc,jc,ic)=(vxo(kp,jc,ic)-vxo(kc,jc,ic))/(xco(kp)-xco(kc))
            enddo
            ! CJH du/dx=0 on boundaries
            !-- Boundary points, enforce continuity
            tpdvo(0,jc,ic) = -tpdvo(1,jc,ic)
            tpdvo(nxo,jc,ic) = -tpdvo(nxmo,jc,ic)
    
        enddo
    enddo
    
    call update_halo(tpdvo,2)

    call interpolate_xyz_old_to_new(tpdvo, tpdv(1:nxm,:,:))

    !-- Integrate vx using interpolated gradients
    vx(1,:,:)=0.d0
    do icr=xstart(3),xend(3)
        do jcr=xstart(2),xend(2)
            do kcr=1,nxm
                lxa = xc(kcr+1)-xc(kcr)
                vx(kcr+1,jcr,icr) = vx(kcr,jcr,icr)+tpdv(kcr,jcr,icr)*lxa
            enddo
        enddo
    enddo

    call DestroyReal3DArray(vxo)


!=========================================================
!     Interpolation of vy

    ! Allocate and read in old vy
    call AllocateReal3DArray(vyo, 0, nxo+1, xstarto(2)-lvlhalo, xendo(2)+lvlhalo, xstarto(3)-lvlhalo, xendo(3)+lvlhalo)

    call HdfReadContinua(nzo, nyo, nxo, xstarto(2), xendo(2), xstarto(3), xendo(3), 2, &
            vyo(1:nxo, xstarto(2)-lvlhalo:xendo(2)+lvlhalo, xstarto(3)-lvlhalo:xendo(3)+lvlhalo))

    call update_halo(vyo, lvlhalo)

    tpdvo(:,:,:) = 0.d0
    tpdv(:,:,:) = 0.d0

    !-- Interpolate dvy/dy to cell center after transpose
    do ic=xstarto(3),xendo(3)
        do jc=xstarto(2),xendo(2)
            jp=jc+1

            !-- Interior points
            do kc=1,nxmo
                tpdvo(kc,jc,ic)=(vyo(kc,jp,ic)-vyo(kc,jc,ic))*dyo
            enddo

            !-- Boundary points, enforce zero velocity
            !CJH e.g. v=0 on x=0 => dv/dy=0 on x=0
            tpdvo(0,jc,ic) = (1-2*inslwS)*tpdvo(1,jc,ic)
            tpdvo(nxo,jc,ic) = (1-2*inslwN)*tpdvo(nxmo,jc,ic)

        enddo
    enddo

    call update_halo(tpdvo,2)

    call interpolate_xyz_old_to_new(tpdvo, tpdv(1:nxm,:,:))

    !-- Interpolate vyr at an arbitrary x-z plane | tpvr | 1st guess
    call AllocateReal2DArray(vyxzc,-1,nxmo+2,xstarto(3)-2,xendo(3)+2)
    call AllocateReal2DArray(vyxzr,1,nxm,xstart(3),xend(3))
    jc0 = 1 !global index
    jcr0 = 1
    vyxzc(:,:)=0.d0
    if (jc0.ge.xstarto(2).and.jc0.le.xendo(2)) then
        do ic=xstarto(3)-2,xendo(3)+2
            do kc=1,nxmo
                vyxzc(kc,ic)=vyo(kc,jc0,ic) !CS Halo updates can be optimised. Otherwise lvlhalo=2 required
            enddo
            ! x boundaries
            vyxzc(0,ic) = (1-2*inslwS)*vyxzc(1,ic)
            vyxzc(nxo,ic) = (1-2*inslwN)*vyxzc(nxmo,ic)
        enddo

        do ic=xstarto(3)-1,xendo(3)
            do kc=0,nxmo
    
                qv2=vyxzc(kc-1:kc+2,ic-1:ic+2)
    
                do icr=max(krangs(ic),xstart(3)),min(krangs(ic+1)-1,xend(3))
                    qv1(:) = qv2(:,1)*czvy(1,icr)+qv2(:,2)*czvy(2,icr) &
                            +qv2(:,3)*czvy(3,icr)+qv2(:,4)*czvy(4,icr)
                    do kcr=max(irangs(kc),1),min(irangs(kc+1)-1,nxm)
                        vy(kcr,jcr0,icr) = sum(qv1*cxvy(:,kcr))
                    enddo
                enddo
    
            enddo
        enddo
    endif

    call MPI_CART_SUB(DECOMP_2D_COMM_CART_X,(/.true.,.false./),comm_col,ierror)
    !call MPI_Comm_rank(comm, row_id, ierror)

    lya = 1./dy
    !-- Integrate along each y-line 
    do icr=xstart(3),xend(3)
        do jcr=xstart(2),xend(2)
            do kcr=1,nxm
                vy(kcr,jcr+1,icr)=vy(kcr,jcr,icr)+tpdv(kcr,jcr,icr)*lya
            enddo
        enddo
    enddo

    call AllocateReal2DArray(dvyloc,1,nxm,xstart(3),xend(3)) !Local
    call AllocateReal2DArray(dvybot,1,nxm,xstart(3),xend(3)) !Scan

    dvyloc=vy(1:nxm,xend(2)+1,xstart(3):xend(3))
    dvybot=0.d0

    !-- Scan velocity in y direction
    call MPI_SCAN(dvyloc,dvybot,nxm*xsize(3),MDP,MPI_SUM,comm_col,ierror)
    do icr=xstart(3),xend(3)
        do jcr=xstart(2),xend(2)+1
            do kcr=1,nxm
                vy(kcr,jcr,icr)=vy(kcr,jcr,icr)+dvybot(kcr,icr)-dvyloc(kcr,icr)
            enddo
        enddo
    enddo
    call MPI_Comm_free(comm_col,ierror) !CS Probably can be avoided by creating comm once
    ! call MPI_BARRIER(comm,ierr)

    call DestroyReal3DArray(vyo)
  
!=========================================================
!     Interpolation of vz

    ! Allocate and read in old vz
    call AllocateReal3DArray(vzo, 0, nxo+1, xstarto(2)-lvlhalo, xendo(2)+lvlhalo, xstarto(3)-lvlhalo, xendo(3)+lvlhalo)

    call HdfReadContinua(nzo, nyo, nxo, xstarto(2), xendo(2), xstarto(3), xendo(3), 3, &
            vzo(1:nxo, xstarto(2)-lvlhalo:xendo(2)+lvlhalo, xstarto(3)-lvlhalo:xendo(3)+lvlhalo))

    call update_halo(vzo, lvlhalo)

    tpdvo(:,:,:) = 0.d0
    tpdv(:,:,:) = 0.d0

      !-- Interpolate dvz/dz to cell center after transpose
    do ic=xstarto(3),xendo(3)
        ip=ic+1
        do jc=xstarto(2),xendo(2)
 
            !-- Interior points
            do kc=1,nxmo
                tpdvo(kc,jc,ic)=(vzo(kc,jc,ip)-vzo(kc,jc,ic))*dzo
            enddo
    
            !-- Boundary points, enforce zero velocity
            tpdvo(0,jc,ic) = -tpdvo(1,jc,ic)
            tpdvo(nxo,jc,ic) = -tpdvo(nxmo,jc,ic)

        enddo
    enddo
 
    call update_halo(tpdvo,2)

    call interpolate_xyz_old_to_new(tpdvo, tpdv(1:nxm,:,:))
    
    !-- Interpolate vzr at an arbitrary x-y plane | tpvr | 1st guess
    call AllocateReal2DArray(vzxyc,-1,nxmo+2,xstarto(2)-2,xendo(2)+2)
    call AllocateReal2DArray(vzxyr,1,nxm,xstart(2),xend(2))
    ic0 = 1 !global index
    icr0 = 1
    vzxyc(:,:)=0.d0
    if (ic0.ge.xstarto(3).and.ic0.le.xendo(3)) then
        do jc=xstarto(2)-2,xendo(2)+2
            do kc=1,nxmo
                vzxyc(kc,jc)=vzo(kc,jc,ic0) !CS Halo updates can be optimised. Otherwise lvlhalo=2 required
            enddo
            ! Boundary points
            vzxyc(0,jc) = -vzxyc(1,jc)
            vzxyc(nxo,jc) = -vzxyc(nxmo,jc)
        enddo

        do jc=xstarto(2)-1,xendo(2)
            do kc=0,nxmo

                qv2=vzxyc(kc-1:kc+2,jc-1:jc+2)

                do jcr=max(jrangs(jc),xstart(2)),min(jrangs(jc+1)-1,xend(2))
                    qv1(:) = qv2(:,1)*cyvz(1,jcr)+qv2(:,2)*cyvz(2,jcr) &
                            +qv2(:,3)*cyvz(3,jcr)+qv2(:,4)*cyvz(4,jcr)
                    do kcr=max(irangs(kc),1),min(irangs(kc+1)-1,nxm)
                        vz(kcr,jcr,icr0) = sum(qv1*cxvz(:,kcr))
                    enddo
                enddo

            enddo
        enddo
    endif

    call MPI_CART_SUB(DECOMP_2D_COMM_CART_X,(/.false.,.true./),comm_row,ierror)
    !call MPI_Comm_rank(comm, row_id, ierror)

    lza = 1./dz
    !-- Integrate along each y-line 
    do icr=xstart(3),xend(3)
        do jcr=xstart(2),xend(2)
            do kcr=1,nxm
                vz(kcr,jcr,icr+1)=vz(kcr,jcr,icr)+tpdv(kcr,jcr,icr)*lza
            enddo
        enddo
    enddo

    call AllocateReal2DArray(dvzloc,1,nxm,xstart(2),xend(2)) !Local
    call AllocateReal2DArray(dvzbot,1,nxm,xstart(2),xend(2)) !Scan

    dvzloc=vz(1:nxm,xstart(2):xend(2),xend(3)+1)
    dvzbot=0.d0

    !-- Scan velocity in z direction
    call MPI_SCAN(dvzloc,dvzbot,nxm*xsize(2),MDP,MPI_SUM,comm_row,ierror)
    do icr=xstart(3),xend(3)+1
        do jcr=xstart(2),xend(2)
            do kcr=1,nxm
                vz(kcr,jcr,icr)=vz(kcr,jcr,icr)+dvzbot(kcr,jcr)-dvzloc(kcr,jcr)
            enddo
        enddo
    enddo
    call MPI_Comm_free(comm_row,ierror) !CS Probably can be avoided by creating comm once
    !call MPI_BARRIER(comm,ierr)

    call DestroyReal3DArray(vzo)

    call DestroyReal3DArray(tpdvo)

!=========================================================
!     Interpolation of T

    ! Allocate and read in old T
    call AllocateReal3DArray(tempo, -1, nxo+1, xstarto(2)-lvlhalo, xendo(2)+lvlhalo, xstarto(3)-lvlhalo, xendo(3)+lvlhalo)

    tempo(:,:,:) = 0.d0

    call HdfReadContinua(nzo, nyo, nxo, xstarto(2), xendo(2), xstarto(3), xendo(3), 4, &
            tempo(1:nxo, xstarto(2)-lvlhalo:xendo(2)+lvlhalo, xstarto(3)-lvlhalo:xendo(3)+lvlhalo))

    ! Temperature BCs
    if (RAYT>0) then
        if (inslwN==0) then ! Single heated wall
            Tup = 0.0
            Tlo = 1.0
        else
            Tup = -0.5
            Tlo = 0.5
        end if
    else
        Tup = 0.5
        Tlo = -0.5
    end if
    if (phasefield) then
        if (RAYT>0) then
            Tup = 0.0
            Tlo = 1.0
        else
            Tup = 1.0
            Tlo = 0.0
        end if
    end if

    do ic=xstarto(3),xendo(3)
        do jc=xstarto(2),xendo(2)
            if (TfixS==1) then
                tempo(0,jc,ic) = 2.0*Tlo - tempo(1,jc,ic)
            else
                tempo(0,jc,ic) = tempo(1,jc,ic)
            end if
            if (TfixN==1) then
                tempo(nxo,jc,ic) = 2.0*Tup - tempo(nxmo,jc,ic)
            else
                tempo(nxo,jc,ic) = tempo(nxmo,jc,ic)
            end if
        end do
    end do

    call update_halo(tempo, lvlhalo)

    call interpolate_xyz_old_to_new(tempo, temp(1:nxm,:,:))

    call DestroyReal3DArray(tempo)

    return
end subroutine InterpInputVel