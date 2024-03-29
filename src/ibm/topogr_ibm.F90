!
!     This routine finds the indices in the computational grid
!     close to the physical location of the body.
!
subroutine topogr
    use param
    use decomp_2d, only: xstart,xend,xstartr,xendr
    use ibm_param
    use afid_salinity, only: RayS
    use mpih
    implicit none
    integer :: i,j,k,l, nc, Npart, ncz
    integer :: km,kp
    
    real    :: xe, xem, xep
    real    :: ye, yem
    real    :: ze, zem
    real    :: delta1x, delta2x, r2, Lhex, radius, porosity
    real    :: solid_temp, rp, tp, amp, rx, ry, rz
    real,allocatable :: xpart(:), ypart(:), zpart(:)
    integer :: ibmask(1:nx,xstart(2):xend(2),xstart(3):xend(3))

    !! Variables for writing details of solid centres to file
    character(len=30) :: dsetname, filename

    logical :: fexist


    allocate(plth1(1:nym,1:nzm))
    allocate(plth2(1:nym,1:nzm))
    plth1=1.d0
    plth2=1.d0
    
    npunx=0
    npuny=0
    npunz=0
    !
    !     IDENTIFICATION OF THE GRID POINTS IN THE BODY
    !     (BOUNDARY + INNER PART)
    !
    !     X_3 vertical direction
    !     X_2 radial direction
    !     X_1 azimuthal direction
    !
    
    allocate(ibmaskx(1:nx,xstart(2):xend(2),xstart(3):xend(3)))
    allocate(ibmasky(1:nx,xstart(2):xend(2),xstart(3):xend(3)))
    allocate(ibmaskz(1:nx,xstart(2):xend(2),xstart(3):xend(3)))
    allocate(ibmaskt(1:nx,xstart(2):xend(2),xstart(3):xend(3)))
    if (salinity) then
        ! allocate(solidr(1:nxmr,xstartr(2)-1:xendr(2)+1,xstartr(3)-1:xendr(3)+1))
        allocate(solidr(-1:nxmr+1,xstartr(2)-2:xendr(2)+2,xstartr(3)-2:xendr(3)+2))
        allocate(ibmaskr(1:nxmr,xstartr(2):xendr(2),xstartr(3):xendr(3)))
    end if

    ! Build array of solid circle positions
    if (solidtype==1) then
        porosity = 0.37
        ! Interpret RayT input as target pore-scale Rayleigh number
        nc = nint((2.0*ylen*(1.0 - porosity)/3.0/pi)**0.5 * (RayS/RayT)**(1.0/3.0))
        if (ismaster) write(*,*) "Number of lattice columns: ",nc
        Lhex = 1.0/nc/ylen
        radius = (ylen*(1.0 - porosity)/6.0/pi)**0.5/nc
        amp = radius*(sqrt(pi/2.0/3.0**0.5/(1.0 - porosity)) - 1)
        Npart = (nc + 1)*(3*nc + 1) + 3*nc*nc
        allocate(xpart(1:Npart))
        allocate(ypart(1:Npart))

        ! Compute solid centres on root process
        if (ismaster) then
            filename = trim("outputdir/solid_centres.h5")
            inquire(file=filename, exist=fexist)
            ! If locations already prescribed, read them from file
            if (fexist) then
                write(*,*) "Reading solid matrix from file."
                dsetname = trim("xpart")
                call HdfSerialReadReal1D(dsetname, filename, xpart, Npart)
                dsetname = trim("ypart")
                call HdfSerialReadReal1D(dsetname, filename, ypart, Npart)
            ! Otherwise, construct the porous matrix
            else
                i = 1
                call random_seed()
                do j=0,nc
                    xe = real(j)/real(nc)
                    do k=1,3*nc
                        ye = Lhex*real(k)
                        call random_number(rp)
                        call random_number(tp)
                        xpart(i) = xe + amp*rp*cos(2*pi*tp)
                        ypart(i) = ye + amp*rp*sin(2*pi*tp)
                        i = i + 1
                        if (k==3*nc) then
                            xpart(i) = xpart(i-1)
                            ypart(i) = ypart(i-1) - ylen
                            i = i + 1
                        end if
                    end do
                    if (j>0) then
                        xe = (real(j) - 0.5)/real(nc)*alx3
                        do k=1,3*nc
                            ye = Lhex*(real(k) - 0.5)
                            call random_number(rp)
                            call random_number(tp)
                            xpart(i) = xe + amp*rp*cos(2*pi*tp)
                            ypart(i) = ye + amp*rp*sin(2*pi*tp)
                            i = i + 1
                        end do
                    end if
                end do
                call HdfCreateBlankFile(filename)
                dsetname = trim("xpart")
                call HdfSerialWriteReal1D(dsetname, filename, xpart, Npart)
                dsetname = trim("ypart")
                call HdfSerialWriteReal1D(dsetname, filename, ypart, Npart)
            end if
        end if

        ! Broadcast centre positions to all processes
        call MPI_BCAST(xpart, Npart, MDP, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(ypart, Npart, MDP, 0, MPI_COMM_WORLD, ierr)
    elseif (solidtype==2) then
        Npart = 14
        allocate(xpart(1:Npart))
        allocate(ypart(1:Npart))
        i = 1
        do j=0,4
            ye = real(j)*ylen/4.0
            do k=0,1
                xe = real(k)
                xpart(i) = xe
                ypart(i) = ye
                i = i + 1
            end do
            if (j>0) then
                ye = (real(j) - 0.5)*ylen/4.0
                xe = 0.5
                xpart(i) = xe
                ypart(i) = ye
                i = i + 1
            end if
        end do

    elseif (solidtype==4) then
        porosity = 0.37
        nc = nint(1.0/sqrt(6.0)*(3*sqrt(2.0)*(1.0 - porosity)*RayS/RayT/pi)**(1.0/3.0))
        if (ismaster) write(*,*) "Number of lattice layers: ",nc
        ncz = nint(sqrt(3.0)*zlen/ylen*nc)
        rx = 0.5/sqrt(6.0)/nc
        ry = 0.5/sqrt(3.0)/nc*ylen
        rz = 0.5*zlen/ncz
        radius = (ylen*zlen/8.0/pi/nc**2/ncz*(1.0 - porosity))**(1.0/3.0)
        if (ismaster) write(*,*) "Bead radius: ",radius
        Npart = (ncz + 1)*((nc + 1)**2 + 2*nc**2) + ncz*3*nc*(nc + 1)
        allocate(xpart(1:Npart))
        allocate(ypart(1:Npart))
        allocate(zpart(1:Npart))

        ! Compute solid centres on root process
        if (ismaster) then
            ! Check if a file prescribing the solid centres already exists
            filename = trim("outputdir/solid_centres.h5")
            inquire(file=filename, exist=fexist)
            ! Read the solid centres if already existing
            if (fexist) then
                write(*,*) "Reading solid matrix from file."
                dsetname = trim("xpart")
                call HdfSerialReadReal1D(dsetname, filename, xpart, Npart)
                dsetname = trim("ypart")
                call HdfSerialReadReal1D(dsetname, filename, ypart, Npart)
                dsetname = trim("zpart")
                call HdfSerialReadReal1D(dsetname, filename, zpart, Npart)
            ! Construct the pore matrix otherwise
            else
                n = 1
                ! First layer
                do k=0,nc
                    xe = 2.0*k*sqrt(6.0)*rx
                    do i=0,ncz
                        do j=0,nc
                            ye = 2.0*sqrt(3.0)*j*ry
                            ze = 2.0*i*rz
                            xpart(n) = xe
                            ypart(n) = ye
                            zpart(n) = ze
                            n = n + 1
                        end do
                    end do
                    do i=1,ncz
                        do j=1,nc
                            ye = (2.0*j - 1.0)*sqrt(3.0)*ry
                            ze = (2.0*i - 1.0)*rz
                            xpart(n) = xe
                            ypart(n) = ye
                            zpart(n) = ze
                            n = n + 1
                        end do
                    end do
                end do
                ! Second layer
                do k=1,nc
                    xe = 2.0*(3.0*k - 2.0)*sqrt(6.0)/3.0*rx
                    do i=0,ncz
                        do j=1,nc
                            ye = sqrt(3.0)*(2.0*j - 2.0/3.0)*ry
                            ze = 2.0*i*rz
                            xpart(n) = xe
                            ypart(n) = ye
                            zpart(n) = ze
                            n = n + 1
                        end do
                    end do
                    do i=1,ncz
                        do j=0,nc
                            ye = sqrt(3.0)*(2.0*j + 1.0/3.0)*ry
                            ze = (2.0*i - 1.0)*rz
                            xpart(n) = xe
                            ypart(n) = ye
                            zpart(n) = ze
                            n = n + 1
                        end do
                    end do
                end do

                ! Third layer
                do k=1,nc
                    xe = 2.0*(3.0*k - 1.0)*sqrt(6.0)/3.0*rx
                    do i=0,ncz
                        do j=1,nc
                            ye = sqrt(3.0)*(2.0*j - 4.0/3.0)*ry
                            ze = 2.0*i*rz
                            xpart(n) = xe
                            ypart(n) = ye
                            zpart(n) = ze
                            n = n + 1
                        end do
                    end do
                    do i=1,ncz
                        do j=0,nc
                            ye = sqrt(3.0)*(2.0*j - 1.0/3.0)*ry
                            ze = (2.0*i - 1.0)*rz
                            xpart(n) = xe
                            ypart(n) = ye
                            zpart(n) = ze
                            n = n + 1
                        end do
                    end do
                end do

                call HdfCreateBlankFile(filename)
                dsetname = trim("xpart")
                call HdfSerialWriteReal1D(dsetname, filename, xpart, Npart)
                dsetname = trim("ypart")
                call HdfSerialWriteReal1D(dsetname, filename, ypart, Npart)
                dsetname = trim("zpart")
                call HdfSerialWriteReal1D(dsetname, filename, zpart, Npart)
            end if
        end if

        ! Broadcast centre positions to all processes
        call MPI_BCAST(xpart, Npart, MDP, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(ypart, Npart, MDP, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(zpart, Npart, MDP, 0, MPI_COMM_WORLD, ierr)

    end if

    do l = 1,3 !{ start do over the 3 velocity components
        n=0
        
        ibmask(:,:,:) = 2 ! Initialise mask as all liquid

        !     l = 1   Q_1 vel. component (VZ)
        !     l = 2   Q_2 vel. component (VY)
        !     l = 3   Q_3 vel. component (VX)
        !

        ! Construct surface height for scallop shape
        if (solidtype==2) then
            plth1(:,:) = 1.0
            do i=xstart(3),xend(3)
                do j=xstart(2),xend(2)
                    ye = ym(j)
                    ze = zm(i)
                    if (l==1) ze = zc(i)
                    if (l==2) ye = yc(j)
                    do k=1,Npart
                        yem = ypart(k)
                        zem = xpart(k)
                        plth1(j,i) = min(plth1(j,i), 0.2*(8.0/ylen)**2*((ye - yem)**2 + 3*(ze - zem)**2))
                    end do
                    ! plth2(j,i) = plth1(j,i)
                end do
            end do
        ! Streamwise riblets
        elseif (solidtype==3) then
            plth1(:,:) = 0.0
            do i=xstart(3),xend(3)
                do j=xstart(2),xend(2)
                    ye = ym(j)
                    ze = zm(i)
                    if (l==1) ze = zc(i)
                    if (l==2) ye = yc(j)
                    r2 = 0.05
                    plth1(j,i) = (2.0d0)*MIN((MOD(ze,r2)),(r2-MOD(ze,r2)))
                end do
            end do
        end if

        ! Track location of solid on each grid:
        do i=xstart(3),xend(3)
            ze = zm(i)
            if (l==1) ze = zc(i)
            do j=xstart(2),xend(2)
                ye = ym(j)
                if (l==2) ye = yc(j)
                do k=1,nxm
                    xe = xm(k)
                    if (l==3) xe = xc(k)
                    if (mod(solidtype,3) == 1) then
                        do nc=1,Npart
                            ! x-position of solid centre
                            xem = xpart(nc)
                            yem = ypart(nc)
                            if (solidtype==4) then
                                zem = zpart(nc)
                                r2 = (xe - xem)**2 + (ye - yem)**2 + (ze - zem)**2
                            else
                                r2 = (xe - xem)**2 + (ye - yem)**2
                            end if
                            if (r2<radius**2) then
                                ibmask(k,j,i) = 0
                            end if
                        end do
                    else
                        if (xe < plth1(j,i)) then
                            ibmask(k,j,i) = 0
                        end if
                    end if
                end do
            end do
        end do

        do i=xstart(3),xend(3)
            ze = zm(i)
            if (l==1) ze = zc(i)
            do j=xstart(2),xend(2)
                ye = ym(j)
                if (l==2) ye = yc(j)
                do k=1,nxm
                    km=kmv(k)
                    kp=kpv(k)
                    xe=xm(k)
                    xem=xm(km)
                    xep=xm(kp)
                    if(l.eq.3) then
                        xe=xc(k)
                        xem=xc(km)
                        xep=xc(kp)
                    end if

                    !    LOWER FLUID/PLATE BOUNDARY
                    !
                    if ((ibmask(k,j,i)==2) .and. (ibmask(km,j,i)==0)) then
                        n=n+1
                        delta1x=(xep-xe)
                        if (solidtype==1) then
                            delta2x = 1.0
                            do nc=1,Npart
                                amp = radius**2 - (ye - ypart(nc))**2
                                if (amp > 0) then
                                    xem = xpart(nc) + sqrt(amp)
                                    if (xe > xem) delta2x = min(delta2x, xe - xem)
                                end if
                            end do
                        elseif (solidtype==4) then
                            delta2x = 1.0
                            do nc=1,Npart
                                amp = radius**2 - (ye - ypart(nc))**2 - (ze - zpart(nc))**2
                                if (amp>0) then
                                    xem = xpart(nc) + sqrt(amp)
                                    if (xe > xem) delta2x = min(delta2x, xe - xem)
                                end if
                            end do
                        else
                            delta2x = xe - plth1(j,i)
                        end if
                        distb(n,l) = delta2x/(delta1x+delta2x)
                        ibmask(k,j,i) = 1
                            
                    !
                    !    UPPER FLUID/PLATE BOUNDARY
                    !
                    elseif ((ibmask(k,j,i)==2) .and. (ibmask(kp,j,i)==0)) then
                        n=n+1
                        delta1x=(xe-xem)
                        if (solidtype==1) then
                            delta2x = 1.0
                            do nc=1,Npart
                                amp = radius**2 - (ye - ypart(nc))**2
                                if (amp > 0) then
                                    xep = xpart(nc) - sqrt(amp)
                                    if (xep > xe) delta2x = min(delta2x, xep - xe)
                                end if
                            end do
                        elseif (solidtype==4) then
                            delta2x = 1.0
                            do nc=1,Npart
                                amp = radius**2 - (ye - ypart(nc))**2 - (ze - zpart(nc))**2
                                if (amp>0) then
                                    xep = xpart(nc) - sqrt(amp)
                                    if (xep > xe) delta2x = min(delta2x, xep - xe)
                                end if
                            end do
                        end if
                        ! delta2x=((alx3-plth2(j,i))-xe)
                        distb(n,l) = delta2x/(delta1x+delta2x)
                        ibmask(k,j,i) = -1
                        
                    end if
                                            
                end do
            end do
        end do
        ! iF(ANY(IsNaN(distb))) write(*,*) 'NaN in distb'
        
        if(l.eq.1) then
            if(n.gt.mpun)  &
            write(*,*) 'Dim max di indgeot e'' stata superata n=',n
            npunz= n
            if (ismaster) write(*,*) 'npunz: ',npunz
            !        write(6,332)npunz
            ! 332  format(5x,'For Q_1 N ='i7)
            ibmaskz(:,:,:) = ibmask
            distz(:) = distb(:,1)
        end if
        if(l.eq.2) then
            if(n.gt.mpun)   &
            write(*,*) 'Dim max di indgeor e'' stata superata n=',n
            npuny= n
            if (ismaster) write(*,*) 'npuny: ',npuny
            !        write(6,331)npuny
            ! 331  format(5x,'For Q_2 N ='i7)
            ibmasky(:,:,:) = ibmask
            disty(:) = distb(:,2)
        end if
        if(l.eq.3) then
            if(n.gt.mpun)  &
            write(*,*) 'Dim max di indgeoz e'' stata superata n=',n
            npunx= n
            if (ismaster) write(*,*) 'npunx: ',npunx
            !        write(6,330)npunx
            ! 330  format(5x,'For Q_3 N ='i7)
            ibmaskx(:,:,:) = ibmask
            distx(:) = distb(:,3)
        end if
    end do   !} end do over the 3 velocity components
                            
                            

    !
    !     INDICES FOR TEMPERATURE
    !

    ! Specify topography shape on temperature grid
    if (solidtype==2) then
        plth1(:,:) = 1.0
        do i=xstart(3),xend(3)
            do j=xstart(2),xend(2)
                ye = ym(j)
                ze = zm(i)
                do k=1,Npart
                    yem = ypart(k)
                    zem = xpart(k)
                    plth1(j,i) = min(plth1(j,i), 0.2*(8.0/ylen)**2*((ye - yem)**2 + 3*(ze - zem)**2))
                end do
                ! plth2(j,i) = plth1(j,i)
            end do
        end do
    ! Streamwise riblets
    elseif (solidtype==3) then
        plth1(:,:) = 0.0
        do i=xstart(3),xend(3)
            do j=xstart(2),xend(2)
                ye = ym(j)
                ze = zm(i)
                plth1(j,i) = (2.0d0)*MIN((MOD(ze,0.025d0)),(0.025d0-MOD(ze,0.025d0)))
            end do
        end do
    end if

    n=0
    ibmask(:,:,:) = 2
    !
    ! Track location of solid on temperature grid:
    do i=xstart(3),xend(3)
        ze = zm(i)
        do j=xstart(2),xend(2)
            ye = ym(j)
            do k=1,nxm
                xe = xm(k)
                if (mod(solidtype,3)==1) then
                    do nc=1,Npart
                        ! x-position of solid centre
                        xem = xpart(nc)
                        yem = ypart(nc)
                        if (solidtype==4) then
                            zem = zpart(nc)
                            r2 = (xe - xem)**2 + (ye - yem)**2 + (ze - zem)**2
                        else
                            r2 = (xe - xem)**2 + (ye - yem)**2
                        end if
                        if (r2<radius**2) then
                            ibmask(k,j,i) = 0
                        end if
                    end do
                else
                    if (xe < plth1(j,i)) then
                        ibmask(k,j,i) = 0
                    end if
                end if
            end do
        end do
    end do

    solid_temp = 1.0 ! Modify this to set fixed temperature value in solid
    do i=xstart(3),xend(3)
        ze = zm(i)
        do j=xstart(2),xend(2)
            ye = ym(j)
            do k=1,nxm
                km=kmv(k)
                kp=kpv(k)
                xe=xm(k)
                xem=xm(km)
                xep=xm(kp)
            !
            !    SOLID PART
            !           

                if ((ibmask(k,j,i)==2) .and. (ibmask(km,j,i)==0)) then
                    n=n+1
                    delta1x=(xep-xe)
                    if (solidtype==1) then
                        delta2x = 1.0
                        do nc=1,Npart
                            amp = radius**2 - (ye - ypart(nc))**2
                            if (amp > 0) then
                                xem = xpart(nc) + sqrt(amp)
                                if (xe > xem) delta2x = min(delta2x, xe - xem)
                            end if
                        end do
                    elseif (solidtype==4) then
                        delta2x = 1.0
                        do nc=1,Npart
                            amp = radius**2 - (ye - ypart(nc))**2 - (ze - zpart(nc))**2
                            if (amp > 0) then
                                xem = xpart(nc) + sqrt(amp)
                                if (xe > xem) delta2x = min(delta2x, xe - xem)
                            end if
                        end do
                    else
                        delta2x = xe - plth1(j,i)
                    end if
                    distbt(n) = delta2x/(delta1x+delta2x)
                    ibmask(k,j,i) = 1

                elseif ((ibmask(k,j,i)==2) .and. (ibmask(kp,j,i)==0)) then
                    n=n+1
                    delta1x=(xe-xem)
                    if (solidtype==1) then
                        delta2x = 1.0
                        do nc=1,Npart
                            amp = radius**2 - (ye - ypart(nc))**2
                            if (amp > 0) then    
                                xep = xpart(nc) - sqrt(amp)
                                if (xep > xe) delta2x = min(delta2x, xep - xe)
                            end if
                        end do
                    elseif (solidtype==4) then
                        delta2x = 1.0
                        do nc=1,Npart
                            amp = radius**2 - (ye - ypart(nc))**2 - (ze - zpart(nc))**2
                            if (amp > 0) then
                                xep = xpart(nc) - sqrt(amp)
                                if (xep > xe) delta2x = min(delta2x, xep - xe)
                            end if
                        end do
                    end if
                    ! delta2x=((alx3-plth2(j,i))-xe)
                    distbt(n) = delta2x/(delta1x+delta2x)
                    ibmask(k,j,i) = -1
                    ! temb(n) = solid_temp
                end if
                                
            end do
        end do
    end do
    if(n.gt.mpun) &
        write(*,*) 'Dim max di indgeote e'' stata superata n=',n
    npunte= n
    !        write(6,329)npunte
    ! 329  format(5x,'For Temperature N ='i7)
    ibmaskt(:,:,:) = ibmask
    if(allocated(plth1)) deallocate(plth1)
    if(allocated(plth2)) deallocate(plth2)
    
    if (salinity) then
        
        allocate(plth1(1:nymr,1:nzmr))
        allocate(plth2(1:nymr,1:nzmr))
        plth1=1.d0
        plth2=1.d0

        do i=1,nzmr
            do j=1,nymr
                do k=1,4
                    ye = ylen/8.0 + real(k - 1)*ylen/4.0
                    plth1(j,i) = min(plth1(j,i), 0.2*(8.0/ylen)**2*(ymr(j) - ye)**2)
                end do
                plth2(j,i) = plth1(j,i)
            end do
        end do
        !
        !     INDICES FOR SALINITY
        !
        n = 0
        ibmaskr(:,:,:) = 2
        solidr(:,:,:) = .false.

        do i=xstartr(3)-1,xendr(3)+1
            ze = zmr(i)
            do j=xstartr(2)-1,xendr(2)+1
                ye = ymr(j)
                do k=1,nxmr
                    xe = xmr(k)
                    if (mod(solidtype,3)==1) then
                        do nc=1,Npart
                            ! x-position of solid centre
                            xem = xpart(nc)
                            yem = ypart(nc)
                            if (solidtype==4) then
                                zem = zpart(nc)
                                r2 = (xe - xem)**2 + (ye - yem)**2 + (ze - zem)**2
                            else
                                r2 = (xe - xem)**2 + (ye - yem)**2
                            end if
                            solidr(k,j,i) = solidr(k,j,i) .or. (r2<radius**2)
                        end do
                    else
                        if (xe < plth1(j,i)) then
                            solidr(k,j,i) = .true.
                        end if
                    end if
                end do
            end do
        end do

        do i=xstartr(3),xendr(3)
            ze = zmr(i)
            do j=xstartr(2),xendr(2)
                ye = ymr(j)
                do k=1,nxmr
                    xe = xmr(k)
                    if (mod(solidtype,3)==1) then
                        do nc=1,Npart
                            ! x-position of solid centre
                            xem = xpart(nc)
                            yem = ypart(nc)
                            if (solidtype==4) then
                                zem = zpart(nc)
                                r2 = (xe - xem)**2 + (ye - yem)**2 + (ze - zem)**2
                            else
                                r2 = (xe - xem)**2 + (ye - yem)**2
                            end if
                            if (r2<radius**2) then
                                ibmaskr(k,j,i) = 0
                            end if
                        end do
                    else
                        if (xe < plth1(j,i)) then
                            ibmaskr(k,j,i) = 0
                        end if
                    end if
                end do
            end do
        end do

        ! klo = 1
        ! kup = 2
        do i=xstartr(3),xendr(3)
            do j=xstartr(2),xendr(2)
                do k=1,nxmr
                    km=max(1,k-1)
                    kp=min(k+1,nxmr)
                    xe=xmr(k)
                    xem=xmr(km)
                    xep=xmr(kp)
                    !
                    !    SOLID PART
                    !           

                    if ((ibmaskr(k,j,i)==0) .and. (ibmaskr(kp,j,i)==2)) then
                        n = n + 1
                        ibmaskr(k,j,i) = 1
                        distbr(n) = 1.0
                    elseif ((ibmaskr(k,j,i)==0) .and. (ibmaskr(km,j,i)==2)) then
                        n = n + 1
                        ibmaskr(k,j,i) = -1
                        distbr(n) = 1.0
                    end if
                        ! salfix(n) = 0.0
                end do
            end do
        end do

        if(n.gt.mpunr) &
            write(*,*) 'Dim max di indgeore e'' stata superata n=',n
        npuntr = n
            if (ismaster) write(6,329)npuntr
        329  format(5x,'For Salinity N ='i7)
        if(allocated(plth1)) deallocate(plth1)
        if(allocated(plth2)) deallocate(plth2)
    end if
    if(allocated(xpart)) deallocate(xpart)
    if(allocated(ypart)) deallocate(ypart)
    if(allocated(zpart)) deallocate(zpart)
    return
end subroutine topogr