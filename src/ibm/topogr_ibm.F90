!
!     This routine finds the indices in the computational grid
!     close to the physical location of the body.
!
subroutine topogr
    use param
    use decomp_2d, only: xstart,xend,xstartr,xendr
    use ibm_param
    use mgrd_arrays, only: sal
    implicit none
    integer :: i,j,k,l,kstartp, nc, Npart
    integer :: km,kp,jm,jp,im,ip,mm
    
    real    :: xe, xem, xep
    real    :: ye, yem, yep
    real    :: ze, zem, zep
    real    :: delta1x, delta2x, r2, Lhex, radius, porosity
    integer,allocatable :: ind1(:), ind2(:)
    real,allocatable :: xpart(:), ypart(:)
    ! infig=1
    q1bo=0.d0
    q2bo=0.d0
    q3bo=0.d0
    densb=0.d0 
    allocate(plth1(1:nym,1:nzm))
    allocate(plth2(1:nym,1:nzm))
    plth1=1.d0
    plth2=1.d0
    
    ! if(ifrough.eq.1) then
    do i=1,nzm
        do j=1,nym
            !            plth1(j,i)=abs(dcos(10*pi*zm(i)/zlen) &
            !                      +dcos(10*pi*ym(j)/ylen))*0.01
            ! plth1(j,i)=dsin(2*pi/rlambda*ym(j)-pi/2)*rheight/2 + rheight/2
            !            plth2(j,i)=abs(dcos(10*pi*zm(i)/zlen)  &
            !                      +dcos(10*pi*ym(j)/ylen))*0.01
            ! plth2(j,i)=dsin(2*pi/rlambda*ym(j)-pi/2)*rheight/2 + rheight/2
            do k=1,4
                ye = ylen/8.0 + real(k - 1)*ylen/4.0
                plth1(j,i) = min(plth1(j,i), 0.1*(8.0/ylen)**2*(ym(j) - ye)**2)
            end do
            plth2(j,i) = plth1(j,i)
        end do
    end do
    ! endif
    
    ! if(ifrough.eq.2) then         ! ratchet
    !     do i=1,nzm
    !         do k=1,nr 
    !             do j=nym/nr*(k-1)+1,nym/nr*k
    !                 !              plth1(j,i)=0.025-(ym(j)-int(ym(j)/0.05)*0.05)/2.0
    !                 !              plth2(j,i)=(ym(j)-int(ym(j)/0.05)*0.05)/2.0
    !                 plth1(j,i)=0.025-(ym(j)-0.05d0*(k-1))/2.d0
    !                 plth2(j,i)=(ym(j)-0.05d0*(k-1))/2.d0
    !             enddo
    !         enddo
    !     enddo
    ! endif
    
    
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
    
    
    allocate(forclo(1:nx,xstart(2):xend(2),xstart(3):xend(3)))
    if (salinity) then
        allocate(forclor(1:nxmr,xstartr(2):xendr(2),xstartr(3):xendr(3)))
        allocate(solidr(1:nxmr,xstartr(2)-1:xendr(2)+1,xstartr(3)-1:xendr(3)+1))
    end if

    porosity = 0.426
    ! Interpret RayT input as target pore-scale Rayleigh number
    nc = nint((2.0*ylen*(1.0 - porosity)/3.0/pi)**0.5 * (RayS/RayT)**(1.0/3.0))
    if (ismaster) write(*,*) "Number of lattice columns: ",nc
    Lhex = 1.0/nc/ylen
    radius = (ylen*(1.0 - porosity)/6.0/pi)**0.5/nc
    Npart = (nc + 1)*(3*nc + 1) + 3*nc*nc
    allocate(xpart(1:Npart))
    allocate(ypart(1:Npart))

    i = 1
    ! do j=0,nc
    do j=0,nc
        xe = real(j)/real(nc)
        do k=0,3*nc
            ye = Lhex*real(k)
            xpart(i) = xe
            ypart(i) = ye
            i = i + 1
        end do
        if (j>0) then
            xe = (real(j) - 0.5)/real(nc)*alx3
            do k=1,3*nc
                ye = Lhex*(real(k) - 0.5)
                xpart(i) = xe
                ypart(i) = ye
                i = i + 1
            end do
        end if
    end do

    
    do l = 1,3 !{ start do over the 3 velocity components
        n=0
        
        forclo = 0.0
        !     l = 1   Q_1 vel. component (VZ)
        !     l = 2   Q_2 vel. component (VY)
        !     l = 3   Q_3 vel. component (VX)
        !
        ! Track location of solid circles:
        do i=xstart(3),xend(3)
            do j=xstart(2),xend(2)
                do k=1,nxm
                    xe = xm(k)
                    ye = ym(j)
                    if (l==3) xe = xc(k)
                    if (l==2) ye = yc(j)
                    do nc=1,Npart
                        ! x-position of solid centre
                        xem = xpart(nc)
                        yem = ypart(nc)
                        r2 = (xe - xem)**2 + (ye - yem)**2
                        if (r2<radius**2) then
                            forclo(k,j,i) = 1.0
                        end if
                    end do
                end do
            end do
        end do

        do i=xstart(3),xend(3)
            do j=xstart(2),xend(2)
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
                    ye = ym(j)
                    if (l==2) ye = yc(j)
                    !
                    !    SOLID PART
                    !
                    
                    ! if ((xe.lt.plth1(j,i))) then
                    if (forclo(k,j,i)>0.9) then
                        n=n+1
                        indgeo(l,n,1)=i
                        indgeo(l,n,2)=j
                        indgeo(l,n,3)=k
                        indgeoe(l,n,1)=i
                        indgeoe(l,n,2)=j
                        indgeoe(l,n,3)=k
                        distb(l,n)= 0.
                    ! elseif(xe.gt.(alx3-plth2(j,i))) then
                        
                    !     n=n+1
                    !     indgeo(l,n,1)=i
                    !     indgeo(l,n,2)=j
                    !     indgeo(l,n,3)=k
                    !     indgeoe(l,n,1)=i
                    !     indgeoe(l,n,2)=j
                    !     indgeoe(l,n,3)=k
                    !     distb(l,n)= 0.
                        
                    !    LOWER FLUID/PLATE BOUNDARY
                    !
                    ! elseif((xe.ge.plth1(j,i)).and.(xem.lt.plth1(j,i))) then
                    elseif (forclo(k,j,i)<forclo(km,j,i)) then
                        n=n+1
                        indgeo(l,n,1)=i
                        indgeo(l,n,2)=j
                        indgeo(l,n,3)=k
                        indgeoe(l,n,1)=i
                        indgeoe(l,n,2)=j 
                        indgeoe(l,n,3)=kp
                        delta1x=(xep-xe)
                        ! delta2x=(xe-plth1(j,i))
                        delta2x = 1.0
                        do nc=1,Npart
                            xem = xpart(nc) + sqrt(radius**2 - (ye - ypart(nc))**2)
                            if (xe > xem) delta2x = min(delta2x, xe - xem)
                        end do
                        distb(l,n) = delta2x/(delta1x+delta2x)
                                
                    !
                    !    UPPER FLUID/PLATE BOUNDARY
                    !
                    ! elseif((xe.le.(alx3-plth2(j,i))).and.(xep.gt.(alx3-plth2(j,i)))) &
                    !     then
                    elseif (forclo(k,j,i)<forclo(kp,j,i)) then
                        n=n+1
                        indgeo(l,n,1)=i
                        indgeo(l,n,2)=j
                        indgeo(l,n,3)=k
                        indgeoe(l,n,1)=i
                        indgeoe(l,n,2)=j 
                        indgeoe(l,n,3)=km
                        delta1x=(xe-xem)
                        ! delta2x=((alx3-plth2(j,i))-xe)
                        delta2x = 1.0
                        do nc=1,Npart
                            xep = xpart(nc) - sqrt(radius**2 - (ye - ypart(nc))**2)
                            if (xep > xe) delta2x = min(delta2x, xep - xe)
                        end do
                        distb(l,n) = delta2x/(delta1x+delta2x)
                                        
                    ! elseif (j==1.or.j==nym) then
                    !     if(ifnoslipy.eq.1) then
                    !         n=n+1
                    !         indgeo(l,n,1)=i
                    !         indgeo(l,n,2)=j
                    !         indgeo(l,n,3)=k
                    !         indgeoe(l,n,1)=i
                    !         indgeoe(l,n,2)=j
                    !         indgeoe(l,n,3)=k
                    !         distb(l,n)= 0.
                    !     endif
                    ! elseif (i==1.or.i==nzm) then
                    !     if(ifnoslipz.eq.1) then
                    !         n=n+1
                    !         indgeo(l,n,1)=i
                    !         indgeo(l,n,2)=j
                    !         indgeo(l,n,3)=k
                    !         indgeoe(l,n,1)=i
                    !         indgeoe(l,n,2)=j
                    !         indgeoe(l,n,3)=k
                    !         distb(l,n)= 0.
                    !     endif
                    end if
                                            
                end do
            end do
        end do
        iF(ANY(IsNaN(distb))) write(*,*) 'NaN in distb'
        
        if(l.eq.1) then
            if(n.gt.mpun)  &
            write(*,*) 'Dim max di indgeot e'' stata superata n=',n
            npunz= n
            !        write(6,332)npunz
            ! 332  format(5x,'For Q_1 N ='i7)
        end if
        if(l.eq.2) then
            if(n.gt.mpun)   &
            write(*,*) 'Dim max di indgeor e'' stata superata n=',n
            npuny= n
            !        write(6,331)npuny
            ! 331  format(5x,'For Q_2 N ='i7)
        end if
        if(l.eq.3) then
            if(n.gt.mpun)  &
            write(*,*) 'Dim max di indgeoz e'' stata superata n=',n
            npunx= n
            !        write(6,330)npunx
            ! 330  format(5x,'For Q_3 N ='i7)
        end if
    end do   !} end do over the 3 velocity components 
                            
                            

    !
    !     INDICES FOR TEMPERATURE
    !
    n=0
    forclo =0.0d0
    !

    do i=xstart(3),xend(3)
        do j=xstart(2),xend(2)
            do k=1,nxm
                km=kmv(k)
                kp=kpv(k)
                xe=xm(k)
                xem=xm(km)
                xep=xm(kp)
            !
            !    SOLID PART
            !           


                if(xe.lt.plth1(j,i)) then
                    n=n+1
                    indgeot(n,1)=i
                    indgeot(n,2)=j
                    indgeot(n,3)=k
                    indgeoet(n,1)=i
                    indgeoet(n,2)=j
                    indgeoet(n,3)=k
                    distbt(n)= 0.
                    temb(n) = tempbp(1,j,i)
                    !             forclo(k,j,i) = 1.

                elseif(xe.gt.(alx3-plth2(j,i))) then
                    n=n+1
                    indgeot(n,1)=i
                    indgeot(n,2)=j
                    indgeot(n,3)=k
                    indgeoet(n,1)=i
                    indgeoet(n,2)=j
                    indgeoet(n,3)=k
                    distbt(n)= 0.
                    temb(n) = temptp(1,j,i)
                    !             forclo(k,j,i) = 1.
                    !            end if
                
            !
            !    LOWER FLUID/PLATE BOUNDARY
            !
                elseif((xe.ge.plth1(j,i)).and.(xem.lt.plth1(j,i))) then
                    n=n+1
                    indgeot(n,1)=i
                    indgeot(n,2)=j
                    indgeot(n,3)=k
                    indgeoet(n,1)=i
                    indgeoet(n,2)=j 
                    indgeoet(n,3)=kp
                    delta1x=(xep-xe)
                    delta2x=(xe-plth1(j,i))
                    distbt(n)= delta2x/(delta1x+delta2x)
                    temb(n) = tempbp(1,j,i)
            !
            !    UPPER FLUID/PLATE BOUNDARY
            !
                elseif((xe.le.(alx3-plth2(j,i))).and.(xep.gt.(alx3-plth2(j,i)))) &
                    then
                    n=n+1
                    indgeot(n,1)=i
                    indgeot(n,2)=j
                    indgeot(n,3)=k
                    indgeoet(n,1)=i
                    indgeoet(n,2)=j 
                    indgeoet(n,3)=km
                    delta1x=(xe-xem)
                    delta2x=((alx3-plth2(j,i))-xe)
                    distbt(n)= delta2x/(delta1x+delta2x)
                    temb(n) = temptp(1,j,i)
                            
            ! elseif (j==1.or.j==nym) then
            !     if(ifnoslipy.eq.1) then      ! SL
            !         n=n+1
            !         indgeot(n,1)=i
            !         indgeot(n,2)=j
            !         indgeot(n,3)=k
            !         indgeoet(n,1)=i
            !         indgeoet(n,2)=j
            !         indgeoet(n,3)=k
            !         distbt(n)= 0.
            !     endif
            ! elseif (i==1.or.i==nzm) then
            !     if(ifnoslipz.eq.1) then        ! SL
            !         n=n+1
            !         indgeot(n,1)=i
            !         indgeot(n,2)=j
            !         indgeot(n,3)=k
            !         indgeoet(n,1)=i
            !         indgeoet(n,2)=j
            !         indgeoet(n,3)=k
            !         distbt(n)= 0.
            !     endif
                end if
                                
            end do
        end do
    end do
    if(n.gt.mpun) &
        write(*,*) 'Dim max di indgeote e'' stata superata n=',n
    npunte= n
    !        write(6,329)npunte
    ! 329  format(5x,'For Temperature N ='i7)
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
                    plth1(j,i) = min(plth1(j,i), 0.1*(8.0/ylen)**2*(ymr(j) - ye)**2)
                end do
                plth2(j,i) = plth1(j,i)
            end do
        end do
        !
        !     INDICES FOR SALINITY
        !
        n = 0
        forclor = 0.0d0
        !

        do i=xstartr(3)-1,xendr(3)+1
            do j=xstartr(2)-1,xendr(2)+1
                do k=1,nxmr
                    xe = xmr(k)
                    ye = ymr(j)
                    do nc=1,Npart
                        ! x-position of solid centre
                        xem = xpart(nc)
                        yem = ypart(nc)
                        r2 = (xe - xem)**2 + (ye - yem)**2
                        solidr(k,j,i) = solidr(k,j,i) .or. (r2<radius**2)
                        ! if (solidr(k,j,i)) then
                        !     forclor(k,j,i) = 1.0
                        ! end if
                    end do
                end do
            end do
        end do
        do i=xstartr(3),xendr(3)
            do j=xstartr(2),xendr(2)
                do k=1,nxmr
                    if (solidr(k,j,i)) then
                        forclor(k,j,i) = 1.0
                    end if
                end do
            end do
        end do

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

                    if (forclor(k,j,i) > 0.9) then
                        if (forclor(kp,j,i) < 0.1) then
                            n = n + 1
                            indgeor(n,1) = i
                            indgeor(n,2) = j
                            indgeor(n,3) = k
                            indgeoer(n,1) = i
                            indgeoer(n,2) = j
                            indgeoer(n,3) = kp
                            distbr(n) = 1.0
                        elseif (forclor(km,j,i) < 0.1) then
                            n = n + 1
                            indgeor(n,1) = i
                            indgeor(n,2) = j
                            indgeor(n,3) = k
                            indgeoer(n,1) = i
                            indgeoer(n,2) = j
                            indgeoer(n,3) = km
                            distbr(n) = 1.0
                        else
                            ! indgeoer(n,3) = k
                            ! distbr(n) = 0.0
                            ! sal(k,j,i) = 0.0
                        end if
                        ! salfix(n) = 0.0

                !   FLUID
                !   -----   BOUNDARY
                !   SOLID

                    ! elseif (forclor(k,j,i) < forclor(km,j,i)) then
                    !     n = n+1
                    !     indgeor(n,1) = i
                    !     indgeor(n,2) = j
                    !     indgeor(n,3) = k
                    !     indgeoer(n,1) = i
                    !     indgeoer(n,2) = j 
                    !     indgeoer(n,3) = kp
                    !     delta1x = (xep - xe)
                    !     delta2x = 1.0
                    !     do nc=1,Npart
                    !         xem = xpart(nc) + sqrt(radius**2 - (ye - ypart(nc))**2)
                    !         if (xe > xem) delta2x = min(delta2x, xe - xem)
                    !     end do
                    !     ! distbr(n) = delta2x/(delta1x + delta2x)
                    !     distbr(n) = 1.0
                    !     salfix(n) = 0.0

                !   SOLID
                !   -----   BOUNDARY
                !   FLUID

                    ! elseif (forclor(k,j,i) < forclor(kp,j,i)) then
                    !     n = n+1
                    !     indgeor(n,1) = i
                    !     indgeor(n,2) = j
                    !     indgeor(n,3) = k
                    !     indgeoer(n,1) = i
                    !     indgeoer(n,2) = j 
                    !     indgeoer(n,3) = km
                    !     delta1x = (xe - xem)
                    !     do nc=1,Npart
                    !         xep = xpart(nc) - sqrt(radius**2 - (ye - ypart(nc))**2)
                    !         if (xep > xe) delta2x = min(delta2x, xep - xe)
                    !     end do
                    !     ! distbr(n)= delta2x/(delta1x + delta2x)
                    !     distbr(n) = 1.0
                    !     salfix(n) = 0.0
                    
                    end if
                                    
                end do
            end do
        end do

        if(n.gt.mpun) &
            write(*,*) 'Dim max di indgeore e'' stata superata n=',n
        npuntr = n
        !        write(6,329)npuntr
        ! 329  format(5x,'For Salinity N ='i7)
        if(allocated(plth1)) deallocate(plth1)
        if(allocated(plth2)) deallocate(plth2)
    end if
    if(allocated(xpart)) deallocate(xpart)
    if(allocated(ypart)) deallocate(ypart)
    return
end subroutine topogr