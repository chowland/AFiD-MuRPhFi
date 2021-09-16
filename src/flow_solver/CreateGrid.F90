!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: CreateGrid.F90                                 !
!    CONTAINS: subroutine CreateGrid                      !
!                                                         ! 
!    PURPOSE: Compute the indices, grid, grid metrics     !
!     and coefficients for differentiation                !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine CreateGrid
    use param
    use AuxiliaryRoutines
    use GridModule
    implicit none

    real :: x1,x2,x3
    real :: a33, a33m, a33p

    integer :: i, j, kc, km, kp
    logical :: fexist

    do kc=1,nxm
        kmv(kc)=kc-1
        kpv(kc)=kc+1
        if(kc.eq.1) kmv(kc)=kc
        if(kc.eq.nxm) kpv(kc)=kc
    end do

    do kc=1,nxm
        kpc(kc)=kpv(kc)-kc
        kmc(kc)=kc-kmv(kc)
    end do


!
!     UNIFORM (HORIZONTAL DIRECTIONS) GRID
!

    call uniform_grid(zc(1:nz), zm(1:nzm), nzm, zlen)

    zm(0) = 2.d0*zm(1) - zm(2)
    zm(nz) = 2.d0*zm(nzm) - zm(nzm-1)

    call uniform_grid(yc(1:ny), ym(1:nym), nym, ylen)

    ym(0) = 2.d0*ym(1) - ym(2)
    ym(ny) = 2.d0*ym(nym) - ym(nym-1)

!
!     VERTICAL COORDINATE DEFINITION
!
!     OPTION 0: UNIFORM CLUSTERING
!

    if (istr3==0) call uniform_grid(xc(1:nx),xm(1:nxm),nxm,alx3)

!
!     OPTION 4: HYPERBOLIC TANGENT-TYPE CLUSTERING
!

    if (istr3==4) call tanh_grid(xc(1:nx),xm(1:nxm),nxm,alx3,str3)

!
!     OPTION 6: CLIPPED CHEBYCHEV-TYPE CLUSTERING
!

    if (istr3==6) call cheb_grid(xc(1:nx),xm(1:nxm),nxm,alx3,str3)

!
!     OPTION 7: As option 6, but only for high resolution at one (lower) wall
!

    if (istr3==7) call asym_cheb_grid(xc(1:nx),xm(1:nxm),nxm,alx3,str3)

!m-----------------------------------------
!
!     METRIC FOR UNIFORM DIRECTIONS
!

    dx=real(nxm)/alx3
    dy=real(nym)/ylen
    dz=real(nzm)/zlen

    dxq=dx*dx                                                      
    dyq=dy*dy                                                      
    dzq=dz*dz                                                      

!
!     STAGGERED COORDINATES AND
!     METRIC QUANTITIES FOR NON-UNIFORM 
!     DIRECTIONS
!

    do kc=1,nxm
        xm(kc)=(xc(kc)+xc(kc+1))*0.5d0
    end do
    xm(nx) = 2*xc(nx) - xm(nxm)
    do kc=2,nxm
        g3rc(kc)=(xc(kc+1)-xc(kc-1))*dx*0.5d0
        g3rm(kc)=(xm(kc+1)-xm(kc-1))*dx*0.5d0
    end do
    !CJH virtual xm(0) = 2*xc(1) - xm(1)
    g3rm(1) = (xm(2)-(2*xc(1)-xm(1)))*dx*0.5d0
    g3rc(1)=(xc(2)-xc(1))*dx
    g3rc(nx)= (xc(nx)-xc(nxm))*dx
    
    do kc=1,nxm
        d3xc(kc) = (xc(kc+1) - xc(kc))*dx
        d3xm(kc) = (xm(kc+1) - xm(kc))*dx
    end do
!
!     WRITE GRID INFORMATION
!
    do kc=1,nxm
        udx3m(kc) = dx/d3xc(kc)
        udx3c(kc) = dx/g3rc(kc)
    end do
    udx3c(nx) = dx/g3rc(nx)
!m====================================================
    if(ismaster) then
        open(unit=78,file='outputdir/axicor.out',status='unknown')
        do kc=1,nx
            write(78,345) kc,xc(kc),xm(kc),g3rc(kc),g3rm(kc)
        end do
        close(78)
   345  format(i4,4(2x,e23.15))
!m===================================================
!
!     QUANTITIES FOR DERIVATIVES
!
        open(unit=78,file='outputdir/fact3.out',status='unknown')
        do kc=1,nxm
            write(78,*) kc,udx3m(kc),udx3c(kc)
        end do
        write(78,*) nx,udx3m(nxm),udx3c(nx)
        close(78)
    end if

!
!    COEFFICIENTS FOR DIFFERENTIATION FOR NON-UNIFORM GRID
!
!    Q3 DIFFERENTIATION (CENTERED VARIABLE)
!

    am3ck(1)=0.d0
    ap3ck(1)=0.d0
    ac3ck(1)=1.d0
    am3ck(nx)=0.d0
    ap3ck(nx)=0.d0
    ac3ck(nx)=1.d0

    do kc=2,nxm
       km=kc-1
       a33=dxq/g3rc(kc)
       a33p=1.d0/d3xc(kc)
       a33m=1.d0/d3xc(km)
       ap3ck(kc)=a33*a33p
       am3ck(kc)=a33*a33m
       ac3ck(kc)=-(ap3ck(kc)+am3ck(kc))
    end do

!
!    Q1/Q2 DIFFERENTIATION (STAGGERED VARIABLE)
!
!

    do kc=2,nxm-1
        km=kc-1
        a33=dxq/g3rm(kc)
        a33p=1.d0/d3xm(kc)
        a33m=1.d0/d3xm(km)
        ap3sk(kc)=a33*a33p
        am3sk(kc)=a33*a33m
        ac3sk(kc)=-(ap3sk(kc)+am3sk(kc))
    end do
!    
!    LOWER WALL BOUNDARY CONDITIONS (INSLWS SETS NO-SLIP vs STRESS-FREE WALL)
!    
    kc=1
    a33=dxq/g3rm(kc)
    a33p=1.d0/d3xm(kc)
    a33m=1.d0/g3rc(kc) ! equivalent to virtual 1/d3xm(0)
    ap3sk(kc)=a33*a33p
    am3sk(kc)=0.d0
    ac3sk(kc)=-a33*(a33p+2.d0*inslws*a33m)

!    
!    UPPER WALL BOUNDARY CONDITIONS (INSLWN SETS NO-SLIP vs STRESS-FREE WALL)
!    

    kc=nxm
    km=kc-1
    a33=dxq/g3rm(kc)
    a33p=1.d0/d3xm(kc)
    a33m=1.d0/d3xm(km)
    ap3sk(kc)=0.d0
    am3sk(kc)=a33*a33m
    ac3sk(kc)=-a33*(2.d0*inslwn*a33p+a33m)

!
!    TEMPERATURE DIFFERENTIATION
!CJH (now staggered!)
!

    do kc=2,nxm-1
        ap3ssk(kc)=ap3sk(kc)
        am3ssk(kc)=am3sk(kc)
        ac3ssk(kc)=ac3sk(kc)
    end do

    !CJH Lower wall BC
    kc = 1
    a33 = dxq/g3rm(kc)
    a33p = 1.d0/d3xm(kc)
    a33m = 1.d0/g3rc(kc) ! equivalent to virtual 1/d3xm(0)
    ap3ssk(kc) = a33*a33p
    am3ssk(kc) = 0.d0
    ac3ssk(kc) = -a33*(a33p + 2.d0*TfixS*a33m)

    !CJH Upper wall BC
    kc = nxm
    km = kc - 1
    a33 = dxq/g3rm(kc)
    a33p = 1.d0/d3xm(kc)
    a33m = 1.d0/d3xm(km)
    ap3ssk(kc) = 0.d0
    am3ssk(kc) = a33*a33m
    ac3ssk(kc) = -a33*(a33m + 2.d0*TfixN*a33p)

      return                                                            
end subroutine CreateGrid