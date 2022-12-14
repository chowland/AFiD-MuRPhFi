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
!     OPTION 1: CENTRE-FOCUSED CLUSTERING
!

    if (istr3==1) call centre_focus_grid(xc(1:nx),xm(1:nxm),nxm,alx3,str3)

!
!     OPTION 2: NATURAL TURB BL CLUSTERING
!

    if (istr3==2) call natural_BL_grid(xc(1:nx),xm(1:nxm),nxm,alx3)

!
!     OPTION 3: Symmetric NATURAL TURB BL CLUSTERING
!

    if (istr3==3) call sym_natural_BL_grid(xc(1:nx),xm(1:nxm),nxm,alx3, 1.0)

!
!     OPTION 4: HYPERBOLIC TANGENT-TYPE CLUSTERING
!

    if (istr3==4) call tanh_grid(xc(1:nx),xm(1:nxm),nxm,alx3,str3)

!
!   OPTION 5: SCALLOP-FOCUSED LOWER WALL CLUSTERING
!
    if (istr3==5) call scallop_grid(xc(1:nx), xm(1:nxm), nxm, alx3, dPdy, 0.5)

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

    call second_derivative_coeff(ap3ck, ac3ck, am3ck, xc(1:nx), alx3, 1, 1)

!
!    Q1/Q2 DIFFERENTIATION (STAGGERED VARIABLE)
!

    call second_derivative_coeff(ap3sk, ac3sk, am3sk, xm(1:nxm), alx3, inslwn, inslws)

!
!    TEMPERATURE DIFFERENTIATION
!

    call second_derivative_coeff(ap3ssk, ac3ssk, am3ssk, xm(1:nxm), alx3, TfixN, TfixS)

    return
end subroutine CreateGrid