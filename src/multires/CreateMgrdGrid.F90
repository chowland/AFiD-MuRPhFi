!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: CreateMgrdGrid.F90                             !
!    CONTAINS: subroutine CreateMgrdGrid                  !
!                                                         ! 
!    PURPOSE: Compute the mgrd indices, grid,             !
!     grid metrics and coefficients for differentiation   !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine CreateMgrdGrid
    use param
    use AuxiliaryRoutines
    use GridModule
    implicit none

    real :: x1,x2,x3
    real :: a33, a33m, a33p

    integer :: i, j, kc, km, kp
    logical :: fexist

    do kc=1,nxmr
        kmvr(kc)=kc-1
        kpvr(kc)=kc+1
        if(kc.eq.1) kmvr(kc)=kc
        if(kc.eq.nxmr) kpvr(kc)=kc
    end do

    do kc=1,nxmr
        kpcr(kc)=kpvr(kc)-kc
        kmcr(kc)=kc-kmvr(kc)
    end do


!
!     UNIFORM (HORIZONTAL DIRECTIONS) GRID
!

    call uniform_grid(zcr(1:nzr), zmr(1:nzmr), nzmr, zlen)

    zmr(0)=2.d0*zmr(1)-zmr(2)
    zmr(nzr)=2.d0*zmr(nzmr)-zmr(nzmr-1)

    call uniform_grid(ycr(1:nyr), ymr(1:nymr), nymr, ylen)

    ymr(0)=2.d0*ymr(1)-ymr(2)
    ymr(nyr)=2.d0*ymr(nymr)-ymr(nymr-1)

!
!     VERTICAL COORDINATE DEFINITION
!
!     OPTION 0: UNIFORM CLUSTERING
!

    if (istr3r==0) call uniform_grid(xcr(1:nxr),xmr(1:nxmr),nxmr,alx3)


!
!     OPTION 4: HYPERBOLIC TANGENT-TYPE CLUSTERING
!

    if (istr3r==4) call tanh_grid(xcr(1:nxr),xmr(1:nxmr),nxmr,alx3,str3)

!
!     OPTION 6: CLIPPED CHEBYCHEV-TYPE CLUSTERING
!

    if (istr3r==6) call cheb_grid(xcr(1:nxr),xmr(1:nxmr),nxmr,alx3,str3)

!
!     OPTION 7: As option 6, but only for high resolution at one (lower) wall
!

    if (istr3r==7) call asym_cheb_grid(xcr(1:nxr),xmr(1:nxmr),nxmr,alx3,str3)

!m-----------------------------------------
!
!     METRIC FOR UNIFORM DIRECTIONS
!

      dxr=real(nxmr)/alx3
      dyr=real(nymr)/ylen
      dzr=real(nzmr)/zlen

      dxqr=dxr*dxr
      dyqr=dyr*dyr
      dzqr=dzr*dzr

!
!     STAGGERED COORDINATES AND
!     METRIC QUANTITIES FOR NON-UNIFORM 
!     DIRECTIONS
!

      do kc=1,nxmr
        xmr(kc)=(xcr(kc)+xcr(kc+1))*0.5d0
      enddo
      xmr(0)=2.d0*xcr(1)-xmr(1)
      xmr(nxr)=2.d0*xcr(nxr)-xmr(nxmr)

      do kc=2,nxmr
        g3rcr(kc)=(xcr(kc+1)-xcr(kc-1))*dxr*0.5d0
        g3rmr(kc)=(xmr(kc+1)-xmr(kc-1))*dxr*0.5d0
      enddo
      g3rmr(1) = (xmr(2)-(2*xcr(1)-xmr(1)))*dxr*0.5d0
      g3rcr(1)=(xcr(2)-xcr(1))*dxr
      g3rcr(nxr)= (xcr(nxr)-xcr(nxmr))*dxr

      do kc=1,nxmr
        d3xcr(kc) = (xcr(kc+1) - xcr(kc))*dxr
        d3xmr(kc) = (xmr(kc+1) - xmr(kc))*dxr
      end do
!
!     WRITE GRID INFORMATION
!
      do kc=1,nxmr
        udx3mr(kc) = dxr/d3xcr(kc)
        udx3cr(kc) = dxr/g3rcr(kc)
      end do
      udx3cr(nxr) = dxr/g3rcr(nxr)
!m====================================================
      if(ismaster) then
      open(unit=78,file='outputdir/axicorr.out',status='unknown')
      do kc=1,nxr
        write(78,345) kc,xcr(kc),xmr(kc),g3rcr(kc),g3rmr(kc)
      end do
      close(78)
 345  format(i4,4(2x,e23.15))
!m===================================================
!
!     QUANTITIES FOR DERIVATIVES
!
      open(unit=78,file='outputdir/fact3r.out',status='unknown')
      do kc=1,nxmr
        write(78,*) kc,udx3mr(kc),udx3cr(kc)
      end do
        write(78,*) nxr,udx3mr(nxmr),udx3cr(nxr)
      close(78)
      endif

!
!    COEFFICIENTS FOR SALINITY DIFFERENTIATION FOR NON-UNIFORM GRID
!
!    Q3 DIFFERENTIATION (CENTERED VARIABLE)
!

      ! am3ckr(1)=0.d0
      ! ap3ckr(1)=0.d0
      ! ac3ckr(1)=1.d0
      ! am3ckr(nxr)=0.d0
      ! ap3ckr(nxr)=0.d0
      ! ac3ckr(nxr)=1.d0

      ! do kc=2,nxmr
      !  km=kc-1
      !  kp=kc+1
      !  a33=dxqr/g3rcr(kc)
      !  a33p=1.d0/g3rmr(kc)
      !  a33m=1.d0/g3rmr(km)
      !  ap3ckr(kc)=a33*a33p
      !  am3ckr(kc)=a33*a33m
      !  ac3ckr(kc)=-(ap3ckr(kc)+am3ckr(kc))
      ! enddo

!CS !
!CS !    Q1/Q2 DIFFERENTIATION (STAGGERED VARIABLE)
!CS !
!CS !
!CS 
!CS       do kc=2,nxm-1
!CS       kp=kc+1
!CS       km=kc-1
!CS       a33=dxq/g3rm(kc)
!CS       a33p= +a33/g3rc(kp)
!CS       a33m= +a33/g3rc(kc)
!CS       ap3sk(kc)=a33p
!CS       am3sk(kc)=a33m
!CS       ac3sk(kc)=-(ap3sk(kc)+am3sk(kc))
!CS       enddo
!CS !    
!CS !    LOWER WALL BOUNDARY CONDITIONS (INSLWS SETS NO-SLIP vs STRESS-FREE WALL)
!CS !    
!CS       kc=1
!CS       kp=kc+1
!CS       a33=dxq/g3rm(kc)
!CS       a33p= +a33/g3rc(kp)
!CS       a33m= +a33/g3rc(kc)
!CS       ap3sk(kc)=a33p
!CS       am3sk(kc)=0.d0
!CS       ac3sk(kc)=-(a33p+inslws*a33m*2.d0)
!CS 
!CS !    
!CS !    UPPER WALL BOUNDARY CONDITIONS (INSLWN SETS NO-SLIP vs STRESS-FREE WALL)
!CS !    
!CS 
!CS       kc=nxm
!CS       kp=kc+1
!CS       a33=dxq/g3rm(kc)
!CS       a33p= +a33/g3rc(kp)
!CS       a33m= +a33/g3rc(kc)
!CS       am3sk(kc)=a33m
!CS       ap3sk(kc)=0.d0
!CS       ac3sk(kc)=-(a33m+inslwn*a33p*2.d0)

      ! am3sskr(1)=0.d0
      ! ap3sskr(1)=0.d0
      ! ac3sskr(1)=1.d0

!
!    SALINITY DIFFERENTIATION
      !CJH (now staggered!)
      !CJH now with additional phi stencil to allow for various BCs
!

      do kc=2,nxmr-1
        km=kc-1
        a33=dxqr/g3rmr(kc)
        a33p=1.d0/d3xmr(kc)
        a33m=1.d0/d3xmr(km)
        ap3sskr(kc)=a33*a33p
        am3sskr(kc)=a33*a33m
        ac3sskr(kc)=-(ap3sskr(kc)+am3sskr(kc))
        ap3spkr(kc)=a33*a33p
        am3spkr(kc)=a33*a33m
        ac3spkr(kc)=-(ap3spkr(kc)+am3spkr(kc))
      enddo

      !CJH Lower wall BC
      kc = 1
      a33 = dxqr/g3rmr(kc)
      a33p=1.d0/d3xmr(kc)
      a33m=1.d0/g3rcr(kc) ! equivalent to virtual 1/d3xmr(0)
      ap3sskr(kc)=a33*a33p
      am3sskr(kc)=0.d0
      ac3sskr(kc)=-a33*(a33p + 2.d0*SfixS*a33m)
      ap3spkr(kc)=a33*a33p
      am3spkr(kc)=0.d0
      ac3spkr(kc)=-a33*a33p ! Sets d/dx(phi)=0 on boundary

      !CJH Upper wall BC
      kc = nxmr
      km = kc - 1
      a33 = dxqr/g3rmr(kc)
      a33p=1.d0/d3xmr(kc)
      a33m=1.d0/d3xmr(km)
      ap3sskr(kc)=0.d0
      am3sskr(kc)=a33*a33m
      ac3sskr(kc)=-a33*(2.d0*SfixN*a33p + a33m)
      ap3spkr(kc)=0.d0
      am3spkr(kc)=a33*a33m
      ac3spkr(kc)=-a33*a33m ! Sets d/dx(phi)=0 on boundary

      return                                                            
      end                                                               

