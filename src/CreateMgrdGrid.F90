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
      implicit none

      real :: x1,x2,x3
      real :: a33, a33m, a33p
      real :: delet, etain, tstr3
      real :: z2dp

      real, allocatable, dimension(:) :: etaz, etazm

      integer :: i, j, kc, km, kp
      integer :: nxmo, nclip
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
       do  i=1,nzr
        x1=real(i-1)/real(nzmr)
        zcr(i)= zlen*x1
       end do

       do i=1,nzmr
         zmr(i)=(zcr(i)+zcr(i+1))*0.5d0
       end do
       zmr(0)=2.d0*zmr(1)-zmr(2)
       zmr(nzr)=2.d0*zmr(nzmr)-zmr(nzmr-1)

       do j=1,nyr
        x2=real(j-1)/real(nymr)
        ycr(j)= ylen*x2
       end do

       do j=1,nymr
        ymr(j)=(ycr(j)+ycr(j+1))*0.5d0
       end do
       ymr(0)=2.d0*ymr(1)-ymr(2)
       ymr(nyr)=2.d0*ymr(nymr)-ymr(nymr-1)

!
!     VERTICAL COORDINATE DEFINITION
!
!     OPTION 0: UNIFORM CLUSTERING
!
      call AllocateReal1DArray(etaz,1,nxr+500)
      call AllocateReal1DArray(etazm,1,nxr+500)

      if (istr3r.eq.0) then
        do kc=1,nxr
          x3=real(kc-1)/real(nxmr)
          etaz(kc)=alx3*x3
          xcr(kc)=etaz(kc)
        enddo
      endif

!
!     OPTION 4: HYPERBOLIC TANGENT-TYPE CLUSTERING
!

        tstr3=tanh(str3)

        if (istr3r.eq.4) then
         xcr(1)=0.0d0
         do kc=2,nxr
          z2dp=float(2*kc-nxr-1)/float(nxmr)
          xcr(kc)=(1+tanh(str3*z2dp)/tstr3)*0.5*alx3
          if(xcr(kc).lt.0.or.xcr(kc).gt.alx3)then
           write(*,*)'Refined grid is too streched: ','zc(',kc,')=',xcr(kc)
           stop
          endif
         end do
        end if

!
!     OPTION 6: CLIPPED CHEBYCHEV-TYPE CLUSTERING
!


      if(istr3r.eq.6) then
      nclip = int(str3)
      nxmo = nxr+nclip+nclip
      do kc=1,nxmo
        etazm(kc)=+cos(pi*(float(kc)-0.5)/float(nxmo))
      end do
      do kc=1,nxr
        etaz(kc)=etazm(kc+nclip)
      end do
      delet = etaz(1)-etaz(nxr)
      etain = etaz(1)
      do kc=1,nxr
        etaz(kc)=etaz(kc)/(0.5*delet)
      end do
      xcr(1) = 0.
      do kc=2,nxmr
        xcr(kc) = alx3*(1.-etaz(kc))*0.5
      end do
      xcr(nxr) = alx3
      endif

      call DestroyReal1DArray(etaz)
      call DestroyReal1DArray(etazm)
      
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
!

      do kc=2,nxmr-1
        km=kc-1
        a33=dxqr/g3rmr(kc)
        a33p=1.d0/d3xmr(kc)
        a33m=1.d0/d3xmr(km)
        ap3sskr(kc)=a33*a33p
        am3sskr(kc)=a33*a33m
        ac3sskr(kc)=-(ap3sskr(kc)+am3sskr(kc))
      enddo

      !CJH Lower wall BC
      kc = 1
      a33 = dxqr/g3rmr(kc)
      a33p=1.d0/d3xmr(kc)
      a33m=1.d0/g3rcr(kc) ! equivalent to virtual 1/d3xmr(0)
      ap3sskr(kc)=a33*a33p
      am3sskr(kc)=0.d0
      ac3sskr(kc)=-a33*(a33p + 2.d0*SfixS*a33m)

      !CJH Upper wall BC
      kc = nxmr
      km = kc - 1
      a33 = dxqr/g3rmr(kc)
      a33p=1.d0/d3xmr(kc)
      a33m=1.d0/d3xmr(km)
      ap3sskr(kc)=0.d0
      am3sskr(kc)=a33*a33m
      ac3sskr(kc)=-a33*(2.d0*SfixN*a33p + a33m)

      return                                                            
      end                                                               

