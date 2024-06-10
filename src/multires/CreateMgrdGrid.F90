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
    use afid_salinity, only: SfixS, SfixN, PraS, ap3sskr, ac3sskr, am3sskr
    use afid_phasefield, only: ap3spkr, ac3spkr, am3spkr
    implicit none
    real, dimension(nxmr) :: ap3sskr_D,ac3sskr_D, am3sskr_D, ap3sskr_N,ac3sskr_N, am3sskr_N
    
    integer :: kc,i

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
!     OPTION 1: CENTRE-FOCUSED CLUSTERING
!

    if (istr3r==1) call centre_focus_grid(xcr(1:nxr),xmr(1:nxmr),nxmr,alx3,str3)

!
!     OPTION 3: Symmetric NATURAL TURB BL CLUSTERING
!

    if (istr3r==3) call sym_natural_BL_grid(xcr(1:nxr),xmr(1:nxmr),nxmr,alx3, PraS)

!
!     OPTION 4: HYPERBOLIC TANGENT-TYPE CLUSTERING
!

    if (istr3r==4) call tanh_grid(xcr(1:nxr),xmr(1:nxmr),nxmr,alx3,str3)

!
!   OPTION 5: SCALLOP-FOCUSED LOWER WALL CLUSTERING
!
    if (istr3==5) call scallop_grid(xcr(1:nxr), xmr(1:nxmr), nxmr, alx3, dPdy, 0.5)

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
    end if

!
!    COEFFICIENTS FOR SALINITY DIFFERENTIATION FOR NON-UNIFORM GRID
!

    ! Salinity differentiation
    if (salinity) then
    if ( FixValueBCRegion_Length==0) then

        call second_derivative_coeff(ap3sskr_D, ac3sskr_D, am3sskr_D, xmr(1:nxmr), alx3, SfixN, SfixS)
        call second_derivative_coeff(ap3sskr_N, ac3sskr_N,  am3sskr_N, xmr(1:nxmr), alx3, SfixN, SfixS)

    else if (FixValueBCRegion_Length/=0) then

        if ( FixValueBCRegion_Nord_or_Sud==0) then
        
            call second_derivative_coeff(ap3sskr_D, ac3sskr_D, am3sskr_D, xmr(1:nxmr), alx3, SfixN, 1)
            call second_derivative_coeff(ap3sskr_N, ac3sskr_N, am3sskr_N, xmr(1:nxmr), alx3, SfixN, 0)
        
        else if( FixValueBCRegion_Nord_or_Sud==1) then
  
            call second_derivative_coeff(ap3sskr_D, ac3sskr_D, am3sskr_D, xmr(1:nxmr), alx3, 1, SfixS)
            call second_derivative_coeff(ap3sskr_N, ac3sskr_N, am3sskr_N, xmr(1:nxmr), alx3, 0, SfixS)
           
        end if
   
    end if
    
    do kc = 1, nxmr
        ap3sskr(kc,1) = ap3sskr_D(kc)
        ap3sskr(kc,2) = ap3sskr_N(kc)
        
        ac3sskr(kc,1) = ac3sskr_D(kc)
        ac3sskr(kc,2) = ac3sskr_N(kc)
        
        am3sskr(kc,1) = am3sskr_D(kc)
        am3sskr(kc,2) = am3sskr_N(kc)
        
   end do

   end if 
    
   nyr_Cold = 0
   nyr_Hot = 0
   do i = 1, nymr
       if( FixValueBCRegion_Length/=0 .and. ymr(i) <= 0.01 * FixValueBCRegion_Length * YLEN) then
            nyr_Cold = nyr_Cold+1
       else if  (FixValueBCRegion_Length/=0 .and. ymr(i) >= YLEN - 0.01 * FixValueBCRegion_Length * YLEN) then
            nyr_Hot = nyr_Hot +1       
        end if
    end do 

   
    ! Phase-field differentiation (ensuring zero gradient at boundaries)
    if (phasefield) call second_derivative_coeff(ap3spkr, ac3spkr, am3spkr, xmr(1:nxmr), alx3, 0, 0)

    return
end subroutine CreateMgrdGrid