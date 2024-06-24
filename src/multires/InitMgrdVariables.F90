!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: InitVariables.F90                              !
!    CONTAINS: subroutine InitVariables                   !
!                                                         ! 
!    PURPOSE: Initialization routine. Sets to zero all    !
!     variables used in the multiresolution code          !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine InitMgrdVariables
    use param
    use mgrd_arrays
    use decomp_2d
    use AuxiliaryRoutines
    implicit none
    
!-------------------------------------------------
! Arrays for grid making
!-------------------------------------------------

    call AllocateReal1DArray(zcr,1,nzr)
    call AllocateReal1DArray(zmr,0,nzr+1)

    call AllocateReal1DArray(ycr,1,nyr)
    call AllocateReal1DArray(ymr,0,nyr+1)

    call AllocateReal1DArray(xcr,1,nxr)
    call AllocateReal1DArray(xmr,0,nxr)

    call AllocateReal1DArray(g3rcr,1,nxr)
    call AllocateReal1DArray(g3rmr,1,nxr)

    call AllocateReal1DArray(d3xcr,1,nxr)
    call AllocateReal1DArray(d3xmr,1,nxr)

    call AllocateReal1DArray(udx3cr,1,nxr)
    call AllocateReal1DArray(udx3mr,1,nxr)

    ! call AllocateReal1DArray(ap3ckr,1,nxr)
    ! call AllocateReal1DArray(ac3ckr,1,nxr)
    ! call AllocateReal1DArray(am3ckr,1,nxr)

    ! call AllocateReal1DArray(ap3sskr,1,nxr)
    ! call AllocateReal1DArray(ac3sskr,1,nxr)
    ! call AllocateReal1DArray(am3sskr,1,nxr)

    ! call AllocateReal1DArray(ap3spkr,1,nxr)
    ! call AllocateReal1DArray(ac3spkr,1,nxr)
    ! call AllocateReal1DArray(am3spkr,1,nxr)

    call AllocateInt1dArray(kmcr,1,nxr)
    call AllocateInt1dArray(kpcr,1,nxr)
    call AllocateInt1dArray(kmvr,1,nxr)
    call AllocateInt1dArray(kpvr,1,nxr)

    call AllocateInt1dArray(irangs,0,nx)
    call AllocateInt1dArray(jrangs,0,ny)
    call AllocateInt1dArray(krangs,0,nz)

    call AllocateInt1dArray(irangc,0,nx+1)
    call AllocateInt1dArray(jrangc,0,ny)
    call AllocateInt1dArray(krangc,0,nz)

    call AllocateInt1dArray(irangb,0,nx+1)
    call AllocateInt1dArray(jrangb,0,ny)
    call AllocateInt1dArray(krangb,0,nz)

    call AllocateInt1dArray(irangr,0,nxr)
    call AllocateInt1dArray(jrangr,0,nyr)
    call AllocateInt1dArray(krangr,0,nzr)

    call AllocateInt1dArray(yc_to_ymr,0,ny)
    call AllocateInt1dArray(zc_to_zmr,0,nz)

    call AllocateReal2DArray(cxvx,1,4,1,nxr)
    call AllocateReal2DArray(cxvy,1,4,1,nxmr)
    call AllocateReal2DArray(cxvz,1,4,1,nxmr)

    call AllocateReal2DArray(cyvx,1,4,1,nymr)
    call AllocateReal2DArray(cyvy,1,4,1,nymr)
    call AllocateReal2DArray(cyvz,1,4,1,nymr)

    call AllocateReal2DArray(czvx,1,4,1,nzmr)
    call AllocateReal2DArray(czvy,1,4,1,nzmr)
    call AllocateReal2DArray(czvz,1,4,1,nzmr)

    call AllocateReal2DArray(cxrs,1,4,1,nxmr)
    call AllocateReal2DArray(cyrs,1,4,1,nymr)
    call AllocateReal2DArray(czrs,1,4,1,nzmr)

    call AllocateReal2DArray(cxphic,1,4,1,nxm)
    call AllocateReal2DArray(cyphic,1,4,1,nym)
    call AllocateReal2DArray(czphic,1,4,1,nzm)

    call AllocateReal2DArray(cxsalc,1,4,1,nx)
    call AllocateReal2DArray(cysalc,1,4,1,nym)
    call AllocateReal2DArray(czsalc,1,4,1,nzm)


    if (IBM) then
        if (phasefield) then
            call AllocateReal2DArray(cych,1,4,1,nym)
            call AllocateReal2DArray(czch,1,4,1,nzm)
        end if
    end if

!-------------------------------------------------
! Arrays with ghost cells
!-------------------------------------------------

    call AllocateReal3DArray(tpdv,-1,nx+1,xstart(2)-2,xend(2)+2,xstart(3)-2,xend(3)+2)
    !CS   For tpdvr, larger array needed to prevent memory overflow in InterpVelMgrd
    call AllocateReal3DArray(tpdvr,-1,nxr+1,xstartr(2)-2,xendr(2)+2,xstartr(3)-2,xendr(3)+2)
    call AllocateReal3DArray(temp_pdvr,-1,nxr+1,xstartr(2)-2,xendr(2)+2,xstartr(3)-2,xendr(3)+2)

    ! RHS array without ghost cells
    call AllocateReal3DArray(rhsr,1,nxmr,xstartr(2),xendr(2),xstartr(3),xendr(3))

    return
end subroutine InitMgrdVariables