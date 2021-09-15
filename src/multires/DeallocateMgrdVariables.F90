!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: DeallocateMgrdVariables.F90                    !
!    CONTAINS: subroutine DeallocateMgrdVariables         !
!                                                         ! 
!    PURPOSE: Finalization routine. Deallocates all       !
!     variables used in the multi-res extension           !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine DeallocateMgrdVariables
    use param
    use mgrd_arrays
    use decomp_2d
    use AuxiliaryRoutines
    implicit none

    call DestroyReal1DArray(zcr)
    call DestroyReal1DArray(zmr)
    call DestroyReal1DArray(ycr)
    call DestroyReal1DArray(ymr)
    call DestroyReal1DArray(xcr)
    call DestroyReal1DArray(xmr)
    call DestroyReal1DArray(g3rcr)
    call DestroyReal1DArray(g3rmr)
    call DestroyReal1DArray(d3xcr)
    call DestroyReal1DArray(d3xmr)
    call DestroyReal1DArray(udx3cr)
    call DestroyReal1DArray(udx3mr)

    call DestroyReal1DArray(ap3ckr)
    call DestroyReal1DArray(ac3ckr)
    call DestroyReal1DArray(am3ckr)

    call DestroyReal1DArray(ap3sskr)
    call DestroyReal1DArray(ac3sskr)
    call DestroyReal1DArray(am3sskr)

    call DestroyReal1DArray(ap3spkr)
    call DestroyReal1DArray(ac3spkr)
    call DestroyReal1DArray(am3spkr)

    call DestroyInt1dArray(kmcr)
    call DestroyInt1dArray(kpcr)
    call DestroyInt1dArray(kmvr)
    call DestroyInt1dArray(kpvr)

    call DestroyInt1dArray(irangs)   !CS mgrd
    call DestroyInt1dArray(jrangs)   !CS mgrd
    call DestroyInt1dArray(krangs)   !CS mgrd

    call DestroyInt1dArray(irangc)   !CS mgrd
    call DestroyInt1dArray(jrangc)   !CS mgrd
    call DestroyInt1dArray(krangc)   !CS mgrd

    call DestroyInt1dArray(irangr)   !CS mgrd
    call DestroyInt1dArray(jrangr)   !CS mgrd
    call DestroyInt1dArray(krangr)   !CS mgrd

    call DestroyReal2DArray(cxvx) !CS mgrd
    call DestroyReal2DArray(cxvy) !CS mgrd
    call DestroyReal2DArray(cxvz) !CS mgrd

    call DestroyReal2DArray(cyvx) !CS mgrd
    call DestroyReal2DArray(cyvy) !CS mgrd
    call DestroyReal2DArray(cyvz) !CS mgrd

    call DestroyReal2DArray(czvx) !CS mgrd
    call DestroyReal2DArray(czvy) !CS mgrd
    call DestroyReal2DArray(czvz) !CS mgrd

    call DestroyReal2DArray(cxrs) !CS mgrd
    call DestroyReal2DArray(cyrs) !CS mgrd
    call DestroyReal2DArray(czrs) !CS mgrd

    call DestroyReal2DArray(cxsalc)
    call DestroyReal2DArray(cysalc)
    call DestroyReal2DArray(czsalc)

    call DestroyReal3DArray(vxr)
    call DestroyReal3DArray(vyr)
    call DestroyReal3DArray(vzr)

    call DestroyReal3DArray(tpdv)
    call DestroyReal3DArray(tpdvr)  !CS BUG: ERROR WHILE DEALLOCATING

    return
end subroutine DeallocateMgrdVariables