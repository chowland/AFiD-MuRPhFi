!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: InitSalVariables.F90                           !
!    CONTAINS: subroutine InitSalVariables                !
!                                                         ! 
!    PURPOSE: Initialization routine. Sets to zero all    !
!     variables loaded when using salinity                !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine InitSalVariables
    use param
    use mgrd_arrays
    use decomp_2d
    use AuxiliaryRoutines
    implicit none
    
    ! Boundary conditions
    call AllocateReal3DArray(salbp,1,1,xstartr(2)-lvlhalo,xendr(2)+lvlhalo,xstartr(3)-lvlhalo,xendr(3)+lvlhalo)
    call AllocateReal3DArray(saltp,1,1,xstartr(2)-lvlhalo,xendr(2)+lvlhalo,xstartr(3)-lvlhalo,xendr(3)+lvlhalo)

    ! Main array with ghost cells
    call AllocateReal3DArray(sal,1,nxr,xstartr(2)-lvlhalo,xendr(2)+lvlhalo,xstartr(3)-lvlhalo,xendr(3)+lvlhalo)

    ! Arrays without ghost cells
    call AllocateReal3DArray(rhsr,1,nxr,xstartr(2),xendr(2),xstartr(3),xendr(3))
    call AllocateReal3DArray(rusal,1,nxr,xstartr(2),xendr(2),xstartr(3),xendr(3))
    call AllocateReal3DArray(hsal,1,nxr,xstartr(2),xendr(2),xstartr(3),xendr(3))

    ! Coarse array
    call AllocateReal3DArray(salc,1,nx,xstart(2)-lvlhalo,xend(2)+lvlhalo,xstart(3)-lvlhalo,xend(3)+lvlhalo)

    !CJH Needed for melt boundary condition
    if (melt) then
        call AllocateReal3DArray(tempr,1,1,xstartr(2)-lvlhalo,xendr(2)+lvlhalo,xstartr(3)-lvlhalo,xendr(3)+lvlhalo)
    end if

    return
end subroutine InitSalVariables