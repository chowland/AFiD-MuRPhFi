!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: InitVariables.F90                              !
!    CONTAINS: subroutine InitVariables                   !
!                                                         ! 
!    PURPOSE: Initialization routine. Sets to zero all    !
!     variables used in the code                          !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine InitVariables
      use param
      use local_arrays
      ! use mgrd_arrays
      use stat_arrays
      use decomp_2d
      use AuxiliaryRoutines
      implicit none
      
!-------------------------------------------------
! Arrays for grid making
!-------------------------------------------------

      call AllocateReal1DArray(zc,1,nz)
      call AllocateReal1DArray(zm,0,nz+1)
      call AllocateReal1DArray(ak1,1,nz)
      call AllocateReal1DArray(ao,1,nz)

      call AllocateReal1DArray(yc,1,ny)
      call AllocateReal1DArray(ym,0,ny+1)
      call AllocateReal1DArray(ak2,1,ny)
      call AllocateReal1DArray(ap,1,ny)

      call AllocateReal1DArray(xc,1,nx)
      call AllocateReal1DArray(xm,1,nx)
      call AllocateReal1DArray(g3rc,1,nx)
      call AllocateReal1DArray(g3rm,1,nx)

      !CJH New terms for modified derivatives
      call AllocateReal1DArray(d3xc,1,nx)
      call AllocateReal1DArray(d3xm,1,nx)

      call AllocateReal1DArray(udx3c,1,nx)
      call AllocateReal1DArray(udx3m,1,nx)

      call AllocateReal1DArray(ap3ck,1,nx)
      call AllocateReal1DArray(ac3ck,1,nx)
      call AllocateReal1DArray(am3ck,1,nx)

      call AllocateReal1DArray(ap3sk,1,nx)
      call AllocateReal1DArray(ac3sk,1,nx)
      call AllocateReal1DArray(am3sk,1,nx)

      call AllocateReal1DArray(ap3ssk,1,nx)
      call AllocateReal1DArray(ac3ssk,1,nx)
      call AllocateReal1DArray(am3ssk,1,nx)

      call AllocateReal1DArray(amphk,1,nx)
      call AllocateReal1DArray(acphk,1,nx)
      call AllocateReal1DArray(apphk,1,nx)
 
      call AllocateInt1dArray(kmc,1,nx)
      call AllocateInt1dArray(kpc,1,nx)
      call AllocateInt1dArray(kmv,1,nx)
      call AllocateInt1dArray(kpv,1,nx)

!-------------------------------------------------
! Arrays for mgrd grid making
!-------------------------------------------------

      ! call AllocateReal1DArray(zcr,1,nzr)
      ! call AllocateReal1DArray(zmr,0,nzr+1) !CS mgrd

      ! call AllocateReal1DArray(ycr,1,nyr)
      ! call AllocateReal1DArray(ymr,0,nyr+1) !CS mgrd

      ! call AllocateReal1DArray(xcr,1,nxr)
      ! call AllocateReal1DArray(xmr,0,nxr) !CS mgrd

      ! call AllocateReal1DArray(g3rcr,1,nxr)
      ! call AllocateReal1DArray(g3rmr,1,nxr)

      ! !CJH New terms for modified derivatives
      ! call AllocateReal1DArray(d3xcr,1,nxr)
      ! call AllocateReal1DArray(d3xmr,1,nxr)

      ! call AllocateReal1DArray(udx3cr,1,nxr)
      ! call AllocateReal1DArray(udx3mr,1,nxr)

      ! call AllocateReal1DArray(ap3ckr,1,nxr)
      ! call AllocateReal1DArray(ac3ckr,1,nxr)
      ! call AllocateReal1DArray(am3ckr,1,nxr)

      ! call AllocateReal1DArray(ap3sskr,1,nxr)
      ! call AllocateReal1DArray(ac3sskr,1,nxr)
      ! call AllocateReal1DArray(am3sskr,1,nxr)

      ! call AllocateInt1dArray(kmcr,1,nxr)
      ! call AllocateInt1dArray(kpcr,1,nxr)
      ! call AllocateInt1dArray(kmvr,1,nxr)
      ! call AllocateInt1dArray(kpvr,1,nxr)

      ! call AllocateInt1dArray(irangs,0,nx)
      ! call AllocateInt1dArray(jrangs,0,ny)
      ! call AllocateInt1dArray(krangs,0,nz)

      ! call AllocateInt1dArray(irangc,0,nx)
      ! call AllocateInt1dArray(jrangc,0,ny)
      ! call AllocateInt1dArray(krangc,0,nz)

      ! call AllocateReal2DArray(cxvx,1,4,0,nxr)
      ! call AllocateReal2DArray(cxvy,1,4,0,nxr)
      ! call AllocateReal2DArray(cxvz,1,4,0,nxr)

      ! call AllocateReal2DArray(cyvx,1,4,0,nyr)
      ! call AllocateReal2DArray(cyvy,1,4,0,nyr)
      ! call AllocateReal2DArray(cyvz,1,4,0,nyr)

      ! call AllocateReal2DArray(czvx,1,4,0,nzr)
      ! call AllocateReal2DArray(czvy,1,4,0,nzr)
      ! call AllocateReal2DArray(czvz,1,4,0,nzr)

      ! call AllocateReal2DArray(cxrs,1,4,0,nxr)
      ! call AllocateReal2DArray(cyrs,1,4,0,nyr)
      ! call AllocateReal2DArray(czrs,1,4,0,nzr)

      ! call AllocateReal2DArray(cxsalc,1,4,0,nx)
      ! call AllocateReal2DArray(cysalc,1,4,0,ny)
      ! call AllocateReal2DArray(czsalc,1,4,0,nz)

!-------------------------------------------------
! Arrays for temperature boundary conditions    
!-------------------------------------------------

      call AllocateReal3DArray(tempbp,1,1,xstart(2)-lvlhalo,xend(2)+lvlhalo,xstart(3)-lvlhalo,xend(3)+lvlhalo)
      call AllocateReal3DArray(temptp,1,1,xstart(2)-lvlhalo,xend(2)+lvlhalo,xstart(3)-lvlhalo,xend(3)+lvlhalo)

      !-- For salinity
      ! call AllocateReal3DArray(salbp,1,1,xstartr(2)-lvlhalo,xendr(2)+lvlhalo,xstartr(3)-lvlhalo,xendr(3)+lvlhalo)
      ! call AllocateReal3DArray(saltp,1,1,xstartr(2)-lvlhalo,xendr(2)+lvlhalo,xstartr(3)-lvlhalo,xendr(3)+lvlhalo)

!-------------------------------------------------
! Arrays for statistics    
!-------------------------------------------------

      ! if (statcal) then
      !  call AllocateReal1DArray(vx_me,1,nxm)
      !  call AllocateReal1DArray(vy_me,1,nxm)
      !  call AllocateReal1DArray(vz_me,1,nxm)
 
      !  call AllocateReal1DArray(vx_rms,1,nxm)
      !  call AllocateReal1DArray(vy_rms,1,nxm)
      !  call AllocateReal1DArray(vz_rms,1,nxm)
 
      !  call AllocateReal1DArray(vx_me_buf,1,nxm)
      !  call AllocateReal1DArray(vy_me_buf,1,nxm)
      !  call AllocateReal1DArray(vz_me_buf,1,nxm)
 
      !  call AllocateReal1DArray(vx_msq_buf,1,nxm)
      !  call AllocateReal1DArray(vy_msq_buf,1,nxm)
      !  call AllocateReal1DArray(vz_msq_buf,1,nxm)

      !  call AllocateReal1DArray(temp_me,1,nxm)
      !  call AllocateReal1DArray(temp_rms,1,nxm)
      !  call AllocateReal1DArray(tempvx_me,1,nxm)

      !  call AllocateReal1DArray(vxvy_corr,1,nxm)
      !  if (disscal) then
      !   call AllocateReal1DArray(disste,1,nxm)
      !   call AllocateReal1DArray(dissth,1,nxm)
      !  end if
      ! end if

      !-------------------------------------------------
      ! Arrays with ghost cells
      !-------------------------------------------------
      call AllocateReal3DArray(vy,1,nx,xstart(2)-lvlhalo,xend(2)+lvlhalo,xstart(3)-lvlhalo,xend(3)+lvlhalo)
      call AllocateReal3DArray(vz,1,nx,xstart(2)-lvlhalo,xend(2)+lvlhalo,xstart(3)-lvlhalo,xend(3)+lvlhalo)
      call AllocateReal3DArray(vx,1,nx,xstart(2)-lvlhalo,xend(2)+lvlhalo,xstart(3)-lvlhalo,xend(3)+lvlhalo)
      call AllocateReal3DArray(pr,1,nx,xstart(2)-lvlhalo,xend(2)+lvlhalo,xstart(3)-lvlhalo,xend(3)+lvlhalo)
      call AllocateReal3DArray(temp,1,nx,xstart(2)-lvlhalo,xend(2)+lvlhalo,xstart(3)-lvlhalo,xend(3)+lvlhalo)
      call AllocateReal3DArray(dphhalo,1,nxm,xstart(2)-lvlhalo,xend(2)+lvlhalo,xstart(3)-lvlhalo,xend(3)+lvlhalo)

      !-- For salinity
      ! call AllocateReal3DArray(sal,1,nxr,xstartr(2)-lvlhalo,xendr(2)+lvlhalo,xstartr(3)-lvlhalo,xendr(3)+lvlhalo)

      !-- For Q-criteria
      call AllocateReal3DArray(qtens,1,nx,xstart(2)-lvlhalo,xend(2)+lvlhalo,xstart(3)-lvlhalo,xend(3)+lvlhalo)

      !-----------------------------------------------
      ! Arrays without ghost cells
      !-----------------------------------------------
      call AllocateReal3DArray(rhs,1,nx,xstart(2),xend(2),xstart(3),xend(3))
      call AllocateReal3DArray(dph,1,nxm,xstart(2),xend(2),xstart(3),xend(3))
      call AllocateReal3DArray(dq,1,nx,xstart(2),xend(2),xstart(3),xend(3))
      call AllocateReal3DArray(qcap,1,nx,xstart(2),xend(2),xstart(3),xend(3))
      call AllocateReal3DArray(rux,1,nx,xstart(2),xend(2),xstart(3),xend(3))
      call AllocateReal3DArray(ruy,1,nx,xstart(2),xend(2),xstart(3),xend(3))
      call AllocateReal3DArray(ruz,1,nx,xstart(2),xend(2),xstart(3),xend(3))
      call AllocateReal3DArray(hro,1,nx,xstart(2),xend(2),xstart(3),xend(3))
      call AllocateReal3DArray(rutemp,1,nx,xstart(2),xend(2),xstart(3),xend(3))

      !-- For salinity
      ! call AllocateReal3DArray(rhsr,1,nxr,xstartr(2),xendr(2),xstartr(3),xendr(3))
      ! call AllocateReal3DArray(rusal,1,nxr,xstartr(2),xendr(2),xstartr(3),xendr(3))
      ! call AllocateReal3DArray(hsal,1,nxr,xstartr(2),xendr(2),xstartr(3),xendr(3))

      return 
      end   

! --------------------------------
!       subroutine InitMgrdVariables
!       use param
!       use mgrd_arrays
!       use decomp_2d
!       use AuxiliaryRoutines
!       implicit none

!       call AllocateReal3DArray(vyr,1,nxr,xstartr(2)-lvlhalo,xendr(2)+lvlhalo,xstartr(3)-lvlhalo,xendr(3)+lvlhalo)
!       call AllocateReal3DArray(vzr,1,nxr,xstartr(2)-lvlhalo,xendr(2)+lvlhalo,xstartr(3)-lvlhalo,xendr(3)+lvlhalo)
!       call AllocateReal3DArray(vxr,1,nxr,xstartr(2)-lvlhalo,xendr(2)+lvlhalo,xstartr(3)-lvlhalo,xendr(3)+lvlhalo)

!       call AllocateReal3DArray(tpdv,-1,nx+1,xstart(2)-2,xend(2)+2,xstart(3)-2,xend(3)+2)
! !CS   For tpdvr, larger array needed to prevent memory overflow in InterpVelMgrd
!       call AllocateReal3DArray(tpdvr,1,nxmr,xstartr(2)-2,xendr(2)+2,xstartr(3)-2,xendr(3)+2)

!       call AllocateReal3DArray(salc,1,nx,xstart(2)-lvlhalo,xend(2)+lvlhalo,xstart(3)-lvlhalo,xend(3)+lvlhalo)
!       !CJH Needed for melt boundary condition
!       if (melt) then
!         call AllocateReal3DArray(tempr,1,1,xstartr(2)-lvlhalo,xendr(2)+lvlhalo,xstartr(3)-lvlhalo,xendr(3)+lvlhalo)
!       end if

!       return 
!       end   
