!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: TimeMarcher.F90                                !
!    CONTAINS: subroutine TimeMarcher                     !
!                                                         ! 
!    PURPOSE: Main time integrating routine, which calls  !
!     other subroutines for calculating the Navier-Stokes !
!     equations and advancing velocity and temperature in !
!     time                                                !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine TimeMarcher
      use param
      use local_arrays
      use mgrd_arrays, only: vxr,vyr,vzr,salc,sal
      use mpih
      use decomp_2d
      implicit none
      integer :: ns
      integer :: j,k,i

      beta=dt/ren*0.5d0

      do ns=1,nsst                                                 

!RO     Coefficients for time marching integration (alpha, gamma, rho)

        al=alm(ns)
        ga=gam(ns)
        ro=rom(ns)

        if (melt) call UpdateBCs

        call ExplicitTermsVX
        call ExplicitTermsVY
        call ExplicitTermsVZ
        call ExplicitTermsTemp
        !if(flag_second_scalar) call ExplicitTermsSal !Refined
        call ExplicitTermsSal !Refined

        call ImplicitAndUpdateVX
        call ImplicitAndUpdateVY
        call ImplicitAndUpdateVZ
        call ImplicitAndUpdateTemp
        !if(flag_second_scalar) call ImplicitAndUpdateSal !Refined
        call ImplicitAndUpdateSal !Refined

        call update_halo(vy,lvlhalo)
        call update_halo(vz,lvlhalo)

        call CalcLocalDivergence
        call SolvePressureCorrection

!EP this copy can be avoided by changing transpose_x_to_y_real and
!transpose_y_to_x_real so these routines can handle arrays with
!halo. This copy is a defacto array temporary. Using inferred size
!arrays in the transpose calls results in 5 more of these, and more
!memory usage.  Time spent on this copy is 0.1% for 65^3 grid.

        do i=xstart(3),xend(3)
          do j=xstart(2),xend(2)
            do k=1,nxm
              dphhalo(k,j,i) = dph(k,j,i)
            enddo
          enddo
        enddo

        call update_halo(dphhalo,lvlhalo)

        call CorrectVelocity
        call CorrectPressure

        call update_halo(vx,lvlhalo)
        call update_halo(vy,lvlhalo)
        call update_halo(vz,lvlhalo)
        call update_halo(pr,lvlhalo)
        call update_halo(temp,lvlhalo)
        call update_halo(sal,lvlhalo)

        call InterpVelMgrd !Vel from base mesh to refined mesh
        call update_halo(vxr,lvlhalo)
        call update_halo(vyr,lvlhalo)
        call update_halo(vzr,lvlhalo)
        call InterpSalMgrd !Sal from refined mesh to base mesh
        call update_halo(salc,lvlhalo)

        enddo


      return
      end
