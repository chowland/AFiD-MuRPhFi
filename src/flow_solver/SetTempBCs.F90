!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: SetTempBCs.F90                                 !
!    CONTAINS: subroutine SetTempBCs                      !
!                                                         ! 
!    PURPOSE: Initialization routine. Calcuates the       !
!     temperature boundary conditions at the plates       !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine SetTempBCs
      use param
      use decomp_2d
      implicit none
      integer :: ic,jc

      if (rayt>=0) then ! unstable T gradient
        do ic=xstart(3),xend(3)
          do jc=xstart(2),xend(2)
            temptp(1,jc,ic)=-0.5d0
            tempbp(1,jc,ic)=0.5d0
          end do
        end do
      else              ! stable T gradient
        do ic=xstart(3),xend(3)
          do jc=xstart(2),xend(2)
            temptp(1,jc,ic)=0.5d0
            tempbp(1,jc,ic)=-0.5d0
          end do
        end do
      end if

      if (phasefield) then
        do ic=xstart(3),xend(3)
          do jc=xstart(2),xend(2)
            temptp(1,jc,ic) = 0.d0
            tempbp(1,jc,ic) = 1.d0
          end do
        end do
      end if

      call update_halo(temptp,lvlhalo)
      call update_halo(tempbp,lvlhalo)

      return
      end
!
      ! subroutine SetSalBCs
      ! use param
      ! use decomp_2d
      ! implicit none
      ! integer :: ic,jc

      ! if (rays>=0) then ! unstable S gradient
      !   do ic=xstartr(3),xendr(3)
      !     do jc=xstartr(2),xendr(2)
      !       saltp(1,jc,ic)=0.5d0
      !       salbp(1,jc,ic)=-0.5d0
      !     enddo
      !   enddo
      ! else              ! stable S gradient
      !   do ic=xstartr(3),xendr(3)
      !     do jc=xstartr(2),xendr(2)
      !       saltp(1,jc,ic)=-0.5d0
      !       salbp(1,jc,ic)=0.5d0
      !     enddo
      !   enddo
      ! end if
      ! !CJH Add halo for interpolation routine
      ! call update_halo(saltp,lvlhalo)
      ! call update_halo(salbp,lvlhalo)

      ! return
      ! end
