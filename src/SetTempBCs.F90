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
      implicit none
      integer :: ic,jc

      if (rayt>=0) then ! unstable T gradient
        do ic=1,nzm
          do jc=1,nym
            temptp(jc,ic)=-0.5d0
            tempbp(jc,ic)=0.5d0
          end do
        end do
      else              ! stable T gradient
        do ic=1,nzm
          do jc=1,nym
            temptp(jc,ic)=0.5d0
            tempbp(jc,ic)=-0.5d0
          end do
        end do
      end if


      return
      end
!
      subroutine SetSalBCs
      use param
      implicit none
      integer :: ic,jc

      if (rays>=0) then ! unstable S gradient
        do ic=1,nzmr
          do jc=1,nymr
            saltp(jc,ic)=0.5d0
            salbp(jc,ic)=-0.5d0
          enddo
        enddo
      else              ! stable S gradient
        do ic=1,nzmr
          do jc=1,nymr
            saltp(jc,ic)=-0.5d0
            salbp(jc,ic)=0.5d0
          enddo
        enddo
      end if
      !CJH Add halo for interpolation routine
      saltp(0 ,:) = saltp(nymr,:)
      saltp(-1,:) = saltp(nymr-1,:)
      saltp(:, 0) = saltp(:,nzmr)
      saltp(:,-1) = saltp(:,nzmr-1)
      saltp(nyr  ,:) = saltp(1,:)
      saltp(nyr+1,:) = saltp(2,:)
      saltp(:,nzr  ) = saltp(:,1)
      saltp(:,nzr+1) = saltp(:,2)
      salbp(0 ,:) = salbp(nymr,:)
      salbp(-1,:) = salbp(nymr-1,:)
      salbp(:, 0) = salbp(:,nzmr)
      salbp(:,-1) = salbp(:,nzmr-1)
      salbp(nyr  ,:) = salbp(1,:)
      salbp(nyr+1,:) = salbp(2,:)
      salbp(:,nzr  ) = salbp(:,1)
      salbp(:,nzr+1) = salbp(:,2)

      return
      end

