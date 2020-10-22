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

      do ic=1,nzm
       do jc=1,nym
        temptp(jc,ic)=-0.5d0
        tempbp(jc,ic)=0.5d0
       enddo
      enddo

      return
      end
!
      subroutine SetSalBCs
      use param
      implicit none
      integer :: ic,jc

      do ic=1,nzmr
       do jc=1,nymr
        saltp(jc,ic)=1.d0
        salbp(jc,ic)=0.d0
       enddo
      enddo

      return
      end

