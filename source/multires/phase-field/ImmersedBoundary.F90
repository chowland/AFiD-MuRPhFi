!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: ImmersedBoundary.F90                           !
!    CONTAINS: subroutine ImmersedBoundary                !
!                                                         ! 
!    PURPOSE: Enforce zero velocity in the solid,         !
!     defined by phi>0.5                                  !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ImmersedBoundary
    
    use param
    use local_arrays, only: vx, vy, vz
    use mgrd_arrays, only: phic
    use decomp_2d, only: xstart,xend
    
    implicit none
    
    integer :: k,j,i

    do i=xstart(3),xend(3)
        do j=xstart(2),xend(2)
            do k=1,nxm
                if (phic(k,j,i) .gt. 0.5) then
                    vx(k,j,i) = 0.0
                    vy(k,j,i) = 0.0
                    vz(k,j,i) = 0.0
                end if
            end do
        end do
    end do

    return

end subroutine ImmersedBoundary