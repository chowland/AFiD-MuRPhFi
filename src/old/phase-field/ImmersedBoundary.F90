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
            do k=2,nxm
                if (0.5*(phic(k,j,i) + phic(k-1,j,i)) .gt. 0.5) then
                    vx(k,j,i) = 0.0
                end if
            end do
            do k=1,nxm
                if (0.5*(phic(k,j,i) + phic(k,j-1,i)) .gt. 0.5) then
                    vy(k,j,i) = 0.0
                end if
                if (0.5*(phic(k,j,i) + phic(k,j,i-1)) .gt. 0.5) then
                    vz(k,j,i) = 0.0
                end if
            end do
        end do
    end do

    return

end subroutine ImmersedBoundary