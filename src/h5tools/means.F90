module means
    use mpih
    use decomp_2d
    use param, only: nym, nymr, nzm, nzmr
    implicit none

contains

subroutine ymean(var, plane, comm, yrank)
    real, intent(in) :: var(:,xstart(2):,xstart(3):)
    real, intent(out) :: plane(:,xstart(3):)
    integer, intent(in) :: comm, yrank

    integer :: i, j, k, bufsize, n1

    plane(:,:) = 0.0

    n1 = size(var, dim=1)
    do i=xstart(3),xend(3)
        do j=xstart(2),xend(2)
            do k=1,n1
                plane(k,i) = plane(k,i) + var(k,j,i)
            end do
        end do
    end do

    bufsize = n1*xsize(3)
    if (yrank==0) then
        call MPI_REDUCE( &
            MPI_IN_PLACE, &
            plane(1:n1,xstart(3):xend(3)), &
            bufsize, MDP, MPI_SUM, 0, comm, ierr &
        )
    else
        call MPI_REDUCE( &
            plane(1:n1,xstart(3):xend(3)), &
            plane(1:n1,xstart(3):xend(3)), &
            bufsize, MDP, MPI_SUM, 0, comm, ierr &
        )
    end if

    plane = plane/nym

end subroutine ymean

subroutine zmean(var, plane, comm, zrank)
    real, intent(in) :: var(:,xstart(2):,xstart(3):)
    real, intent(out) :: plane(:,xstart(2):)
    integer, intent(in) :: comm, zrank

    integer :: i, j, k, bufsize, n1

    plane(:,:) = 0.0

    n1 = size(var, dim=1)
    do i=xstart(3),xend(3)
        do j=xstart(2),xend(2)
            do k=1,n1
                plane(k,j) = plane(k,j) + var(k,j,i)
            end do
        end do
    end do

    bufsize = n1*xsize(2)
    if (zrank==0) then
        call MPI_REDUCE( &
            MPI_IN_PLACE, &
            plane(1:n1,xstart(2):xend(2)), &
            bufsize, MDP, MPI_SUM, 0, comm, ierr &
        )
    else
        call MPI_REDUCE( &
            plane(1:n1,xstart(2):xend(2)), &
            plane(1:n1,xstart(2):xend(2)), &
            bufsize, MDP, MPI_SUM, 0, comm, ierr &
        )
    end if

    plane = plane/nzm

end subroutine zmean

subroutine ymeanr(var, plane, comm, yrank)
    real, intent(in) :: var(:,xstartr(2):,xstartr(3):)
    real, intent(out) :: plane(:,xstartr(3):)
    integer, intent(in) :: comm, yrank

    integer :: i, j, k, bufsize, n1
    
    plane(:,:) = 0.0

    n1 = size(var, dim=1)
    do i=xstartr(3),xendr(3)
        do j=xstartr(2),xendr(2)
            do k=1,n1
                plane(k,i) = plane(k,i) + var(k,j,i)
            end do
        end do
    end do

    bufsize = n1*xsizer(3)
    if (yrank==0) then
        call MPI_REDUCE( &
            MPI_IN_PLACE, &
            plane(1:n1,xstartr(3):xendr(3)), &
            bufsize, MDP, MPI_SUM, 0, comm, ierr &
        )
    else
        call MPI_REDUCE( &
            plane(1:n1,xstartr(3):xendr(3)), &
            plane(1:n1,xstartr(3):xendr(3)), &
            bufsize, MDP, MPI_SUM, 0, comm, ierr &
        )
    end if

    plane = plane/nymr

end subroutine ymeanr

subroutine zmeanr(var, plane, comm, zrank)
    real, intent(in) :: var(:,xstartr(2):,xstartr(3):)
    real, intent(out) :: plane(:,xstartr(2):)
    integer, intent(in) :: comm, zrank

    integer :: i, j, k, bufsize, n1
    
    plane(:,:) = 0.0

    n1 = size(var, dim=1)
    do i=xstartr(3),xendr(3)
        do j=xstartr(2),xendr(2)
            do k=1,n1
                plane(k,j) = plane(k,j) + var(k,j,i)
            end do
        end do
    end do

    bufsize = n1*xsizer(2)
    if (zrank==0) then
        call MPI_REDUCE( &
            MPI_IN_PLACE, &
            plane(1:n1,xstartr(2):xendr(2)), &
            bufsize, MDP, MPI_SUM, 0, comm, ierr &
        )
    else
        call MPI_REDUCE( &
            plane(1:n1,xstartr(2):xendr(2)), &
            plane(1:n1,xstartr(2):xendr(2)), &
            bufsize, MDP, MPI_SUM, 0, comm, ierr &
        )
    end if

    plane = plane/nzmr

end subroutine zmeanr

end module means