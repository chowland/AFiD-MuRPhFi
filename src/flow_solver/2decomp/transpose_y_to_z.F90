! -*- mode: f90 -*-
!=======================================================================
! This is part of the 2DECOMP&FFT library
! 
! 2DECOMP&FFT is a software framework for general-purpose 2D (pencil) 
! decomposition. It also implements a highly scalable distributed
! three-dimensional Fast Fourier Transform (FFT).
!
! Copyright (C) 2009-2011 Ning Li, the Numerical Algorithms Group (NAG)
!
!=======================================================================

! This file contains the routines that transpose data from Y to Z pencil

subroutine transpose_y_to_z_real(src, dst, opt_decomp)

  implicit none

  real(mytype), dimension(:,:,:), intent(IN) :: src
  real(mytype), dimension(:,:,:), intent(OUT) :: dst
  TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp

  TYPE(DECOMP_INFO) :: decomp

  integer :: s1,s2,s3,d1,d2,d3
  integer :: ierror

  if (present(opt_decomp)) then
     decomp = opt_decomp
  else
     decomp = decomp_main
  end if

  s1 = SIZE(src,1)
  s2 = SIZE(src,2)
  s3 = SIZE(src,3)
  d1 = SIZE(dst,1)
  d2 = SIZE(dst,2)
  d3 = SIZE(dst,3)

  ! rearrange source array as send buffer
  call mem_split_yz_real(src, s1, s2, s3, work1_r, dims(2), &
       decomp%y2dist, decomp)

  if (decomp%even) then
     call MPI_ALLTOALL(work1_r, decomp%y2count, &
          real_type, dst, decomp%z2count, &
          real_type, DECOMP_2D_COMM_ROW, ierror)
  else
     call MPI_ALLTOALL(work1_r, decomp%y2count, &
          real_type, work2_r, decomp%z2count, &
          real_type, DECOMP_2D_COMM_ROW, ierror)
  end if

  ! rearrange receive buffer
  if (.not. decomp%even) then
     call mem_merge_yz_real(work2_r, d1, d2, d3, dst, dims(2), &
          decomp%z2dist, decomp)
  end if

  return
end subroutine transpose_y_to_z_real


subroutine transpose_y_to_z_complex(src, dst, opt_decomp)

  implicit none

  complex(mytype), dimension(:,:,:), intent(IN) :: src
  complex(mytype), dimension(:,:,:), intent(OUT) :: dst
  TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp

  TYPE(DECOMP_INFO) :: decomp

  integer :: s1,s2,s3,d1,d2,d3
  integer :: ierror

  if (present(opt_decomp)) then
     decomp = opt_decomp
  else
     decomp = decomp_main
  end if

  s1 = SIZE(src,1)
  s2 = SIZE(src,2)
  s3 = SIZE(src,3)
  d1 = SIZE(dst,1)
  d2 = SIZE(dst,2)
  d3 = SIZE(dst,3)

  ! rearrange source array as send buffer
  call mem_split_yz_complex(src, s1, s2, s3, work1_c, dims(2), &
       decomp%y2dist, decomp)

  if (decomp%even) then
     call MPI_ALLTOALL(work1_c, decomp%y2count, &
          complex_type, dst, decomp%z2count, &
          complex_type, DECOMP_2D_COMM_ROW, ierror)
  else
     call MPI_ALLTOALL(work1_c, decomp%y2count, &
          complex_type, work2_c, decomp%z2count, &
          complex_type, DECOMP_2D_COMM_ROW, ierror)
  end if

  ! rearrange receive buffer
  if (.not. decomp%even) then
     call mem_merge_yz_complex(work2_c, d1, d2, d3, dst, dims(2), &
          decomp%z2dist, decomp)
  end if

  return
end subroutine transpose_y_to_z_complex


! pack/unpack ALLTOALL(V) buffers

subroutine mem_split_yz_real(in,n1,n2,n3,out,iproc,dist,decomp)

  implicit none

  integer, intent(IN) :: n1,n2,n3
  real(mytype), dimension(n1,n2,n3), intent(IN) :: in
  real(mytype), dimension(*), intent(OUT) :: out
  integer, intent(IN) :: iproc
  integer, dimension(0:iproc-1), intent(IN) :: dist
  TYPE(DECOMP_INFO), intent(IN) :: decomp

  integer :: i,j,k, m,i1,i2,pos

  do m=0,iproc-1
     if (m==0) then 
        i1 = 1
        i2 = dist(0)
     else
        i1 = i2+1
        i2 = i1+dist(m)-1
     end if

     pos = m * decomp%y2count + 1

     do k=1,n3
        do j=i1,i2
           do i=1,n1
              out(pos) = in(i,j,k)
              pos = pos + 1
           end do
        end do
     end do
  end do

  return
end subroutine mem_split_yz_real


subroutine mem_split_yz_complex(in,n1,n2,n3,out,iproc,dist,decomp)

  implicit none

  integer, intent(IN) :: n1,n2,n3
  complex(mytype), dimension(n1,n2,n3), intent(IN) :: in
  complex(mytype), dimension(*), intent(OUT) :: out
  integer, intent(IN) :: iproc
  integer, dimension(0:iproc-1), intent(IN) :: dist
  TYPE(DECOMP_INFO), intent(IN) :: decomp

  integer :: i,j,k, m,i1,i2,pos

  do m=0,iproc-1
     if (m==0) then 
        i1 = 1
        i2 = dist(0)
     else
        i1 = i2+1
        i2 = i1+dist(m)-1
     end if

     pos = m * decomp%y2count + 1

     do k=1,n3
        do j=i1,i2
           do i=1,n1
              out(pos) = in(i,j,k)
              pos = pos + 1
           end do
        end do
     end do
  end do

  return
end subroutine mem_split_yz_complex


subroutine mem_merge_yz_real(in,n1,n2,n3,out,iproc,dist,decomp)

  implicit none

  integer, intent(IN) :: n1,n2,n3
  real(mytype), dimension(*), intent(IN) :: in
  real(mytype), dimension(n1,n2,n3), intent(OUT) :: out
  integer, intent(IN) :: iproc
  integer, dimension(0:iproc-1), intent(IN) :: dist
  TYPE(DECOMP_INFO), intent(IN) :: decomp

  integer :: i,j,k, m,i1,i2, pos

  do m=0,iproc-1
     if (m==0) then
        i1 = 1
        i2 = dist(0)
     else
        i1 = i2+1
        i2 = i1+dist(m)-1
     end if

     pos = m * decomp%z2count + 1

     do k=i1,i2
        do j=1,n2
           do i=1,n1
              out(i,j,k) = in(pos)
              pos = pos + 1
           end do
        end do
     end do
  end do

  return
end subroutine mem_merge_yz_real


subroutine mem_merge_yz_complex(in,n1,n2,n3,out,iproc,dist,decomp)

  implicit none

  integer, intent(IN) :: n1,n2,n3
  complex(mytype), dimension(*), intent(IN) :: in
  complex(mytype), dimension(n1,n2,n3), intent(OUT) :: out
  integer, intent(IN) :: iproc
  integer, dimension(0:iproc-1), intent(IN) :: dist
  TYPE(DECOMP_INFO), intent(IN) :: decomp

  integer :: i,j,k, m,i1,i2, pos

  do m=0,iproc-1
     if (m==0) then
        i1 = 1
        i2 = dist(0)
     else
        i1 = i2+1
        i2 = i1+dist(m)-1
     end if

     pos = m * decomp%z2count + 1

     do k=i1,i2
        do j=1,n2
           do i=1,n1
              out(i,j,k) = in(pos)
              pos = pos + 1
           end do
        end do
     end do
  end do

  return
end subroutine mem_merge_yz_complex
