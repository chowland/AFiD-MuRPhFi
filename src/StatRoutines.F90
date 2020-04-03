!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: StatRoutines.F90                               !
!    CONTAINS: subroutine CalcStats,WriteStats            !
!                         CalcSpecStats,WriteSpecStats    !
!                                                         ! 
!    PURPOSE: Calculates and writes out statistics for    !
!     the flow field. All quantities are averaged in the  !
!     two horizontal (homogeneous) directions             !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine CalcStats
      use param
      use local_arrays, only: vz,vy,vx,temp
      use decomp_2d, only: xstart,xend
      use stat_arrays
      use mpih
      implicit none
      real :: usnzm,usnym
      integer :: i,j,k
      character(5) :: ipfi
      character(70) :: filnam1,filnam2,filnam3,filnam4,filnam5,filnam6

      vx_me_buf(:)=0.d0; vy_me_buf(:)=0.d0; vz_me_buf(:)=0.d0
      vx_msq_buf(:)=0.d0; vy_msq_buf(:)=0.d0; vz_msq_buf(:)=0.d0

      nstatsamples = nstatsamples + 1

      usnym = 1.0/nym
      usnzm = 1.0/nzm

      do i=xstart(3),xend(3)
        do j=xstart(2),xend(2)
          do k=1,nxm
               vx_me_buf(k) = vx_me_buf(k) + vx(k,j,i)*usnzm*usnym
               vy_me_buf(k) = vy_me_buf(k) + vy(k,j,i)*usnzm*usnym
               vz_me_buf(k) = vz_me_buf(k) + vz(k,j,i)*usnzm*usnym

               vx_me(k) = vx_me(k) + vx(k,j,i)*usnzm*usnym
               vy_me(k) = vy_me(k) + vy(k,j,i)*usnzm*usnym
               vz_me(k) = vz_me(k) + vz(k,j,i)*usnzm*usnym

               temp_me(k) = temp_me(k) + temp(k,j,i)*usnzm*usnym

               vx_rms(k) = vx_rms(k) + vx(k,j,i)**2*usnzm*usnym
               vy_rms(k) = vy_rms(k) + vy(k,j,i)**2*usnzm*usnym
               vz_rms(k) = vz_rms(k) + vz(k,j,i)**2*usnzm*usnym

               temp_rms(k) = temp_rms(k) +  &
     &                               temp(k,j,i)**2*usnzm*usnym
               tempvx_me(k) = tempvx_me(k) +  &
     &                         temp(k,j,i)*vx(k,j,i)*usnzm*usnym
            end do
         end do
      end do

      call MpiAllSumReal1D(vx_me_buf,nxm)
      call MpiAllSumReal1D(vy_me_buf,nxm)
      call MpiAllSumReal1D(vz_me_buf,nxm)

      do i=xstart(3),xend(3)
        do j=xstart(2),xend(2)
          do k=1,nxm-1
               vxvy_corr(k) = vxvy_corr(k) + (vy(k,j,i)-vy_me_buf(k))* &
     &                                       ((vx(k,j,i)+vx(k+1,j,i))- &
     &                                       (vx_me_buf(k)+vx_me_buf(k+1)))* &
     &                                       0.5*usnzm*usnym
               vx_msq_buf(k) = vx_msq_buf(k) + (vx(k,j,i)-vx_me_buf(k))**2*usnzm*usnym
               vy_msq_buf(k) = vy_msq_buf(k) + (vy(k,j,i)-vy_me_buf(k))**2*usnzm*usnym
               vz_msq_buf(k) = vz_msq_buf(k) + (vz(k,j,i)-vz_me_buf(k))**2*usnzm*usnym
            end do
         end do
      end do

      !-- Reduce and write runtime stats
      call MpiSumReal1D(vx_msq_buf,nxm)
      call MpiSumReal1D(vy_msq_buf,nxm)
      call MpiSumReal1D(vz_msq_buf,nxm)

      write(ipfi,'(i5.5)')nint(time/tframe)
      filnam1='outputdir/stst/vx_avg_'//ipfi//'.out'
      filnam2='outputdir/stst/vy_avg_'//ipfi//'.out'
      filnam3='outputdir/stst/vz_avg_'//ipfi//'.out'
      filnam4='outputdir/stst/vxvx_msq_'//ipfi//'.out'
      filnam5='outputdir/stst/vyvy_msq_'//ipfi//'.out'
      filnam6='outputdir/stst/vzvz_msq_'//ipfi//'.out'

      if(ismaster) then
         open(unit=200,file=filnam1,status='unknown')
            do k=1,nxm
               write(200,'(ES20.10)') vx_me_buf(k)
            end do
         close(200)
         open(unit=201,file=filnam2,status='unknown')
            do k=1,nxm
               write(201,'(ES20.10)') vy_me_buf(k)
            end do
         close(201)
         open(unit=202,file=filnam3,status='unknown')
            do k=1,nxm
               write(202,'(ES20.10)') vz_me_buf(k)
            end do
         close(202)
         open(unit=203,file=filnam4,status='unknown')
            do k=1,nxm
               write(203,'(ES20.10)') vx_msq_buf(k)
            end do
         close(203)
         open(unit=204,file=filnam5,status='unknown')
            do k=1,nxm
               write(204,'(ES20.10)') vy_msq_buf(k)
            end do
         close(204)
         open(unit=205,file=filnam6,status='unknown')
            do k=1,nxm
               write(205,'(ES20.10)') vz_msq_buf(k)
            end do
         close(205)
      end if

      return  
      end
!    
!***********************************************************************
      subroutine WriteStats
      use mpih
      use param
      use stat_arrays
      use hdf5

      implicit none

      integer :: nstatsamples_old

      character*30 dsetname_vxme
      character*30 dsetname_vyme
      character*30 dsetname_vzme
      character*30 dsetname_tempme

      character*30 dsetname_vxrms
      character*30 dsetname_vyrms
      character*30 dsetname_vzrms
      character*30 dsetname_temprms

      character*30 dsetname_tempvxme
      character*30 dsetname_vxvycorr
      character*30 dsetname_dissth
      character*30 dsetname_disste
      character*30 filnam,dsetname
      logical :: fexist

      filnam = trim('outputdir/stafield_master.h5')

      dsetname_vxme = trim('vx_mean')
      dsetname_vyme = trim('vy_mean')
      dsetname_vzme = trim('vz_mean')
      dsetname_tempme = trim('temp_mean')

      dsetname_vxrms = trim('vx_rms')
      dsetname_vyrms = trim('vy_rms')
      dsetname_vzrms = trim('vz_rms')
      dsetname_temprms = trim('temp_rms')

      dsetname_tempvxme = trim('tempvx_mean')

      dsetname_vxvycorr = trim('vxvy_corr')

      dsetname_dissth = trim('dissth')
      dsetname_disste = trim('disste')

      dsetname = trim('averaging_time')

      inquire(file=filnam,exist=fexist)
      if (.not.fexist) then 
        if(ismaster) write(6,*) 'Unable to read statistical files'
        if(ismaster) write(6,*) 'Restarting statistics from zero' 
        readstats=.false.
      end if
       

      if (ismaster) then
       if(readstats) then
        call HdfSerialReadIntScalar(dsetname,filnam,nstatsamples_old)
        nstatsamples = nstatsamples + nstatsamples_old
       else 
        call HdfCreateBlankFile(filnam)
       endif
      end if

      call StatReadReduceWrite(vx_me,filnam,dsetname_vxme)
      call StatReadReduceWrite(vy_me,filnam,dsetname_vyme)
      call StatReadReduceWrite(vz_me,filnam,dsetname_vzme)
      call StatReadReduceWrite(temp_me,filnam,dsetname_tempme)

      call StatReadReduceWrite(vx_rms,filnam,dsetname_vxrms)
      call StatReadReduceWrite(vy_rms,filnam,dsetname_vyrms)
      call StatReadReduceWrite(vz_rms,filnam,dsetname_vzrms)
      call StatReadReduceWrite(temp_rms,filnam,dsetname_temprms)
 
      call StatReadReduceWrite(tempvx_me,filnam,dsetname_tempvxme)

      call StatReadReduceWrite(vxvy_corr,filnam,dsetname_vxvycorr)

      if(disscal) then 
       ! call StatReadReduceWrite(dissth,filnam,dsetname_dissth)
       call StatReadReduceWrite(disste,filnam,dsetname_disste)
      end if

      if (ismaster) then

       call HdfSerialWriteIntScalar(dsetname,filnam,nstatsamples)

       dsetname = trim('X_cordin')
       call HdfSerialWriteReal1D(dsetname,filnam,xm,nxm)

       dsetname = trim('Reynolds Number')
       call HdfSerialWriteRealScalar(dsetname,filnam,ren)

       dsetname = trim('Prandtl Number')
       call HdfSerialWriteRealScalar(dsetname,filnam,pra)


      endif

      return  
      end
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine CalcGlobalStats
      use param
      use local_arrays, only: vz,vy,vx
      use decomp_2d, only: xstart,xend
      use stat_arrays
      use mpih
      implicit none
      real :: usnzm,usnym
      integer :: i,j,k

      vx_global=0.d0; vy_global=0.d0; vz_global=0.d0

      usnym = 1.0/nym
      usnzm = 1.0/nzm

      do i=xstart(3),xend(3)
        do j=xstart(2),xend(2)
          do k=1,nxm
               vx_global = vx_global + vx(k,j,i)*usnzm*usnym/udx3m(k)
               vy_global = vy_global + vy(k,j,i)*usnzm*usnym/udx3m(k)
               vz_global = vz_global + vz(k,j,i)*usnzm*usnym/udx3m(k)
            end do
         end do
      end do

      call MpiSumRealScalar(vx_global)
      call MpiSumRealScalar(vy_global)
      call MpiSumRealScalar(vz_global)

      return  
      end
!
