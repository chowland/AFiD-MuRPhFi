!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: ResetLogs.F90                                  !
!    CONTAINS: subroutine ResetLogs                       !
!                                                         ! 
!    PURPOSE: Initialization routine. Reset all log files !
!     if necessary                                        !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine ResetLogs
      use param
      implicit none

      if(resetlogstime) then    

!EP    nusse.out  in GlobalQuantities.F
      open(95,file='outputdir/nu_vol.out',status='unknown',access='sequential', &
        position='append')
      close(95,status='delete')

!EP    nusse2.out  in CalculatePlateNu.F
       open(97,file="outputdir/nu_plate.out",status='unknown',access='sequential', &
        position='append')
      close(97,status='delete')

!EP   nusse3.out in CalcDissipationNu.F
      if (disscal) then
      open(92,file='outputdir/nu_diss.out',status='unknown',access='sequential', &
       position='append')
      close(92,status='delete')
      end if

!EP   rms_vel.out in GlobalQuantities.F
       open(94,file='outputdir/rms_vel.out',status='unknown',position='append', &
        access='sequential')
      close(94,status='delete')

!CS   cfl.out in main.F90                                                    
       open(96,file='outputdir/cfl.out',status='unknown',position='append', &
        access='sequential')                                          
      close(96,status='delete')

!CS    cf.out  in CalculatePlateCf.F
       open(98,file="outputdir/cf_plate.out",status='unknown',access='sequential', &
        position='append')
      close(98,status='delete')

      endif

      return
      end   
      

