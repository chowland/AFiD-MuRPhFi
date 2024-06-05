!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
!    FILE: ReadInputFile.F90                              !
!    CONTAINS: subroutine ReadInputFile                   !
!                                                         !
!    PURPOSE: Read parameters from bou.in file            !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ReadInputFile
    use param
    use ibm_param, only: solidtype
    use afid_salinity, only: RayS, PraS, bycs, PecS, SfixN, SfixS
    use afid_phasefield, only: pf_A, pf_D, pf_eps, pf_Lambda, pf_S, pf_Tm
    implicit none
    character(len=4) :: dummy
    integer :: flagmelt
    integer :: flagMR, flagsal, flagPF
    integer :: FFscaleS

    open(unit=15,file='bou.in',status='old')
        read(15,301) dummy
        read(15,301) dummy
        read(15,*) nxm, nym, nzm, nsst
        read(15,301) dummy
        read(15,301) dummy
        read(15,301) dummy
        read(15,*) flagMR, nxmr, nymr, nzmr
        read(15,301) dummy
        read(15,301) dummy
        read(15,301) dummy
        read(15,*) flagsal, flagPF
        read(15,301) dummy
        read(15,301) dummy
        read(15,301) dummy
        read(15,*) nread, ireset
        read(15,301) dummy
        read(15,301) dummy
        read(15,301) dummy
        read(15,*) ntst, walltimemax, tmax
        read(15,301) dummy
        read(15,301) dummy
        read(15,301) dummy
        read(15,*) tout, tframe, save_3D
        read(15,301) dummy
        read(15,301) dummy
        read(15,301) dummy
        read(15,*) alx3, ylen, zlen
        read(15,301) dummy
        read(15,301) dummy
        read(15,301) dummy
        read(15,*) istr3, str3, istr3r
        read(15,301) dummy
        read(15,301) dummy
        read(15,301) dummy
        read(15,*) rayt, prat, rays, pras, FFscaleS
        read(15,301) dummy
        read(15,301) dummy
        read(15,301) dummy
        read(15,*) idtv, dt, resid, limitCFL, dtmin, dtmax
        read(15,301) dummy
        read(15,301) dummy
        read(15,301) dummy
        read(15,*) inslws, inslwn, TfixS, TfixN, SfixS, SfixN
        read(15,301) dummy
        read(15,301) dummy
        read(15,301) dummy
        read(15,*) active_T, active_S, gAxis
        read(15,301) dummy
        read(15,301) dummy
        read(15,301) dummy
        read(15,*) xplusU, xminusU, dPdy, dPdz, flagmelt
        read(15,301) dummy
        read(15,301) dummy
        read(15,301) dummy
        read(15,*) pf_D, pf_A, pf_S, pf_Tm, solidtype, pf_IC
301     format(a4)
    close(15)

    nx=nxm+1
    ny=nym+1
    nz=nzm+1

    nxr=nxmr+1
    nyr=nymr+1
    nzr=nzmr+1

    ! Prescribe booleans
    if (solidtype.ne.0) IBM = .true.

    ! if(flagstat.ne.0) statcal = .true.
    if(idtv.eq.0) variabletstep = .false.
    ! if(flagbal.ne.0) disscal = .true.
    if(nread.ne.0) readflow = .true.
    if(ireset.ne.0) resetlogstime = .true.
    if(flagmelt.ne.0) melt = .true.
    if(flagMR.ne.0) multires = .true.
    if(flagPF.ne.0) phasefield = .true.
    if(flagsal.ne.0) salinity = .true.

    if (sidewall) call ReadSidewallInput
    if (Non_uniform_BC)      call ReadNon_uniform_BC
    ! if(starea.ne.0) then 
    !   readstats = .true.
    !   if (.not. readflow) write(6,*) 'Warning: Restarting flowfield with statistics read'
    ! endif

!m============================================
!
!   DEFINITIONS FOR THE NATURAL CONVECTION
!
    Rrho = abs(rayt)*pras / (abs(rays)*prat) !CJH inverted to match literature

    if (salinity .and. FFscaleS.eq.1) then    !CJH nondim. velocity with salinity free-fall scale
        ren = dsqrt(abs(rays)/pras)
        byct = Rrho
        bycs = 1.d0
    else                          !CJH nondim. velocity with thermal free-fall scale
        ren = dsqrt(abs(rayt)/prat)
        byct = 1.d0
        bycs = 1.d0/Rrho
    end if

    pect = ren*prat
    pecs = ren*pras
    pi = 2.d0*dasin(1.d0)

    pf_D = 1.2/pect/pf_S/pf_A
    pf_eps = 1.0/nxmr
    pf_Lambda = 2.8e-3/Rrho
    if (ismaster .and. salinity) then
        write(*,*) "Density ratio: ",Rrho
        if (phasefield) then
            write(*,*) "Liquidus slope: ",pf_Lambda
        end if
    end if

    return
end subroutine ReadInputFile

!> Read the sidewall.in input file to specify the boundary conditions in y and z
subroutine ReadSidewallInput
    use param
    implicit none
    integer :: i
    character(len=4) :: dummy
 

    open(unit=15,file='sidewall.in',status='old')
        do i=1,4
          read(15,301) dummy
        end do
        read(15,*) periodic_bc_z_direction
        read(15,301) dummy
        read(15,301) dummy
        read(15,*) bc_vx_y_fix_lo, bc_vx_y_val_lo, bc_vx_y_fix_up,  bc_vx_y_val_up
        read(15,301) dummy
        read(15,*) bc_vx_z_fix_lo, bc_vx_z_val_lo, bc_vx_z_fix_up,  bc_vx_z_val_up
        read(15,301) dummy
        read(15,301) dummy
        read(15,301) dummy
        read(15,*) bc_vy_z_fix_lo, bc_vy_z_val_lo, bc_vy_z_fix_up,  bc_vy_z_val_up
        read(15,301) dummy
        read(15,301) dummy
        read(15,301) dummy
        read(15,*) bc_vz_y_fix_lo, bc_vz_y_val_lo, bc_vz_y_fix_up,  bc_vz_y_val_up
        read(15,301) dummy
        read(15,301) dummy
        read(15,301) dummy
        read(15,*) bc_temp_y_fix_lo, bc_temp_y_val_lo, bc_temp_y_fix_up,  bc_temp_y_val_up
        read(15,301) dummy
        read(15,*) bc_temp_z_fix_lo, bc_temp_z_val_lo, bc_temp_z_fix_up,  bc_temp_z_val_up
        read(15,301) dummy
        read(15,301) dummy
        read(15,301) dummy
        read(15,*) bc_sal_y_fix_lo, bc_sal_y_val_lo, bc_sal_y_fix_up,  bc_sal_y_val_up
        read(15,301) dummy
        read(15,*) bc_sal_z_fix_lo, bc_sal_z_val_lo, bc_sal_z_fix_up,  bc_sal_z_val_up
        read(15,301) dummy
        read(15,301) dummy
        read(15,301) dummy
        read(15,*) bc_phi_y_fix_lo, bc_phi_y_val_lo, bc_phi_y_fix_up,  bc_phi_y_val_up
        read(15,301) dummy
        read(15,*) bc_phi_z_fix_lo, bc_phi_z_val_lo, bc_phi_z_fix_up,  bc_phi_z_val_up
        read(15,301) dummy
        read(15,301) dummy
        read(15,301) dummy
        read(15,*) bc_humid_y_fix_lo, bc_humid_y_val_lo, bc_humid_y_fix_up,  bc_humid_y_val_up
        read(15,301) dummy
        read(15,*) bc_humid_z_fix_lo, bc_humid_z_val_lo, bc_humid_z_fix_up,  bc_humid_z_val_up
    301     format(a4)
    close(15)
    


    if(zlen==0.01 .and. nz == 2)then
        write(*,*) 
        if(bc_vx_z_fix_lo)then
          bc_vx_z_val_lo = 0  
          ErrorSetSideWallBC = .true.
        end if 
  
        if(bc_vx_z_fix_up)then
          bc_vx_z_val_up = 0
          ErrorSetSideWallBC = .true.
        end if 
  
  
        if(bc_temp_z_fix_lo)then
          bc_temp_z_val_lo = 0
          ErrorSetSideWallBC = .true.
        end if 
  
  
        if(bc_temp_z_fix_up)then
          bc_temp_z_val_up = 0
          ErrorSetSideWallBC = .true.
        end if 
  
  
        if(bc_sal_z_fix_lo )then
          bc_sal_z_val_lo = 0
          ErrorSetSideWallBC = .true.
        end if 
  
        if(bc_sal_z_fix_up )then
          bc_sal_z_val_up = 0
          ErrorSetSideWallBC = .true.
        end if
      end if 
  

end subroutine ReadSidewallInput

subroutine ReadNon_uniform_BC
    use param
    implicit none
    integer :: i
    character(len=4) :: dummy

    open(unit=15,file='Non_uniform_BC.in',status='old')
        read(15,301) dummy
        read(15,*) FixValueBCRegion_Length
        read(15,301) dummy
        read(15,*) FixValueBCRegion_Nord_or_Sud
        
    301     format(a4)
    close(15)

end subroutine ReadNon_uniform_BC