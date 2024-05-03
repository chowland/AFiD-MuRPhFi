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
    use afid_sides, only: ReadSidewallInput
    implicit none
    character(len=4) :: dummy
    integer :: flagmelt, io
    integer :: flagMR, flagsal, flagPF
    integer :: FFscaleS

    open(newunit=io,file='bou.in',status='old')
        read(io,301) dummy
        read(io,301) dummy
        read(io,*) nxm, nym, nzm, nsst
        read(io,301) dummy
        read(io,301) dummy
        read(io,301) dummy
        read(io,*) flagMR, nxmr, nymr, nzmr
        read(io,301) dummy
        read(io,301) dummy
        read(io,301) dummy
        read(io,*) flagsal, flagPF
        read(io,301) dummy
        read(io,301) dummy
        read(io,301) dummy
        read(io,*) nread, ireset
        read(io,301) dummy
        read(io,301) dummy
        read(io,301) dummy
        read(io,*) ntst, walltimemax, tmax
        read(io,301) dummy
        read(io,301) dummy
        read(io,301) dummy
        read(io,*) tout, tframe, save_3D
        read(io,301) dummy
        read(io,301) dummy
        read(io,301) dummy
        read(io,*) alx3, ylen, zlen
        read(io,301) dummy
        read(io,301) dummy
        read(io,301) dummy
        read(io,*) istr3, str3, istr3r
        read(io,301) dummy
        read(io,301) dummy
        read(io,301) dummy
        read(io,*) rayt, prat, rays, pras, FFscaleS
        read(io,301) dummy
        read(io,301) dummy
        read(io,301) dummy
        read(io,*) idtv, dt, resid, limitCFL, dtmin, dtmax
        read(io,301) dummy
        read(io,301) dummy
        read(io,301) dummy
        read(io,*) inslws, inslwn, TfixS, TfixN, SfixS, SfixN
        read(io,301) dummy
        read(io,301) dummy
        read(io,301) dummy
        read(io,*) active_T, active_S, gAxis
        read(io,301) dummy
        read(io,301) dummy
        read(io,301) dummy
        read(io,*) xplusU, xminusU, dPdy, dPdz, flagmelt
        read(io,301) dummy
        read(io,301) dummy
        read(io,301) dummy
        read(io,*) pf_D, pf_A, pf_S, pf_Tm, solidtype, pf_IC
301     format(a4)
    close(io)

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
