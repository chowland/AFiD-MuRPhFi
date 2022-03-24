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
    implicit none
    character(len=4) :: dummy
    integer flagstat,flagbal,flagmelt
    integer :: flagMR, flagsal, flagPF
    integer :: FFscaleS, pf_IBM
    logical fexist

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
        read(15,*) pf_D, pf_A, pf_S, pf_Tm, pf_IBM, pf_IC
301     format(a4)
    close(15)

    nx=nxm+1
    ny=nym+1
    nz=nzm+1

    nxr=nxmr+1
    nyr=nymr+1
    nzr=nzmr+1

    ! Prescribe booleans
    if (pf_IBM.ne.0) IBM = .true.

    ! if(flagstat.ne.0) statcal = .true.
    if(idtv.eq.0) variabletstep = .false.
    ! if(flagbal.ne.0) disscal = .true.
    if(nread.ne.0) readflow = .true.
    if(ireset.ne.0) resetlogstime = .true.
    if(flagmelt.ne.0) melt = .true.
    if(flagMR.ne.0) multires = .true.
    if(flagPF.ne.0) phasefield = .true.
    if(flagsal.ne.0) salinity = .true.

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
    write(*,*) "Density ratio: ",Rrho
    write(*,*) "Liquidus slope: ",pf_Lambda

    return
end subroutine ReadInputFile