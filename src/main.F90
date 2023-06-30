program AFiD
    use mpih
    use param
    use local_arrays, only: vx,vy,vz,temp,pr
    ! use mgrd_arrays, only: phic,tempr,phi!,vxr,vyr,vzr,salc,sal
    use AuxiliaryRoutines
    use hdf5
    use decomp_2d
    use decomp_2d_fft
    use afid_pressure
    use afid_moisture
    use afid_salinity
    use afid_phasefield
    use h5_tools, only: InitSliceCommunicators
    ! use stat_arrays, only: nstatsamples,vx_global,vy_global,vz_global

!$    use omp_lib
    implicit none
    integer :: errorcode!, nthreads, i, j, k
    real    :: instCFL,CFLmr,dmax,dmaxr
    real    :: ti(2), tin(3), minwtdt
    real :: ts!, varptb,chksum
    real :: td(2)   !< debugging time measure
    integer :: prow=0,pcol=0
    integer :: lfactor,lfactor2
    character(100) :: arg
    logical :: nanexist, write_mean_planes=.true.
    ! real,allocatable,dimension(:,:) :: dummy,dscan,dbot
    ! integer :: comm,ierror,row_id,row_coords(2),ic,jc,kc

!*******************************************************
!******* Read input file bou.in by all processes********
!*******************************************************
!

    ! Set `moist` to .true. if humid.in exists
    inquire(file="humid.in", exist=moist)

    call ReadInputFile
    if (nzm==1 .or. nym==1) write_mean_planes = .false.

    if (command_argument_count().eq.2) then
        call get_command_argument(1,arg)
        read(arg,'(i10)') prow
        call get_command_argument(2,arg)
        read(arg,'(i10)') pcol
    endif

    call decomp_2d_init(nxm ,nym ,nzm ,&
                        nxmr,nymr,nzmr,&
                        prow,pcol,&
                        (/ .false.,.true.,.true. /))

    ts=MPI_WTIME()
    tin(1) = MPI_WTIME()

    call MpiBarrier

    call HdfStart

    if (nrank.eq.master) ismaster = .true.

    if (ismaster) write(6,*) 'MPI tasks=', nproc

!$    if (ismaster) then
!$OMP PARALLEL
!$OMP MASTER
!$        nthreads = omp_get_num_threads()
!$OMP END MASTER
!$OMP END PARALLEL
!$        write(6,*) 'OMP threads=', nthreads
!$    end if

    if(ismaster) then
        call system('mkdir outputdir outputdir/flowmov outputdir/fields')
!m==========================================
        call ResetLogs
!m====================================================
        write(6,112)ylen/alx3,zlen/alx3
    112 format(//,20x,'Double Diffusive Convection ',//,10x, &
        '3D Cell with aspect-ratio:  D1/H = ',f5.2,' D2/H = ',f5.2)
        write(6,142)
    142 format(//,8x,'Periodic lateral wall boundary condition')
        write(6,202) rayt,prat,rays,pras
    202 format(/,5x,'Parameters: ',' RaT=',es10.3,' PrT= ',es10.3,' RaS=',es10.3,' PrS= ',es10.3)
        if(variabletstep) then
            write(6,204) limitCFL
        204 format(/,5x,'Variable dt and fixed cfl= ', es11.4,/ )
        else
            write(6,205) dtmax,limitCFL
        205 format(/,5x,'Fixed dt= ',es11.4,' and maximum cfl=', es11.4,/ )
        endif
    endif

    call InitTimeMarchScheme

    call InitVariables
    call InitPressureVars
    if (multires) call InitMgrdVariables  !CS mgrd
    if (salinity) call InitSalVariables
    if (phasefield) call InitPFVariables
    if (moist) call InitMoistVariables

    call InitSliceCommunicators

    call CreateGrid
    if (multires) call CreateMgrdGrid     !CS mgrd

    call WriteGridInfo

    inquire(file=trim("spectra.in"),exist=specwrite)
    if (specwrite) call InitPowerSpec

!m===================================
!m===================================
!m===================================

    if(ismaster) then
        write(6,754)nx,ny,nz
    754 format(/,5x,'grid resolution: ',' nx = ',i5,' ny = ',i5,' nz = ',i5)
        write(6,756)nxr,nyr,nzr
    756 format(5x,'grid resolution: ',' nxr= ',i5,' nyr= ',i5,' nzr= ',i5)
        write(6,755) 1.d0/dx,1.d0/dy,1.d0/dz,dt,ntst
    755 format(/,2x,' dx=',es10.3,' dy=',es10.3,' dz=',es10.3, &
                  ' dt=',es10.3,' ntst=',i7,/)
        if (multires) then
            write(6,757) 1.d0/dxr,1.d0/dyr,1.d0/dzr
        757 format(/,2x,' dxr=',es10.3,' dyr=',es10.3,' dzr=',es10.3,/)
        end if
    endif

!m===================================
!m===================================
    if (ismaster .and. multires) then
        if ((modulo(nym,prow) /= 0) .or. (modulo(nzm,pcol) /= 0) .or. &
            (modulo(nymr,prow)/= 0) .or. (modulo(nzmr,pcol)/= 0)) then
            write(*,*) "********** WARNING **********"
            write(*,*) "Grid size not a perfect factor of the pencil decomposition"
            write(*,*) "This would cause severe issues with the multi-grid interpolation"
            write(*,*) "Terminating simulation..."
            tin(2) = MPI_WTIME()
            errorcode = 666
            call QuitRoutine(tin, .false., errorcode)
        end if
    end if

    time=0.d0
    ! if(statcal) nstatsamples = 0

    call InitPressureSolver
    call SetTempBCs
    if (salinity) call SetSalBCs
    if (moist) call SetHumidityBCs

    if(readflow) then

        if(ismaster) write(6,*) 'Reading initial condition from file'

        call ReadFlowInterp(prow,pcol)

    else

        if(ismaster) write(6,*) 'Creating initial condition'

        ntime=0
        time=0.d0
        instCFL=0.d0

        call CreateInitialConditions
        if (salinity) call CreateInitialSalinity
        if (phasefield) call CreateInitialPhase
        if (moist) call CreateInitialHumidity

    endif

!CS   Create multigrid stencil for interpolation
    if (multires) call CreateMgrdStencil
    if (phasefield .and. IBM) call CreatePFStencil

    if (phasefield) call update_halo(phi,lvlhalo)

    if (IBM) then
        call topogr
        ! if (phasefield) call UpdateIBMLocation
    end if

!EP   Update all relevant halos
    call update_halo(vx,lvlhalo)
    call update_halo(vy,lvlhalo)
    call update_halo(vz,lvlhalo)
    call update_halo(temp,lvlhalo)
    if (salinity) call update_halo(sal,lvlhalo)
    call update_halo(pr,lvlhalo)
    if (moist) call update_halo(humid,lvlhalo)


!CS   Interpolate initial values
    if (salinity) then
        call InterpVelMgrd
        call InterpSalMultigrid
    end if
    if (phasefield) then
        call InterpTempMultigrid
        call InterpPhiMultigrid
    end if

!EP   Update all relevant halos
    if (salinity) then
        call update_halo(vxr,lvlhalo)
        call update_halo(vyr,lvlhalo)
        call update_halo(vzr,lvlhalo)
        call update_halo(salc,lvlhalo)
    end if
    if (phasefield) then
        call update_halo(tempr,lvlhalo)
        call update_halo(phic,lvlhalo)
    end if

    call CalcMeanProfiles
    if (specwrite) call WritePowerSpec
    if(ismaster)  write(6,*) 'Write plane slices'
    call Mkmov_xcut
    call Mkmov_ycut
    call Mkmov_zcut
    if (write_mean_planes) then
        call mean_yplane
        call mean_zplane
    end if
    
    if (ismaster) write(*,*) "Writing 3D fields"
    call MpiBarrier
    td(1) = MPI_WTIME()
    call WriteFlowField(.false.)
    call MpiBarrier
    td(2) = MPI_WTIME()
    if (ismaster) write(*,*) "Flow field writing took: ",td(2)-td(1)

!EP   Check divergence. Should be reduced to machine precision after the first
!phcalc. Here it can still be high.

    call CheckDivergence(dmax,dmaxr)

    if(ismaster) write(6,*)' Initial maximum divergence: ',dmax,dmaxr

!!EP   Write some values
!      if(variabletstep) then
!       if(ismaster) write(6,*)ntime,time,dt,dmax,tempm,tempmax,tempmin
!      else
!       if(ismaster) write(6,*)ntime,time,dt,dmax,tempm,tempmax,tempmin !RO Fix??
!      end if
!

    if(ismaster) then
       tin(2) = MPI_WTIME()
       write(6,'(a,f6.2,a)') 'Initialization Time = ', tin(2) -tin(1), ' sec.'
    endif

!  ********* starts the time dependent calculation ***
    errorcode = 0 !EP set errocode to 0 (OK)
    minwtdt = huge(0.0d0) !EP initialize minimum time step walltime

    ! Check input for efficient FFT
    ! factorize input FFT directions. The largest factor should
    ! be relatively small to have efficient FFT's
    lfactor=2 ! initialize
    call Factorize(nym,lfactor2) ! check nym
    lfactor=max(lfactor,lfactor2)
    call Factorize(nzm,lfactor2)
    lfactor=max(lfactor,lfactor2)
    ! if largest factor larger than 7 quit the simulation
    ! if (lfactor>7) errorcode=444

    do ntime=0,ntst
        ti(1) = MPI_WTIME()

!EP   Determine timestep size
        call CalcMaxCFL(instCFL,CFLmr)

        if(variabletstep) then
            if(ntime.gt.1) then
                if(CFLmr.lt.1.0d-8) then !EP prevent fp-overflow
                    dt=dtmax
                else
                    dt=limitCFL/CFLmr
                endif
                if(dt.gt.dtmax) dt=dtmax
            else
                dt=dt !CJH First time step: use dt defined in bou.in
            endif
            if(dt.lt.dtmin) errorcode = 166
        else
!RO    fixed time-step
            CFLmr=CFLmr*dt
            if(CFLmr.gt.limitCFL) errorcode = 165
            ! if (ismaster) write(*,*) "CFL value  ",CFLmr
        endif

        call TimeMarcher

        time=time+dt

        if(mod(time,tout).lt.dt) then
            if(ismaster) then
                write(6,*) ' -------------------------------------------------- '
                write(6,'(a,ES11.4,a,i9,a,ES11.4)') '  T = ',time,' NTIME = ',ntime,' DT = ',dt
            endif
            call CalcMeanProfiles
            if (specwrite) then
                if (ismaster) write(*,*) "Writing power spectra"
                call WritePowerSpec
                if (ismaster) write(*,*) "Done writing power spectra"
            end if
            if(ismaster) then
                open(96,file='outputdir/cfl.out',status='unknown',position='append',access='sequential')
                write(96,769) ntime,time,dt,instCFL*dt!,vx_global,vy_global,vz_global
                close(96)
            endif
        769 format(1x,i12,3(1x,ES20.8))
        endif

        if((mod(time,tframe).lt.dt) .and. (floor(time/tframe).ne.0)) then
            if(ismaster)  write(6,*) 'Write slice ycut and zcut'
            !call CalcWriteQ
            call Mkmov_xcut
            call Mkmov_ycut
            call Mkmov_zcut
            if (write_mean_planes) then
                call mean_yplane
                call mean_zplane
            end if
        endif

        if(ntime.eq.1.or.mod(time,tout).lt.dt) then
            !call GlobalQuantities
            !if(vmax(1).gt.limitVel.and.vmax(2).gt.limitVel) errorcode = 266

            call CalcMaxCFL(instCFL,CFLmr)
            call CheckDivergence(dmax,dmaxr)
            !  call CalcPlateNu
            !call CalcPlateCf

            if(.not.variabletstep) instCFL=instCFL*dt

            if(abs(dmax).gt.resid) errorcode = 169

        endif

        if(time.gt.tmax) errorcode = 333

        ti(2) = MPI_WTIME()
        minwtdt = min(minwtdt,ti(2) - ti(1))
        if(mod(time,tout).lt.dt) then
            if(ismaster) then
                write(6,*) ' Maximum divergence = ', dmax, dmaxr
                !write(6,*)'ntime - time - vmax(1) - vmax(2) - vmax(3)  -&
                !           tempm - tempmax - tempmin'
                !write(6,*)ntime,time,vmax(1),vmax(2),vmax(3),tempm,tempmax,tempmin
                write(6,'(a,f8.3,a)') '  Minimum Iteration Time = ', minwtdt,' sec.'
            endif
            minwtdt = huge(0.0d0)
        endif

        if ((mod(time,save_3D).lt.dt) .and. (floor(time/tframe).ne.0)) then
            if(ismaster) write(6,*) '*** Writing 3D fields ***'
            call WriteFlowField(.false.)
            call MpiBarrier
            if(ismaster) write(6,*) '********* Done **********'
        end if

        if( (ti(2) - tin(1)) .gt. walltimemax) errorcode = 334

        if( ntime .eq. ntst ) errorcode = 555

        call MpiBcastInt(errorcode)

!EP   Conditional exits
        if(errorcode.ne.0) then

!EP    dt too small
            if(errorcode.eq.166) call QuitRoutine(tin,.false.,errorcode)

!EP   cfl too high
            if(errorcode.eq.165) call QuitRoutine(tin,.false.,errorcode)

!EP   velocities diverged
            if(errorcode.eq.266) call QuitRoutine(tin,.false.,errorcode)

!EP   mass not conserved
            if(errorcode.eq.169) call QuitRoutine(tin,.false.,errorcode)

!EP   Physical time exceeded tmax, no error; normal quit
            if(errorcode.eq.333) call QuitRoutine(tin,.true.,errorcode)

!EP   walltime exceeded walltimemax, no error; normal quit
            if(errorcode.eq.334) call QuitRoutine(tin,.true.,errorcode)

!RS   FFT input not correct
            if(errorcode.eq.444) call QuitRoutine(tin,.false.,errorcode)

!RS   maximum number of timesteps reached, no error; normal quit
            if(errorcode.eq.555) call QuitRoutine(tin,.true.,errorcode)

            errorcode = 100 !EP already finalized

            exit

        endif

    enddo !EP main loop

    call QuitRoutine(tin,.true.,errorcode)

end