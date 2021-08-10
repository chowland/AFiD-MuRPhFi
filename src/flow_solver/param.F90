
! Declaration of global variables
!***********************************************************
module param
    implicit none
    !==========================================================
    !       read from input file bou.in
    !==========================================================
    integer   :: nx, ny, nz
    integer   :: nxr, nyr, nzr, istr3r  !CS mgrd
    integer   :: nsst, nread, ntst, ireset
    real      :: walltimemax,tout,tmax
    real      :: alx3,str3
    integer   :: istr3
    real      :: ylen,zlen
    real      :: rayt,prat,rays,pras,dt,resid
    integer   :: inslws,inslwn
    integer   :: TfixS,TfixN,SfixS,SfixN        !CJH option for fixed T/S BCs
    integer   :: gAxis      !CJH option to choose gravity axis
    real      :: xminusU,xplusU,dPdz,dPdy
    real      :: pf_A, pf_eps, pf_C, pf_S, pf_Tm, pf_Lambda
    logical   :: IBM
    ! integer   :: starea,tsta
    real      :: dtmin,dtmax,limitCFL
    integer   :: nson,idtv
    real      :: tframe, save_3D
    integer   :: active_T, active_S, pf_IC !CJH Option for passive scalars
    !=================================================
    !       end of input file
    !=================================================
    real :: time
    !******* Grid parameters**************************
    real :: dx,dy,dz,dxq,dyq,dzq
    real :: dxr,dyr,dzr,dxqr,dyqr,dzqr      !CS mgrd
    !
    real, allocatable, dimension(:) :: xc,xm
    real, allocatable, dimension(:) :: yc,ym
    real, allocatable, dimension(:) :: zc,zm
    real, allocatable, dimension(:) :: g3rc,g3rm
    real, allocatable, dimension(:) :: d3xc, d3xm    !CJH modified derivatives
    real, allocatable, dimension(:) :: xcr,xmr       !CS mgrd
    real, allocatable, dimension(:) :: ycr,ymr       !CS mgrd
    real, allocatable, dimension(:) :: zcr,zmr       !CS mgrd
    real, allocatable, dimension(:) :: g3rcr,g3rmr   !CS mgrd
    real, allocatable, dimension(:) :: d3xcr, d3xmr  !CJH modified derivatives
    !====================================================
    !******* QUANTITIES FOR DERIVATIVES******************
    real, allocatable, dimension(:) :: udx3c,udx3m
    real, allocatable, dimension(:) :: udx3cr,udx3mr   !CS mgrd
    !==========================================================
    !******* Grid indices**************************************
    integer, allocatable, dimension(:) :: kmc,kpc,kmv,kpv
    integer, allocatable, dimension(:) :: kmcr,kpcr,kmvr,kpvr  !CS mgrd
    !===========================================================
    !******* Metric coefficients *******************************
    real, allocatable, dimension(:) :: ap3ck,ac3ck,am3ck
    real, allocatable, dimension(:) :: ap3ckr,ac3ckr,am3ckr      !CS mgrd
    real, allocatable, dimension(:) :: ap3sk,ac3sk,am3sk
    real, allocatable, dimension(:) :: ap3ssk,ac3ssk,am3ssk
    real, allocatable, dimension(:) :: ap3sskr,ac3sskr,am3sskr   !CS mgrd
    real, allocatable, dimension(:) :: ap3spkr,ac3spkr,am3spkr   !CJH phase-field
    !============================================================
    !******* Variables for FFTW and Poisson solver****************
    real, allocatable, dimension(:) :: ak2,ap
    real, allocatable, dimension(:) :: ak1,ao
    real, allocatable, dimension(:) :: amphk,acphk,apphk

    !===========================================================
    !******* Other variables ***********************************
    integer  :: nxm, nym, nzm
    integer  :: nxmr, nymr, nzmr   !CS mgrd
    real :: ren, pect, pecs
    real :: Rrho, byct, bycs
    real :: pi
    real :: al,ga,ro
    real :: beta
    real :: qqmax,qqtot
    real :: re
    real :: tempmax,tempmin,tempm
    integer :: ntime
    integer, parameter:: ndv=3
    real, dimension(1:ndv) :: vmax
    real, dimension(1:3) :: gam,rom,alm
    real, allocatable, dimension(:,:,:) :: tempbp,temptp !CJH make BCs 3D arrays so we can use update_halo
    real, allocatable, dimension(:,:,:) :: salbp,saltp
    integer, dimension(5) :: spec_idx

    logical :: dumpslabs=.false.
    ! logical :: statcal=.false.
    ! logical :: disscal=.false.
    logical :: readflow=.false.
    logical :: readstats=.false.
    logical :: ismaster=.false.
    logical :: resetlogstime=.false.
    logical :: variabletstep=.true.
    logical :: melt=.false.
    logical :: multires=.false.
    logical :: phasefield=.false.
    logical :: salinity=.false.
    logical :: writespec=.false.

    integer :: lvlhalo=2

end module param

!************* End of param module******************************
!===============================================================
!******* 2D arrays, dynamically allocated by each process*******
module local_arrays
    use param
    implicit none
    real,allocatable,dimension(:,:,:) :: vx,vy,vz
    real,allocatable,dimension(:,:,:) :: pr,temp,rhs
    real,allocatable,dimension(:,:,:) :: rux,ruy,ruz,rutemp
    real,allocatable,dimension(:,:,:) :: dph,qcap,dq,hro,dphhalo
    real,allocatable,dimension(:,:,:) :: qtens
    ! real,allocatable,dimension(:,:,:) :: sal,rhsr,rusal,hsal
end module local_arrays

module mgrd_arrays
    use param
    implicit none
    integer,allocatable,dimension(:) :: irangs,jrangs,krangs
    integer,allocatable,dimension(:) :: irangc,jrangc,krangc
    real,allocatable,dimension(:,:,:) :: sal,rhsr,rusal,hsal
    real,allocatable,dimension(:,:) :: cxvx, cxvy, cxvz, cxrs, cxsalc
    real,allocatable,dimension(:,:) :: cyvx, cyvy, cyvz, cyrs, cysalc
    real,allocatable,dimension(:,:) :: czvx, czvy, czvz, czrs, czsalc
    real,allocatable,dimension(:,:,:) :: vxr,vyr,vzr !CS mgrd
    real,allocatable,dimension(:,:,:) :: tpdv,tpdvr  !CS mgrd
    real,allocatable,dimension(:,:,:) :: salc
    real,allocatable,dimension(:,:,:) :: tempr,Tplaner
    real,allocatable,dimension(:,:,:) :: phi,phic,ruphi,hphi
end module mgrd_arrays
!===============================================================
module stat_arrays
    implicit none
    real,allocatable, dimension(:) :: vz_me,vz_rms
    real,allocatable, dimension(:) :: vy_me,vx_me,vy_rms,vx_rms
    real,allocatable, dimension(:) :: temp_me,temp_rms
    real,allocatable, dimension(:) :: disste,dissth,tempvx_me
    real,allocatable, dimension(:) :: vxvy_corr
    real,allocatable, dimension(:) :: vz_me_buf,vy_me_buf,vx_me_buf
    real,allocatable, dimension(:) :: vz_msq_buf,vy_msq_buf,vx_msq_buf
    real :: vx_global, vy_global, vz_global
    integer :: nstatsamples
end module stat_arrays
!=====================================================
module input_grids !CJH Grids for input flow field
    implicit none
    integer :: nxo, nyo, nzo
    integer :: nxmo, nymo, nzmo
    integer :: nxro, nyro, nzro
    integer :: nxmro, nymro, nzmro
    integer :: istr3o, istr3ro
    integer, dimension(3) :: xstarto, xendo, xsizeo
    real :: str3o
    real, allocatable, dimension(:) :: xco, xmo
    real, allocatable, dimension(:) :: xcro, xmro
    real, allocatable, dimension(:) :: yco, ymo
    real, allocatable, dimension(:) :: ycro, ymro
    real, allocatable, dimension(:) :: zco, zmo
    real, allocatable, dimension(:) :: zcro, zmro
    integer,allocatable,dimension(:) :: irangsr,jrangsr,krangsr
end module input_grids
!=====================================================
module stat3_param
    implicit none
    integer :: kslab(1:9)
    real    :: xslab(1:9)
end module stat3_param
!=====================================================
module mpih
    implicit none
    include 'mpif.h'
    integer :: ierr
    integer, parameter :: master=0
    integer :: MDP = MPI_DOUBLE_PRECISION
end module mpih
!====================================================
module fftw_params
    !        use param, only: m2m,m2mh,m1m
    use iso_c_binding

    type, bind(C) :: fftw_iodim
        integer(C_INT) n, is, os
    end type fftw_iodim

    interface
        type(C_PTR) function fftw_plan_guru_dft(rank,dims, &
        howmany_rank,howmany_dims,in,out,sign,flags) &
        bind(C, name='fftw_plan_guru_dft')
        import
        integer(C_INT), value :: rank
        type(fftw_iodim), dimension(*), intent(in) :: dims
        integer(C_INT), value :: howmany_rank
        type(fftw_iodim), dimension(*), intent(in) :: howmany_dims
        complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: in
        complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: out
        integer(C_INT), value :: sign
        integer(C_INT), value :: flags
        end function fftw_plan_guru_dft

        type(C_PTR) function fftw_plan_guru_dft_r2c(rank,dims, &
            howmany_rank,howmany_dims,in,out,flags) &
            bind(C, name='fftw_plan_guru_dft_r2c')
            import
            integer(C_INT), value :: rank
            type(fftw_iodim), dimension(*), intent(in) :: dims
            integer(C_INT), value :: howmany_rank
            type(fftw_iodim), dimension(*), intent(in) :: howmany_dims
            real(C_DOUBLE), dimension(*), intent(out) :: in
            complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: out
            integer(C_INT), value :: flags
        end function fftw_plan_guru_dft_r2c

        type(C_PTR) function fftw_plan_guru_dft_c2r(rank,dims, &
            howmany_rank,howmany_dims,in,out,flags)  &
            bind(C, name='fftw_plan_guru_dft_c2r')
            import
            integer(C_INT), value :: rank
            type(fftw_iodim), dimension(*), intent(in) :: dims
            integer(C_INT), value :: howmany_rank
            type(fftw_iodim), dimension(*), intent(in) :: howmany_dims
            complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: in
            real(C_DOUBLE), dimension(*), intent(out) :: out
            integer(C_INT), value :: flags
        end function fftw_plan_guru_dft_c2r


    end interface

    integer FFTW_PATIENT, FFTW_FORWARD, FFTW_BACKWARD,FFTW_ESTIMATE
    parameter (FFTW_PATIENT=32)
    parameter (FFTW_ESTIMATE=64)
    parameter (FFTW_FORWARD=-1)
    parameter (FFTW_BACKWARD=1)
    type(C_PTR) :: fwd_guruplan_y,bwd_guruplan_y
    type(C_PTR) :: fwd_guruplan_z,bwd_guruplan_z
    logical :: planned=.false.

    real,allocatable,dimension(:,:,:) :: ry1,rz1
    complex,allocatable,dimension(:,:,:) :: cy1,cz1,dphc,fouvar1,fouvar2

end module fftw_params