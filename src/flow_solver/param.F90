
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
    real      :: rayt,prat,dt,resid
    integer   :: inslws,inslwn
    integer   :: TfixS,TfixN    !CJH option for fixed T/S BCs
    integer   :: gAxis      !CJH option to choose gravity axis
    real      :: xminusU,xplusU,dPdz,dPdy
    logical   :: IBM, moist
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
    real, allocatable, dimension(:) :: ap3sk,ac3sk,am3sk
    real, allocatable, dimension(:) :: ap3ssk,ac3ssk,am3ssk
    !============================================================

    !===========================================================
    !******* Other variables ***********************************
    integer  :: nxm, nym, nzm
    integer  :: nxmr, nymr, nzmr   !CS mgrd
    real :: ren, pect!, pecs
    real :: Rrho, byct!, bycs
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
    logical :: specwrite=.false.

    integer :: lvlhalo=2

    logical :: sidewall = .false.     !! Flag to determine whether to impose sidewalls in y and z (using a DCT)

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
end module local_arrays

module mgrd_arrays
    use param
    implicit none
    integer,allocatable,dimension(:) :: irangs,jrangs,krangs
    integer,allocatable,dimension(:) :: irangc,jrangc,krangc
    integer,allocatable,dimension(:) :: irangb,jrangb,krangb
    integer,allocatable,dimension(:) :: irangr,jrangr,krangr
    integer,allocatable,dimension(:) :: yc_to_ymr, zc_to_zmr
    real,allocatable,dimension(:,:,:) :: rhsr
    real,allocatable,dimension(:,:) :: cxvx, cxvy, cxvz, cxrs, cxsalc, cxphic
    real,allocatable,dimension(:,:) :: cyvx, cyvy, cyvz, cyrs, cysalc, cyphic
    real,allocatable,dimension(:,:) :: czvx, czvy, czvz, czrs, czsalc, czphic
    real,allocatable,dimension(:,:) :: cych, czch
    real,allocatable,dimension(:,:,:) :: tpdv,tpdvr  !CS mgrd
    real,allocatable,dimension(:,:,:) :: Tplaner
    real,allocatable,dimension(:,:) :: solid_height, height_vx, height_vy, height_vz
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
    use mpi
    implicit none
    integer :: ierr
    integer, parameter :: master=0
    integer :: MDP = MPI_DOUBLE_PRECISION
    integer :: comm_yz  !> MPI communicator across y and z (2D decomposition)
    integer :: comm_xy  !> MPI communicator across x and y (1D decomposition)
    integer :: comm_xz  !> MPI communicator across x and z (1D decomposition)
end module mpih
!====================================================