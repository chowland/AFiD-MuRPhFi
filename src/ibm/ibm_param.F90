module ibm_param
      implicit none
      integer :: npunx,npuny,npunz,npunte,mpun,npuntr
      integer :: n, solidtype
      real :: aldto
      parameter (mpun=50000)
      real,allocatable,dimension(:,:,:) :: forclo, forclor
      real,allocatable,dimension(:,:) :: plth1, plth2
      logical,allocatable,dimension(:,:,:) :: solidr
      integer,dimension(3,mpun,3) :: indgeo, indgeoe
      integer,dimension(mpun,3) :: indgeot, indgeoet
      integer,dimension(mpun,3) :: indgeor, indgeoer
      real,dimension(3,mpun) :: distb
      real,dimension(mpun) :: distbt,distbr
      real,dimension(mpun) :: temb,salfix
      real,dimension(mpun) :: q1bo,q2bo,q3bo,densb,salb
      
end module ibm_param

! SL =======================================================
! module param_particle
!       implicit none
!       integer :: Nparticle,NL,Ne,npp,nll,nee,partshape,thermal_couple,iangle
!       integer :: imonitor,MonitorNumber,iJointPDF,HSectPosi
!       integer, allocatable, dimension(:,:) :: MonitorPosition
!       integer, allocatable, dimension(:,:,:) :: allsdind
!       real :: radius,diameter,kparticle,friction,friction2,TOUTHSect,HSectCoor,rhocpparticle
!       real :: dVL,rx,ry,Am,delta,angular_velocity
!       real, allocatable, dimension(:)  :: vec_angle,part_surf_vel_x,part_surf_vel_y
!       real, allocatable, dimension(:)  :: HSectTemp,HSectVx,HSectVy
!       real, allocatable, dimension(:,:)    :: ParticleCenter
!       real, allocatable, dimension(:,:,:)  :: phaindx,phaindy,phaindz,kcp,insidex,insidey,insidez,rhocpcp
!       real, allocatable, dimension(:,:,:)  :: vxcp,vycp,vzcp,vxcp_mean,vycp_mean
!       real, allocatable, dimension(:,:,:)  :: vxc_ibm,vyc_ibm,vzc_ibm
!       real, allocatable, dimension(:,:)    :: lagnodx,lagnody,ffx,ffy
!       real, allocatable, dimension(:,:,:)  :: mls_transfer
!       real, allocatable, dimension(:,:)    :: mls_h,mls_c
!       !        real, allocatable, dimension(:,:,:)  :: delta_vx_index,delta_vy_index,delta_vz_index
! end module param_particle
      
      
! module param_tracer
!       implicit none
!       real, allocatable, dimension(:,:,:) :: vxc,vyc,vzc,tempc
!       real, allocatable, dimension(:,:,:) :: str_11,str_22,str_33
!       real, allocatable, dimension(:,:,:) :: str_12,str_13,str_23
!       real, allocatable, dimension(:,:,:) :: vor_y,vor_z,vor_x
      
!       real, allocatable, dimension(:,:) :: xp,xpo
!       real, allocatable, dimension(:) :: vxp,vyp,vzp,Tp 
!       real, allocatable, dimension(:) :: vxop,vyop,vzop
      
!       real, allocatable, dimension(:) :: str11p,str22p,str33p
!       real, allocatable, dimension(:) :: str12p,str13p,str23p
!       real, allocatable, dimension(:) :: voryp,vorzp,vorxp
      
!       integer iftracer,Ntracer,ONtracer,tracer_write,tracer_read,if_tracer_pair
!       real tracer_init_sep
!       real x_i,y_i,z_i,x_f,y_f,z_f,z_c
!       real timeONtracer,TOUT_tracer
!       integer dimen
! end module param_tracer
! ==========================================================
