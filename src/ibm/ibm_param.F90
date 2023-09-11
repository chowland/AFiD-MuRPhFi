module ibm_param
      implicit none
      integer :: npunx,npuny,npunz,npunte,mpun,mpunr,npuntr
      integer :: n, solidtype
      real :: aldto
      parameter (mpun=50000)
      parameter (mpunr=300000)
      real,allocatable,dimension(:,:) :: plth1, plth2
      logical,allocatable,dimension(:,:,:) :: solidr
      integer,allocatable,dimension(:,:,:) :: ibmaskx,ibmasky,ibmaskz
      integer,allocatable,dimension(:,:,:) :: ibmaskt, ibmaskr
      real,dimension(mpun,3) :: distb
      real,dimension(mpun) :: distx, disty, distz
      real,dimension(mpun) :: distbt
      real,dimension(mpunr) :: distbr
end module ibm_param
