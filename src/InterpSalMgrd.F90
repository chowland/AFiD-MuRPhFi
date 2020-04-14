      subroutine InterpSalMgrd

      use param
      use local_arrays, only: sal
      use mgrd_arrays, only: salc,cxvxc,cyvxc,czvxc,irangc,jrangs,krangs
      use mpih
      use decomp_2d
      use AuxiliaryRoutines
      implicit none

      integer  :: ic,jc,kc, icr,jcr,kcr

      real,dimension(4,4,4) :: qv3 
      real,dimension(4,4) :: qv2
      real,dimension(4) :: qv1

      ! Interpolate to coarse grid here. A better option is to apply
      ! a box filter.
      salc(:,:,:) = 0.d0

      do ic=xstart(3),xend(3)
       icr = irangc(ic)-1

       do jc=xstart(2),xend(2)
         jcr = jrangs(jc)-1

         do kc=1,nxm
          kcr = krangs(kc)-1

          qv3      = sal(kcr-1:kcr+2,jcr-1:jcr+2,icr-1:icr+2)
          qv2(:,:) = qv3(:,:,1)*czvxc(1,ic)+qv3(:,:,2)*czvxc(2,ic)&
                    +qv3(:,:,3)*czvxc(3,ic)+qv3(:,:,4)*czvxc(4,ic)
          qv1(:)   = qv2(:,1)*cyvxc(1,jc)+qv2(:,2)*cyvxc(2,jc) &
                    +qv2(:,3)*cyvxc(3,jc)+qv2(:,4)*cyvxc(4,jc)
          salc(kc,jc,ic) = sum(qv1(1:4)*cxvxc(1:4,kc))

         enddo
       enddo
      enddo

      return
      end subroutine InterpSalMgrd
