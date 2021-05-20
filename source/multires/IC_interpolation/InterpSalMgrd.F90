      subroutine InterpSalMgrd

      use param
      use mgrd_arrays, only: sal,salc,cxsalc,cysalc,czsalc,irangs,jrangs,krangs
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
       icr = krangs(ic)-1

       do jc=xstart(2),xend(2)
         jcr = jrangs(jc)-1

         do kc=1,nxm
          kcr = irangs(kc)-1

          if (kcr==1) then
            if (SfixS==1) then    !CJH apply lower fixed value BC
              qv3(1,:,:) = 2.d0*salbp(1,jcr-1:jcr+2,icr-1:icr+2) &
                            - sal(kcr,jcr-1:jcr+2,icr-1:icr+2)
            else    !CJH apply no flux BC
              qv3(1,:,:) = sal(kcr,jcr-1:jcr+2,icr-1:icr+2)
            end if
            qv3(2:,:,:) = sal(kcr:kcr+2,jcr-1:jcr+2,icr-1:icr+2)
          else if (kcr==nxmr) then
            qv3(1:3,:,:) = sal(kcr-1:kcr+1,jcr-1:jcr+2,icr-1:icr+2)
            if (SfixN==1) then    !CJH apply upper fixed value BC
              qv3(4,:,:) = 2.d0*saltp(1,jcr-1:jcr+2,icr-1:icr+2) &
                            - sal(kcr,jcr-1:jcr+2,icr-1:icr+2)
            else    !CJH apply no flux BC
              qv3(4,:,:) = sal(kcr,jcr-1:jcr+2,icr-1:icr+2)
            end if
          else
            qv3 = sal(kcr-1:kcr+2,jcr-1:jcr+2,icr-1:icr+2)
          end if
          qv2(:,:) = qv3(:,:,1)*czsalc(1,ic)+qv3(:,:,2)*czsalc(2,ic)&
                    +qv3(:,:,3)*czsalc(3,ic)+qv3(:,:,4)*czsalc(4,ic)
          qv1(:)   = qv2(:,1)*cysalc(1,jc)+qv2(:,2)*cysalc(2,jc) &
                    +qv2(:,3)*cysalc(3,jc)+qv2(:,4)*cysalc(4,jc)
          salc(kc,jc,ic) = sum(qv1(1:4)*cxsalc(1:4,kc))

         enddo
       enddo
      enddo

      !CJH Not appropriate for staggered grid
      ! do ic=xstart(3),xend(3)
      !  do jc=xstart(2),xend(2)
      !     salc(1,jc,ic) = 0.d0
      !  enddo
      ! enddo

      return
      end subroutine InterpSalMgrd
