!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: CalcWriteQ.F90                                 !
!    CONTAINS: subroutine CalcWriteQ                      !
!                                                         ! 
!    PURPOSE: Compute and write the 3D q-criteria field   !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
      subroutine CalcWriteQ
      use param
      use local_arrays, only: vz,vy,vx,qtens
      use decomp_2d
      use AuxiliaryRoutines
      implicit none

!      real, allocatable, dimension(:,:,:) :: qtens
      real :: dvxx1,dvxx2,dvxx3
      real :: dvyx1,dvyx2,dvyx3
      real :: dvzx1,dvzx2,dvzx3
      real :: strn, omeg, tprfi

      integer :: ic,jc,kc,ip,jp,kp,im,jm,km,itime

      character(30) :: filnam,filnamxdm,path
      character(5) :: ipfi

      call update_halo(vx,lvlhalo)
      call update_halo(vy,lvlhalo)
      call update_halo(vz,lvlhalo)

!      call AllocateReal3DArray(qtens,1,nx,xstart(2),xend(2),xstart(3),xend(3))
      qtens(:,:,:) = 0.d0

      do ic=xstart(3),xend(3)
      im=ic-1
      ip=ic+1
      do jc=xstart(2),xend(2)
      jm=jc-1
      jp=jc+1
      do kc=2,nxm
      km=kc-1
      kp=kc+1

      dvxx1=(vx(kp,jc,ic)-vx(kc,jc,ic))*udx3m(kc)

      dvxx2=(vx(kc,jp,ic)+vx(kc,jp,ip))- &
            (vx(kc,jm,ic)+vx(kc,jm,ip))*0.25d0*dy

      dvxx3=(vx(kp,jc,ic)+vx(kp,jc,ip))- &
            (vx(km,jc,ic)+vx(km,jc,ip))*0.25d0*dz


      dvyx1=(vy(kc,jc,ip)+vy(kc,jp,ip))- &
            (vy(kc,jc,im)+vy(kc,jp,im))*0.25d0*udx3m(kc)

      dvyx2=(vy(kc,jp,ic)-vy(kc,jc,ic))*dy

      dvyx3=(vy(kp,jc,ic)+vy(kp,jp,ic))- &
            (vy(km,jc,ic)+vy(km,jp,ic))*0.25d0*dz


      dvzx1=(vz(kc,jc,ip)+vz(kp,jc,ip))- &
            (vz(kc,jc,im)+vz(kp,jc,im))*0.25d0*udx3m(kc)

      dvzx2=(vz(kc,jp,ic)+vz(kp,jp,ic))- &
            (vz(kc,jm,ic)+vz(kp,jm,ic))*0.25d0*dy

      dvzx3=(vz(kp,jc,ic)-vz(kc,jc,ic))*dz

      strn=dvxx1**2.0+dvyx2**2.0+dvzx3**2.0+ &
           0.5*((dvyx1+dvxx2)**2.0+ &
                (dvzx1+dvxx3)**2.0+ &
                (dvyx3+dvzx2)**2.0)

      omeg=0.5*((dvyx1-dvxx2)**2.0+ &
                (dvzx1-dvxx3)**2.0+ &
                (dvzx2-dvyx3)**2.0)

      qtens(kc,jc,ic)=0.5*(omeg-strn)

      end do
      end do
      end do

      !-- file name      
      tprfi = 1.d0/tframe
      itime=nint(time*tprfi)
      write(ipfi,'(i5.5)')itime

      path=trim('outputdir/3d/')
      filnam='q'//ipfi//'.h5'
      filnamxdm='outputdir/3d/q'//ipfi//'.xmf'

      !-- Write single precision
      call HdfWriteSingleHalo3D(trim(path)//filnam,real(qtens,selected_real_kind(6, 37)))
      ! call HdfWriteRealHalo3D(trim(path)//filnam,qtens)
!      if(allocated(qtens)) deallocate(qtens)
      
      if (ismaster) then
        open(45,file=filnamxdm,status='unknown')
        rewind(45)
        write(45,'("<?xml version=""1.0"" ?>")')
        write(45,'("<!DOCTYPE Xdmf SYSTEM ""Xdmf.dtd"" []>")')
        write(45,'("<Xdmf Version=""2.0"">")')
        write(45,'("<Domain>")')
        write(45,'("<Grid Name=""basegrid"" GridType=""Uniform"">")')

        write(45,'("<Topology TopologyType=""3DRectMesh"" NumberOfElements=""",i4," ",i4," ",i4,"""/>")') nzm,nym,nx

        write(45,'("<Geometry GeometryType=""VXVYVZ"">")')
        write(45,'("<DataItem Dimensions=""",i4,""" NumberType=""Float"" Precision=""4"" Format=""HDF"">")')nx
        write(45,'("cordin_info.h5:/xc")')
        write(45,'("</DataItem>")')
        write(45,'("<DataItem Dimensions=""",i4,""" NumberType=""Float"" Precision=""4"" Format=""HDF"">")')nym
        write(45,'("cordin_info.h5:/ym")')
        write(45,'("</DataItem>")')
        write(45,'("<DataItem Dimensions=""",i4,""" NumberType=""Float"" Precision=""4"" Format=""HDF"">")')nzm
        write(45,'("cordin_info.h5:/zm")')
        write(45,'("</DataItem>")')
        write(45,'("</Geometry>")')

        write(45,'("<Attribute Name=""Q"" AttributeType=""Scalar"" Center=""Node"">")')
        write(45,'("<DataItem Dimensions=""",i4," ",i4," ",i4,""" NumberType=""Double"" Precision=""2"" Format=""HDF"">")') nzm,nym,nx
        write(45,*) trim(filnam)//':/var'
        write(45,'("</DataItem>")')
        write(45,'("</Attribute>")')

       write(45,'("<Time Value=""",e12.5,""" />")') time

       write(45,'("</Grid>")')
       write(45,'("</Domain>")')
       write(45,'("</Xdmf>")')
       close(45) 
      end if

      end subroutine CalcWriteQ

