!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: CalcMaxCFL.F90                                 !
!    CONTAINS: subroutine CalcMaxCFL                      !
!                                                         ! 
!    PURPOSE: Compute the maximum value of the local CFL  !
!     stability condition for the explicit terms          !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine CalcMaxCFL(cflm,cflmr)
    use param
    use local_arrays, only: vx,vy,vz,vxr,vyr,vzr
    use decomp_2d
    use mpih
    implicit none
    real,intent(out) :: cflm, cflmr
    integer :: j,k,jp,kp,i,ip
    real :: qcf
    cflm=0.00000001d0
  
    !
    !$OMP  PARALLEL DO &
    !$OMP   DEFAULT(none) &
    !$OMP   SHARED(xstart,xend,nxm,vz,vy,vx) &
    !$OMP   SHARED(dz,dy,udx3m) &
    !$OMP   PRIVATE(i,j,k,ip,jp,kp,qcf) &
    !$OMP   REDUCTION(max:cflm)
    do i=xstart(3),xend(3)
        ip=i+1
        do j=xstart(2),xend(2)
            jp=j+1
            do k=1,nxm
                kp=k+1
                qcf=( abs((vz(k,j,i)+vz(k,j,ip))*0.5d0*dz) &
                     +abs((vy(k,j,i)+vy(k,jp,i))*0.5d0*dy) &
                     +abs((vx(k,j,i)+vx(kp,j,i))*0.5d0*udx3m(k)))

                cflm = max(cflm,qcf)

            end do
        end do
    end do

    !$OMP END PARALLEL DO
    call MpiAllMaxRealScalar(cflm)

  
    if (salinity .or. multiRes_Temp) then
        ! Refined mesh
        cflmr = 1.d-8
        do i=xstartr(3),xendr(3)
            ip = i + 1
            do j=xstartr(2),xendr(2)
                jp = j + 1
                do k=1,nxmr
                    kp = k + 1
                    qcf=( abs((vzr(k,j,i) + vzr(k,j,ip))*0.5d0*dzr) &
                    +abs((vyr(k,j,i) + vyr(k,jp,i))*0.5d0*dyr) &
                    +abs((vxr(k,j,i) + vxr(kp,j,i))*0.5d0*udx3mr(k)))
                    
                    cflmr = max(cflmr,qcf)

                end do
            end do
        end do
        call MpiAllMaxRealScalar(cflmr)

    else

        cflmr = cflm
    end if

    return
end subroutine CalcMaxCFL