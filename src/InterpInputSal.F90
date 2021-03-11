subroutine InterpInputSal
    
    use param
    use input_grids
    use local_arrays, only: sal
    use mgrd_arrays
    use mpih
    use decomp_2d
    use AuxiliaryRoutines
    implicit none

    integer  :: ic,jc,kc, ip,jp,kp, icr,jcr,kcr

    real,dimension(4,4,4) :: qv3 
    real,dimension(4,4) :: qv2
    real,dimension(4) :: qv1

    real, allocatable, dimension(:,:,:) :: salo
    
    sal(:,:,:) = 0.d0

!=========================================================
!     Interpolation of S

    ! Allocate and read in old S
    call AllocateReal3DArray(salo, -1, nxro+1, xs2o-lvlhalo, xe2o+lvlhalo, xs3o-lvlhalo, xe3o+lvlhalo)

    salo(:,:,:) = 0.d0

    call HdfReadContinua(nzro, nyro, nxro, xs2o, xe2o, xs3o, xe3o, 5, &
            salo(1:nxro, xs2o-lvlhalo:xe2o+lvlhalo, xs3o-lvlhalo:xe3o+lvlhalo))

    do ic=xs3o,xe3o
        do jc=xs2o,xe2o
            !-- Boundary points
            !CJH ONLY WORKS FOR CONSTANT SAL BCs CURRENTLY
            salo(0,jc,ic) = salbp(1,jc,ic)
            salo(nxro,jc,ic) = saltp(1,jc,ic)
        end do
    end do

    call update_halo(salo, lvlhalo)

    !-- Interpolate salinity to refined grid
    do ic=xs3o-1,xe3o
        do jc=xs2o-1,xe2o
            do kc=0,nxmro
    
                qv3=salo(kc-1:kc+2,jc-1:jc+2,ic-1:ic+2)
        
                do icr=max(krangs(ic),1),min(krangs(ic+1)-1,nzmr)
                    qv2(:,:) = qv3(:,:,1)*czrs(1,icr)+qv3(:,:,2)*czrs(2,icr) &
                                +qv3(:,:,3)*czrs(3,icr)+qv3(:,:,4)*czrs(4,icr)
                    do jcr=max(jrangs(jc),1),min(jrangs(jc+1)-1,nymr)
                        qv1(:) = qv2(:,1)*cyrs(1,jcr)+qv2(:,2)*cyrs(2,jcr) &
                                +qv2(:,3)*cyrs(3,jcr)+qv2(:,4)*cyrs(4,jcr)
                        do kcr=max(irangs(kc),1),min(irangs(kc+1)-1,nxmr)
                            sal(kcr,jcr,icr) = sum(qv1*cxrs(:,kcr))
                        enddo
                    enddo
                enddo
 
            enddo
        enddo
    enddo

    call DestroyReal3DArray(salo)

    return

end subroutine InterpInputSal