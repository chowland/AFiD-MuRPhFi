!> Module providing the ability to store time-averaged data that
!! is updated at each time step of the simulation
module afid_averaging
    use param, only: nx, nxm, time, dt
    use local_arrays, only: vx, vy, vz, temp
    use decomp_2d, only: xstart, xend

    real, allocatable, dimension(:,:,:) :: vx_tav   !! Time average of x-component of velocity
    real, allocatable, dimension(:,:,:) :: vy_tav   !! Time average of y-component of velocity
    real, allocatable, dimension(:,:,:) :: vz_tav   !! Time average of z-component of velocity
    real, allocatable, dimension(:,:,:) :: te_tav   !! Time average of temperature field

    real :: tav_start   !! Time at which time-averaging begins

contains

!> Allocate the memory for the arrays storing temporal averages
subroutine InitAveragingVariables
    use AuxiliaryRoutines, only: AllocateReal3DArray

    ! Note we don't need ghost cells for these array
    call AllocateReal3DArray(vx_tav,1,nx,xstart(2),xend(2),xstart(3),xend(3))
    call AllocateReal3DArray(vy_tav,1,nx,xstart(2),xend(2),xstart(3),xend(3))
    call AllocateReal3DArray(vz_tav,1,nx,xstart(2),xend(2),xstart(3),xend(3))
    call AllocateReal3DArray(te_tav,1,nx,xstart(2),xend(2),xstart(3),xend(3))

    vx_tav = 0.0
    vy_tav = 0.0
    vz_tav = 0.0
    te_tav = 0.0

    tav_start = time

end subroutine InitAveragingVariables

!> Free memory from arrays used to store temporal averages
!! (at end of simulation)
subroutine DeallocateAveragingVariables
    use AuxiliaryRoutines, only: DestroyReal3DArray

    call DestroyReal3DArray(vx_tav)
    call DestroyReal3DArray(vy_tav)
    call DestroyReal3DArray(vz_tav)
    call DestroyReal3DArray(te_tav)

end subroutine DeallocateAveragingVariables

!> Add dt*f(t) to each array used for time-averaging
subroutine UpdateTemporalAverages
    integer :: i, j, k

    do i=xstart(3),xend(3)
        do j=xstart(2),xend(2)
            do k=1,nxm
                vx_tav(k,j,i) = vx_tav(k,j,i) + dt*vx(k,j,i)
                vy_tav(k,j,i) = vy_tav(k,j,i) + dt*vy(k,j,i)
                vz_tav(k,j,i) = vz_tav(k,j,i) + dt*vz(k,j,i)
                te_tav(k,j,i) = te_tav(k,j,i) + dt*temp(k,j,i)
            end do
        end do
    end do

end subroutine UpdateTemporalAverages

!> Normalize the arrays to complete the time-averaging
!! and write them out to fields/var_mean.h5
subroutine WriteTemporalAverages
    use h5_tools, only: write_3D_array
    real :: tav_length  !! Duration of time averaging interval
    real :: i_tav       !! Inverse of averaging duration
    character(len=30) :: filename

    tav_length = time - tav_start
    i_tav = 1.0/tav_length

    ! Normalize the arrays
    vx_tav = vx_tav*i_tav
    vy_tav = vy_tav*i_tav
    vz_tav = vz_tav*i_tav
    te_tav = te_tav*i_tav

    ! Save the arrays
    filename = 'outputdir/fields/vx_mean.h5'
    call write_3D_array(filename, vx_tav)
    filename = 'outputdir/fields/vy_mean.h5'
    call write_3D_array(filename, vy_tav)
    filename = 'outputdir/fields/vz_mean.h5'
    call write_3D_array(filename, vz_tav)
    filename = 'outputdir/fields/te_mean.h5'
    call write_3D_array(filename, te_tav)
    
end subroutine WriteTemporalAverages

end module afid_averaging