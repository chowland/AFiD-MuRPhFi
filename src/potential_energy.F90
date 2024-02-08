!> Module used to compute the available and background potential energy
!! on the fly
module afid_APE
    use param
    use afid_salinity
    implicit none

    integer :: nb=nxmr  !< Number of bins for the sorted density
    real :: zref(nb)    !< Reference vertical coordinate profile Z(rho)
    real :: edges(nb+1) !< Locations of bin edges for density sorting
    real :: bmids(nb)   !< Midpoints of the bins for density
    real :: bmin = -1.0_dp      !< Minimum density
    real :: bmax = 1.0_dp       !< Maximum density
    real :: bw          !< Bin width
    
contains

!> Define the bin edges
subroutine init_bin_edges
    integer :: n

    ! Set the bin edges, spanning the range [-1,1]
    bw = (bmax - bmin)/nb
    edges(1) = bmin
    do n=1,nb
        edges(n+1) = bmin + n*bw
        bmids(n) = 0.5*(edges(n) + edges(n+1))
    end do
end subroutine init_bin_edges

!> Compute a PDF of the salinity field to reconstruct a sorted profile
subroutine sort_salinity
    real :: bpdf(nb)    !< PDF for density
    real :: dxc(nxmr)   !< Vertical grid spacing of refined grid
    integer :: i, j, k, n

    ! Compute grid spacing for volume weighting
    do k=1,nxmr
        dxc(k) = xcr(k+1) - xcr(k)
    end do

    ! Compute the volume-weighted histogram
    do i=xstartr(3),xendr(3)
        do j=xstartr(2),xendr(2)
            do k=1,nxmr
                do n=1,nb
                    if (sal(k,j,i) < edges(n+1) .and. sal(k,j,i) > edges(n-1)) then
                        bpdf(n) = bpdf(n) + dxc(k)
                    end if
                end do
            end do
        end do
    end do

    call MpiAllSumReal1D(bpdf, nb)

    ! Normalize the pdf such that its integral is one
    bpdf = bpdf/sum(bpdf)*nb

    ! Construct a reference vertical coordinate
    zref(1) = bpdf(1)*bw
    do n=2,nb
        zref(n) = zref(n-1) + bpdf(n)*bw
    end do

end subroutine sort_salinity

!> Interpolate the profile Z(rho) stored in `zref`
!! to the value given by vsal
function Zstar(vsal) result(zz)
    real, intent(in) :: vsal !< Salinity value input
    real :: zz
    integer :: n

    if (vsal < bmids(1)) then
        zz = (vsal - bmin)*0.5/bw*zref(1)
    elseif (vsal > bmids(nb)) then
        zz = zref(nb) + (vsal - bmids(nb))*0.5/bw*(1.0 - zref(nb))
    else
        do n=1,nb-1
            if (vsal < bmids(n+1) .and. vsal > bmids(n)) then
                zz = zref(n) + (vsal - bmids(n))/bw*(zref(nb+1) - zref(nb))
            end if
        end do
    end if
end function Zstar


!> Compute the potential energy of the salinity field, along with the
!! background potential energy and the available potential energy
subroutine compute_potential_energy
    real :: rho         !< dimensionless density perturbation
    real :: E_P = 0.0       !< potential energy <rho*g*z>
    real :: E_B = 0.0       !< background potential energy <rho*g*Z(rho)>
    real :: E_A         !< available potential energy E_A = E_P - E_B
    real :: dxc         !< vertical grid spacing (for integration in x)
    real :: inv         !< 1/Lx/ny/nz for normalisation of integration
    logical :: exists
    integer :: i, j, k, io, stat

    do i=xstartr(3),xendr(3)
        do j=xstartr(2),xendr(2)
            do k=1,nxmr
                dxc = xcr(k+1) - xcr(k)
                rho = sal(k,j,i)
                E_P = E_P + rho*xmr(k)*dxc
                E_B = E_B + rho*Zstar(rho)*dxc
            end do
        end do
    end do

    call MpiSumRealScalar(E_P)
    call MpiSumRealScalar(E_B)

    inv = 1.0_dp/alx3/nymr/nzmr
    E_P = E_P*inv
    E_B = E_B*inv
    E_A = E_P - E_B

    if (ismaster) then
        inquire(file="PE.txt", exist=exists)
        if (exists) then
            open(file="PE.txt", newunit=io, iostat=stat, &
                status="old", action="write", position="append")    
        else
            open(file="PE.txt", newunit=io, iostat=stat, &
            status="new", action="write")
        end if
        write(io, *) E_P, E_B, E_A
        close(io)
    end if

end subroutine compute_potential_energy

end module afid_APE