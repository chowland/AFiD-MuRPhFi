!> Module used to compute the available and background potential energy
!! on the fly
module afid_APE
    use param
    use afid_salinity
    implicit none

    integer :: nb       !< Number of bins for the sorted density
    real, allocatable :: zref(:)    !< Reference vertical coordinate profile Z(rho)
    real, allocatable :: edges(:)   !< Locations of bin edges for density sorting
    real, allocatable :: bmids(:)   !< Midpoints of the bins for density
    real :: bmin = -1.d0      !< Minimum density
    real :: bmax = 1.d0       !< Maximum density
    real :: bw          !< Bin width
    integer :: kmin     !< Minimum x-index to measure PE from
    integer :: kmax     !< Maximum x-index to measure PE from
    
contains

!> Define the bin edges
subroutine init_bin_edges
    integer :: n
    logical :: exists
    integer :: io, stat

    nb = nxmr
    allocate(zref(nb+1))
    allocate(edges(nb+1))
    allocate(bmids(nb))

    kmin = nxmr/6 + 1
    kmax = 5*nxmr/6

    if (ismaster) write(*,*) 'kmin, kmax: ', kmin, kmax
    if (ismaster) write(*,*) 'zmin: ', xcr(kmin)
    if (ismaster) write(*,*) 'zmax: ', xcr(kmax+1)

    ! Set the bin edges, spanning the range [-1,1]
    bw = (bmax - bmin)/nb
    edges(1) = bmin
    do n=1,nb
        edges(n+1) = bmin + n*bw
        bmids(n) = 0.5*(edges(n) + edges(n+1))
    end do

    if (ismaster) then
        inquire(file="outputdir/edges.txt", exist=exists)
        if (.not. exists) then
            open(file="outputdir/edges.txt", newunit=io, iostat=stat, &
            status="new", action="write")
            write(io, *) edges
            close(io)
        end if
    end if

end subroutine init_bin_edges

!> Compute a PDF of the salinity field to reconstruct a sorted profile
subroutine sort_salinity
    real :: bpdf(nb)    !< PDF for density
    real :: dxc(nxmr)   !< Vertical grid spacing of refined grid
    integer :: i, j, k, n
    real :: zscale, zmin, zmax
    logical :: exists
    character(30) :: nstat
    character(30) :: fname

    ! Compute grid spacing for volume weighting
    do k=1,nxmr
        dxc(k) = xcr(k+1) - xcr(k)
    end do

    bpdf(:) = 0.0

    ! Compute the volume-weighted histogram
    do i=xstartr(3),xendr(3)
        do j=xstartr(2),xendr(2)
            do k=kmin,kmax
                do n=1,nb
                    if (sal(k,j,i) < edges(n+1) .and. sal(k,j,i) > edges(n)) then
                        bpdf(n) = bpdf(n) + dxc(k)
                    end if
                end do
            end do
        end do
    end do

    call MpiAllSumReal1D(bpdf, nb)

    ! Normalize the pdf such that its integral is one
    bpdf = bpdf/sum(bpdf)*nb/(bmax - bmin)

    ! Construct a reference vertical coordinate
    zmax = xcr(kmax+1)
    zmin = xcr(kmin)
    zscale = bw*(zmax - zmin)
    zref(1) = zmax
    do n=1,nb
        zref(n+1) = zref(n) - bpdf(n)*zscale
    end do

    if (ismaster) then
        write(nstat,"(i5.5)") nint(time/tout)
        fname = "outputdir/Zr.h5"
        inquire(file=trim(fname), exist=exists)
        if (exists) then
            call HdfSerialWriteReal1D(trim(nstat), trim(fname), zref, nb+1)
        else
            call HdfCreateBlankFile(trim(fname))
            call HdfSerialWriteReal1D(trim(nstat), trim(fname), zref, nb+1)
        end if
    end if

end subroutine sort_salinity

!> Numerically integrate the function f over the grid x
!! using the trapezoidal rule and return the result as I
function trapz(f, x) result(I)
    real, intent(in) :: f(:)    !< function to integrate
    real, intent(in) :: x(:)    !< x-coordinate vector
    real :: I       !< Integral result
    integer :: n, k

    n = size(x)
    I = 0.0
    do k=1,n-1
        I = I + 0.5*(f(k) + f(k+1))*(x(k+1) - x(k))
    end do

end function trapz


!> Compute the potential energy of the salinity field, along with the
!! background potential energy and the available potential energy
subroutine compute_potential_energy
    real :: rho         !< dimensionless density perturbation
    real :: E_P = 0.0       !< potential energy <rho*g*z>
    real :: E_B = 0.0       !< background potential energy <rho*g*Z(rho)>
    real :: E_B2(nb+1)       !< background potential energy <rho_b*g*z>
    real :: E_A         !< available potential energy E_A = E_P - E_B
    real :: dxc         !< vertical grid spacing (for integration in x)
    real :: inv         !< 1/Lx/ny/nz for normalisation of integration
    logical :: exists
    integer :: i, j, k, io, stat

    call sort_salinity

    E_P = 0.0

    do i=xstartr(3),xendr(3)
        do j=xstartr(2),xendr(2)
            do k=kmin,kmax
                dxc = xcr(k+1) - xcr(k)
                rho = sal(k,j,i)
                E_P = E_P + rho*xmr(k)*dxc
            end do
        end do
    end do

    do k=1,nb+1
        E_B2(k) = edges(k)*zref(k)
    end do

    E_B = -trapz(E_B2, zref)

    call MpiSumRealScalar(E_P)

    inv = 1.d0/nymr/nzmr
    E_P = E_P*inv
    E_A = E_P - E_B

    if (ismaster) then
        inquire(file="outputdir/PE.txt", exist=exists)
        if (exists) then
            open(file="outputdir/PE.txt", newunit=io, iostat=stat, &
                status="old", action="write", position="append")    
        else
            open(file="outputdir/PE.txt", newunit=io, iostat=stat, &
            status="new", action="write")
        end if
        write(io, *) E_P, E_B, E_A
        close(io)
    end if

end subroutine compute_potential_energy

end module afid_APE