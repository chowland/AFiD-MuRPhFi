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
    
contains

!> Define the bin edges
subroutine init_bin_edges
    integer :: n

    nb = nxmr
    allocate(zref(nb+1))
    allocate(edges(nb+1))
    allocate(bmids(nb))

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
    real :: zscale
    logical :: exists
    integer :: io, stat

    ! Compute grid spacing for volume weighting
    do k=1,nxmr
        dxc(k) = xcr(k+1) - xcr(k)
    end do

    bpdf(:) = 0.0

    ! Compute the volume-weighted histogram
    do i=xstartr(3),xendr(3)
        do j=xstartr(2),xendr(2)
            do k=1,nxmr
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
    zscale = bw*alx3
    zref(1) = alx3
    do n=1,nb
        zref(n+1) = zref(n) - bpdf(n)*zscale
    end do

    if (ismaster) then
        ! call HdfSerialWriteReal1D()
        inquire(file="Zr.txt", exist=exists)
        if (.not. exists) then
            open(file="Zr.txt", newunit=io, iostat=stat, &
            status="new", action="write")
            write(io, *) zref
            close(io)
        end if
        inquire(file="pdf.txt", exist=exists)
        if (.not. exists) then
            open(file="pdf.txt", newunit=io, iostat=stat, &
            status="new", action="write")
            write(io, *) bpdf
            close(io)
        end if
    end if

end subroutine sort_salinity

!> Interpolate the profile Z(rho) stored in `zref`
!! to the value given by vsal
function Zstar(vsal) result(zz)
    real, intent(in) :: vsal !< Salinity value input
    real :: zz
    integer :: n

    do n=1,nb
        if (vsal < edges(n+1) .and. vsal > edges(n)) then
            zz = zref(n) + (vsal - edges(n))*(zref(n+1) - zref(n))/bw
        end if
    end do

end function Zstar

!> Interpolate the profile C(xx) stored implicitly by `zref`
!! to the height given by xx
function Cstar(xx) result(C)
    real, intent(in) :: xx !< Vertical position
    real :: C
    integer :: n

    do n=1,nb
        if (xx > zref(n+1) .and. xx < zref(n)) then
            C = edges(n) + (xx - zref(n))*bw/(zref(n+1) - zref(n))
        end if
    end do
end function Cstar


!> Compute the potential energy of the salinity field, along with the
!! background potential energy and the available potential energy
subroutine compute_potential_energy
    real :: rho         !< dimensionless density perturbation
    real :: E_P = 0.0       !< potential energy <rho*g*z>
    real :: E_B = 0.0       !< background potential energy <rho*g*Z(rho)>
    real :: E_B2       !< background potential energy <rho_b*g*z>
    real :: E_A         !< available potential energy E_A = E_P - E_B
    real :: E_A2         !< available potential energy E_A = E_P - E_B2
    real :: dxc         !< vertical grid spacing (for integration in x)
    real :: inv         !< 1/Lx/ny/nz for normalisation of integration
    real :: Cs(nxmr), Zs(nxmr)
    logical :: exists
    integer :: i, j, k, io, stat
    character(5) :: nstat

    call sort_salinity

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

    Cs(:) = 0.0

    E_B2 = 0.0
    do k=1,nxmr
        Cs(k) = Cstar(xmr(k))
        Zs(k) = Zstar(-tanh(xmr(k) - 0.5*alx3))
        E_B2 = E_B2 + xmr(k)*(xcr(k+1) - xcr(k))*Cs(k)
    end do
    E_B2 = E_B2/alx3
    if (ismaster) then
        write(nstat,"(i5.5)") nint(time/tout)
        inquire(file="Cs.h5", exist=exists)
        if (exists) then
            call HdfSerialWriteReal1D(nstat, "Cs.h5", Cs, nxmr)
        else
            call HdfCreateBlankFile("Cs.h5")
            call HdfSerialWriteReal1D(nstat, "Cs.h5", Cs, nxmr)
        end if

        inquire(file="Cs.txt", exist=exists)
        if (.not. exists) then
            open(file="Cs.txt", newunit=io, iostat=stat, &
            status="new", action="write")
            write(io, *) Cs(:)
            close(io)
        end if

        inquire(file="Zs.txt", exist=exists)
        if (.not. exists) then
            open(file="Zs.txt", newunit=io, iostat=stat, &
            status="new", action="write")
            write(io, *) Zs
            close(io)
        end if
    end if

    call MpiSumRealScalar(E_P)
    call MpiSumRealScalar(E_B)

    inv = 1.d0/alx3/nymr/nzmr
    E_P = E_P*inv
    E_B = E_B*inv
    ! E_B2 = E_B2*inv
    E_A = E_P - E_B
    E_A2 = E_P - E_B2

    if (ismaster) then
        inquire(file="PE.txt", exist=exists)
        if (exists) then
            open(file="PE.txt", newunit=io, iostat=stat, &
                status="old", action="write", position="append")    
        else
            open(file="PE.txt", newunit=io, iostat=stat, &
            status="new", action="write")
        end if
        write(io, *) E_P, E_B, E_B2, E_A, E_A2
        close(io)
    end if

end subroutine compute_potential_energy

end module afid_APE