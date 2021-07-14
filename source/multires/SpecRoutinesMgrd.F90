subroutine CalcFourierCoef(var,fouvar)
    use, intrinsic :: iso_c_binding
    use param
    use fftw_params
    use decomp_2d
    use decomp_2d_fft
    use mpih
    implicit none
    integer :: i,j,k,info
    integer :: phpiv(nxmr)
    integer :: nymh
    real,intent(in),dimension(1:nxmr,xstartr(2):xendr(2),xstartr(3):xendr(3)) :: var
    complex,intent(inout),dimension(sp%xst(1):sp%xen(1), &
                                    sp%xst(2):sp%xen(2), &
                                    sp%xst(3):sp%xen(3)) :: fouvar
    
    type(fftw_iodim),dimension(1) :: iodim
    type(fftw_iodim),dimension(2) :: iodim_howmany
    
    !RO   Allocate variables for FFT transform
    
    call CreateFFTTmpArrays
    
    nymh=nymr/2+1
    
    call transpose_x_to_y(var,ry1,ph)
    
    !RO   Plan FFT transforms if not planned previously
    
    if (.not.plannedr) then
        iodim(1)%n=nzmr
        iodim(1)%is=(sp%zen(1)-sp%zst(1)+1)*(sp%zen(2)-sp%zst(2)+1)
        iodim(1)%os=(sp%zen(1)-sp%zst(1)+1)*(sp%zen(2)-sp%zst(2)+1)
        iodim_howmany(1)%n=(sp%zen(1)-sp%zst(1)+1)
        iodim_howmany(1)%is=1
        iodim_howmany(1)%os=1
        iodim_howmany(2)%n=(sp%zen(2)-sp%zst(2)+1)
        iodim_howmany(2)%is=(sp%zen(1)-sp%zst(1)+1)
        iodim_howmany(2)%os=(sp%zen(1)-sp%zst(1)+1)
        fwd_guruplan_zr=fftw_plan_guru_dft(1,iodim, &
                2,iodim_howmany,cz1,cz1, &
                FFTW_FORWARD,FFTW_ESTIMATE)
        iodim(1)%n=nzmr
        bwd_guruplan_zr=fftw_plan_guru_dft(1,iodim, &
                2,iodim_howmany,cz1,cz1, &
                FFTW_BACKWARD,FFTW_ESTIMATE)
    
        if (.not.c_associated(bwd_guruplan_z)) then
            if (ismaster) print*,'Failed to create guru plan. You should'
            if (ismaster) print*,'link with FFTW3 before MKL'
            if (ismaster) print*,'Please check linking order.'
            call MPI_Abort(MPI_COMM_WORLD,1,info)
        endif
    
        iodim(1)%n=nymr
        iodim(1)%is=ph%yen(1)-ph%yst(1)+1
        iodim(1)%os=sp%yen(1)-sp%yst(1)+1
        iodim_howmany(1)%n=(ph%yen(1)-ph%yst(1)+1)
        iodim_howmany(1)%is=1
        iodim_howmany(1)%os=1
        iodim_howmany(2)%n=(ph%yen(3)-ph%yst(3)+1)
        iodim_howmany(2)%is=(ph%yen(1)-ph%yst(1)+1) &
                *(ph%yen(2)-ph%yst(2)+1)
        iodim_howmany(2)%os=(sp%yen(1)-sp%yst(1)+1) &
                *(sp%yen(2)-sp%yst(2)+1)
        fwd_guruplan_yr=fftw_plan_guru_dft_r2c(1,iodim, &
                2,iodim_howmany,ry1,cy1, &
                FFTW_ESTIMATE)
    
        iodim(1)%n=nymr
        iodim(1)%is=sp%yen(1)-sp%yst(1)+1
        iodim(1)%os=ph%yen(1)-ph%yst(1)+1
        iodim_howmany(1)%n=(sp%yen(1)-sp%yst(1)+1)
        iodim_howmany(1)%is=1
        iodim_howmany(1)%os=1
        iodim_howmany(2)%n=(sp%yen(3)-sp%yst(3)+1)
        iodim_howmany(2)%is=(sp%yen(1)-sp%yst(1)+1) &
                *(sp%yen(2)-sp%yst(2)+1)
        iodim_howmany(2)%os=(ph%yen(1)-ph%yst(1)+1) &
                *(ph%yen(2)-ph%yst(2)+1)
        bwd_guruplan_yr=fftw_plan_guru_dft_c2r(1,iodim, &
                2,iodim_howmany,cy1,ry1, &
                FFTW_ESTIMATE)
        plannedr=.true.
    endif
    
    call dfftw_execute_dft_r2c(fwd_guruplan_y,ry1,cy1)
    
    call transpose_y_to_z(cy1,cz1,sp)
    
    call dfftw_execute_dft(fwd_guruplan_z,cz1,cz1)
    
    !EP   Normalize. FFT does not do this
    cz1 = cz1 / (nzmr*nymr)
    
    call transpose_z_to_x(cz1,fouvar,sp)
    
    call DestroyFFTTmpArrays
    
    return
end subroutine CalcFourierCoef
!
!***********************************************************************
subroutine CalcPowerSpecMgrd(var, idx, specy, specz)
    use, intrinsic :: iso_c_binding
    use param
    use fftw_params
    use decomp_2d
    use decomp_2d_fft
    use mpih
    implicit none

    real,intent(in),dimension(1:nxmr,xstartr(2):xendr(2),xstartr(3):xendr(3)) :: var
    integer,intent(in),dimension(5) :: idx
    real,intent(inout),dimension(5,1:int(nymr/2+1)) :: specy
    real,intent(inout),dimension(5,1:int(nzmr/2+1)) :: specz
    real :: re_uhat,re_vhat,im_uhat,im_vhat
    real :: re_spec,im_spec
    integer :: i,j,k,izk,jyk,niz,nizk,kk

!-- output is specy and specz
    specy(:,:)=0.d0
    specz(:,:)=0.d0
    
    allocate(fouvar1(sp%xst(1):sp%xen(1), &
                    sp%xst(2):sp%xen(2), &
                    sp%xst(3):sp%xen(3)))
    
    !-- Calculate DFT for var
    call CalcFourierCoef(var,fouvar1)
        
    !-- Calculate the spectra in y and z (periodic dirs)
    do i=sp%xst(3),sp%xen(3)
        do j=sp%xst(2),sp%xen(2)
            do kk=1,5
                k = idx(kk)

                re_uhat= real(fouvar1(k,j,i))
                im_uhat=aimag(fouvar1(k,j,i))

                re_spec = re_uhat**2 + im_uhat**2
    
                izk = i-1 !-- i counter for specz
                jyk = j-1 !-- j counter for specy
                niz = mod(nzmr-(i-1),nzmr)+1
                nizk = niz-1

                !-- specy
                if ( i.eq.1 .and. j.eq.1 ) then
                else
                    specy(k,jyk) = specy(k,jyk) + re_spec
                end if

                !-- specz
                if ( j.eq.1 ) then
                    if ( i.gt.1 .and. i.le.(nzmr/2+1) ) then
                        specz(k,izk) = specz(k,izk) + re_spec
                    end if
                else if( j.eq.(nymr/2+1) ) then
                    if ( i.ge.1 .and. i.le.(nzmr/2+1) ) then
                        specz(k,izk) = specz(k,izk) + re_spec
                    end if
                else
                    if ( i.ge.1 .and. i.le.(nzmr/2+1) ) then
                        specz(k,izk) = specz(k,izk) + re_spec
                    end if
                    if ( niz.ge.1 .and. niz.le.(nzmr/2+1) ) then
                        specz(k,nizk) = specz(k,nizk) + re_spec
                    end if
                end if
            end do
        end do
    end do

    if(allocated(fouvar1)) deallocate(fouvar1)

!-- Save the files
    do k=1,5
        call MpiSumReal1D(specy(k,:),(nymr/2+1))
        call MpiSumReal1D(specz(k,:),(nzmr/2+1))
    end do
    
    return
end subroutine CalcPowerSpecMgrd