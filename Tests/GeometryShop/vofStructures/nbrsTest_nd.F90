module nbrsTest_nd_module

  use amrex_fort_module, only : amrex_real, dim=>bl_spacedim
  use nbrs_test_module, only : nbr_sten, face_sten

  implicit none

contains

  subroutine fill_redist_stencil(lo, hi, slo, shi, sten, Nsten, vf,vf_lo,vf_hi) bind(C,name="fill_redist_stencil")

    integer,          intent(in) ::  lo(0:2),  hi(0:2)
    integer,          intent(in) :: slo(0:2), shi(0:2)
    integer,          intent(in) :: Nsten
    type(nbr_sten),intent(inout) :: sten(0:Nsten-1)
    integer,          intent(in) :: vf_lo(0:2), vf_hi(0:2)
    real(amrex_real), intent(in) :: vf( vf_lo(0):vf_hi(0),vf_lo(1):vf_hi(1),vf_lo(2):vf_hi(2))
    integer :: i,j,k,n,ii,jj,kk,iii,jjj,kkk
    real(amrex_real) :: kappa_tot

    do n = 0, Nsten-1

       i = sten(n) % iv(0)
       j = sten(n) % iv(1)
       k = sten(n) % iv(2)

       if (i.ge.lo(0) .and. i.le.hi(0) &
            .and. j.ge.lo(1) .and. j.le.hi(1) &
            .and. k.ge.lo(2) .and. k.le.hi(2) ) then

          kappa_tot = 0.d0
          do kk=slo(2),shi(2)
             kkk = k+kk
             do jj=slo(1),shi(1)
                jjj = j+jj
                do ii=slo(0),shi(0)
                   iii = i+ii
                   if (sten(n)%val(ii,jj,kk) .gt. 0) then
                      kappa_tot = kappa_tot + vf(iii,jjj,kkk)
                   endif
                enddo
             enddo
          enddo

          kappa_tot = 1.d0 / kappa_tot

          do kk=slo(2),shi(2)
             kkk = k+kk
             do jj=slo(1),shi(1)
                jjj = j+jj
                do ii=slo(0),shi(0)
                   iii = i+ii
                   if (sten(n)%val(ii,jj,kk) .gt. 0) then
                      sten(n) % val(ii,jj,kk) = vf(iii,jjj,kkk) * kappa_tot
                   else
                      sten(n) % val(ii,jj,kk) = 0.d0
                   endif
                enddo
             enddo
          enddo

       endif

    enddo

  end subroutine fill_redist_stencil

  subroutine apply_redist_stencil(lo, hi, slo, shi, sten, Nsten, vin, vin_lo, vin_hi, &
    vout, vout_lo, vout_hi) bind(C,name="apply_redist_stencil")

    integer,          intent(in   ) ::  lo(0:2),  hi(0:2)
    integer,          intent(in   ) :: slo(0:2), shi(0:2)
    integer,          intent(in   ) :: Nsten
    type(nbr_sten),   intent(in   ) :: sten(0:Nsten-1)
    integer,          intent(in   ) ::  vin_lo(0:2),  vin_hi(0:2)
    integer,          intent(in   ) :: vout_lo(0:2), vout_hi(0:2)
    real(amrex_real), intent(in   ) ::  vin( vin_lo(0):vin_hi(0),  vin_lo(1):vin_hi(1),  vin_lo(2):vin_hi(2))
    real(amrex_real), intent(inout) :: vout(vout_lo(0):vout_hi(0),vout_lo(1):vout_hi(1),vout_lo(2):vout_hi(2))
    integer :: i,j,k,n,ii,jj,kk,iii,jjj,kkk
    real(amrex_real) :: vtot

    do n = 0, Nsten-1

       i = sten(n) % iv(0)
       j = sten(n) % iv(1)
       k = sten(n) % iv(2)

       if (i.ge.lo(0) .and. i.le.hi(0) &
            .and. j.ge.lo(1) .and. j.le.hi(1)&
            .and. k.ge.lo(2) .and. k.le.hi(2) ) then

          vtot = 0.d0
          do kk=slo(2),shi(2)
             kkk = k+kk
             do jj=slo(1),shi(1)
                jjj = j+jj
                do ii=slo(0),shi(0)
                   iii = i+ii
                   vtot = vtot + sten(n)%val(ii,jj,kk) * vin(iii,jjj,kkk)
                enddo
             enddo
          enddo
          vout(i,j,k) = vtot

       endif

    enddo
    
  end subroutine apply_redist_stencil

  subroutine fill_flux_interp_stencil(lo, hi, slo, shi, sten, Nsten, idir, cent, cent_lo, cent_hi) &
       bind(C,name="fill_flux_interp_stencil")

    integer,          intent(in)  ::  lo(0:2),  hi(0:2)
    integer,          intent(in)  :: slo(0:2), shi(0:2)
    integer,          intent(in)  :: Nsten, idir
    type(face_sten),intent(inout) :: sten(0:Nsten-1)
    integer,          intent(in)  :: cent_lo(0:2), cent_hi(0:2)
    real(amrex_real), intent(in)  :: cent( cent_lo(0):cent_hi(0),cent_lo(1):cent_hi(1),cent_lo(2):cent_hi(2),0:2)
    integer :: i,j,k,n,in, jn, kn
    real(amrex_real) :: cx,cy,cz,acx,acy,acz

    if (idir .eq. 0) then
       do n = 0, Nsten-1

          i = sten(n) % iv(0)
          j = sten(n) % iv(1)
          k = sten(n) % iv(2)

          if (i.ge.lo(0) .and. i.le.hi(0) &
               .and. j.ge.lo(1) .and. j.le.hi(1) &
               .and. k.ge.lo(2) .and. k.le.hi(2) ) then

#if BL_SPACEDIM == 2
             sten(n)%val(:) = 0.d0
             cy = cent(i,j,k,1)
             jn = SIGN(1.d0, cy)
             acy = ABS(cy)

             sten(n)%val(0) = 1.d0 - acy
             sten(n)%val(jn) = acy
#else
             sten(n)%val(:,:) = 0.d0
             cy = cent(i,j,k,1)
             jn = SIGN(1.d0, cy)
             acy = ABS(cy)
             cz = cent(i,j,k,2)
             kn = SIGN(1.d0, cz)
             acz = ABS(cz)

             sten(n)%val( 0, 0) = (1.d0 - acy) * (1.d0 - acz)
             sten(n)%val(jn, 0) =     acy      * (1.d0 - acz)
             sten(n)%val( 0,kn) = (1.d0 - acy) *     acz
             sten(n)%val(jn,kn) =     acy      *     acz
#endif

          endif

       enddo
    else if (idir .eq. 1) then
       do n = 0, Nsten-1

          i = sten(n) % iv(0)
          j = sten(n) % iv(1)
          k = sten(n) % iv(2)

          if (i.ge.lo(0) .and. i.le.hi(0) &
               .and. j.ge.lo(1) .and. j.le.hi(1) &
               .and. k.ge.lo(2) .and. k.le.hi(2) ) then
             
#if BL_SPACEDIM == 2
             sten(n)%val(:) = 0.d0
             cx = cent(i,j,k,0)
             in = SIGN(1.d0, cx)
             acx = ABS(cx)

             sten(n)%val(0) = 1.d0 - acx
             sten(n)%val(in) = acx
#else
             sten(n)%val(:,:) = 0.d0
             cx = cent(i,j,k,0)
             in = SIGN(1.d0, cx)
             acx = ABS(cx)
             cz = cent(i,j,k,2)
             kn = SIGN(1.d0, cz)
             acz = ABS(cz)

             sten(n)%val( 0, 0) = (1.d0 - acx) * (1.d0 - acz)
             sten(n)%val(in, 0) =     acx      * (1.d0 - acz)
             sten(n)%val( 0,kn) = (1.d0 - acx) *     acz
             sten(n)%val(in,kn) =     acx      *     acz
#endif
          endif

       enddo
#if BL_SPACEDIM == 3
    else if (idir .eq. 2) then
       do n = 0, Nsten-1

          i = sten(n) % iv(0)
          j = sten(n) % iv(1)
          k = sten(n) % iv(2)

          if (i.ge.lo(0) .and. i.le.hi(0) &
               .and. j.ge.lo(1) .and. j.le.hi(1) &
               .and. k.ge.lo(2) .and. k.le.hi(2) ) then
             
             sten(n)%val(:,:) = 0.d0
             cx = cent(i,j,k,0)
             in = SIGN(1.d0, cx)
             acx = ABS(cx)
             cy = cent(i,j,k,1)
             jn = SIGN(1.d0, cy)
             acy = ABS(cy)

             sten(n)%val( 0, 0) = (1.d0 - acx) * (1.d0 - acy)
             sten(n)%val(in, 0) =     acx      * (1.d0 - acy)
             sten(n)%val( 0,jn) = (1.d0 - acx) *     acy
             sten(n)%val(in,jn) =     acx      *     acy
          endif

       enddo
#endif
    endif

  end subroutine fill_flux_interp_stencil

  subroutine apply_flux_interp_stencil(lo, hi, slo, shi, sten, Nsten, idir, vin, vin_lo, vin_hi, &
    vout, vout_lo, vout_hi) bind(C,name="apply_flux_interp_stencil")

    integer,          intent(in   ) ::  lo(0:2),  hi(0:2)
    integer,          intent(in   ) :: slo(0:2), shi(0:2)
    integer,          intent(in   ) :: Nsten, idir
    type(face_sten),  intent(in   ) :: sten(0:Nsten-1)
    integer,          intent(in   ) ::  vin_lo(0:2),  vin_hi(0:2)
    integer,          intent(in   ) :: vout_lo(0:2), vout_hi(0:2)
    real(amrex_real), intent(in   ) ::  vin( vin_lo(0):vin_hi(0),  vin_lo(1):vin_hi(1),  vin_lo(2):vin_hi(2))
    real(amrex_real), intent(inout) :: vout(vout_lo(0):vout_hi(0),vout_lo(1):vout_hi(1),vout_lo(2):vout_hi(2))
    integer :: i,j,k,n,ii,jj,kk,iii,jjj,kkk
    real(amrex_real) :: vtot


    if (idir.eq.0) then

       do n = 0, Nsten-1
          i = sten(n) % iv(0)
          j = sten(n) % iv(1)
          k = sten(n) % iv(2)

          if (i.ge.lo(0) .and. i.le.hi(0) &
               .and. j.ge.lo(1) .and. j.le.hi(1)&
               .and. k.ge.lo(2) .and. k.le.hi(2) ) then
             
             vtot = 0.d0
#if BL_SPACEDIM == 2
             do jj=slo(1),shi(1)
                jjj = j+jj
                vtot = vtot + sten(n)%val(jj) * vin(i,jjj,k)
             enddo
#else
             do kk=slo(2),shi(2)
                kkk = k+kk
                do jj=slo(1),shi(1)
                   jjj = j+jj
                   vtot = vtot + sten(n)%val(jj,kk) * vin(i,jjj,kkk)
                enddo
             enddo
#endif
             vout(i,j,k) = vtot
          endif
       enddo
    
    else if (idir.eq.1) then

       do n = 0, Nsten-1
          i = sten(n) % iv(0)
          j = sten(n) % iv(1)
          k = sten(n) % iv(2)

          if (i.ge.lo(0) .and. i.le.hi(0) &
               .and. j.ge.lo(1) .and. j.le.hi(1)&
               .and. k.ge.lo(2) .and. k.le.hi(2) ) then
             
             vtot = 0.d0
#if BL_SPACEDIM == 2
             do ii=slo(0),shi(0)
                iii = i + ii
                vtot = vtot + sten(n)%val(ii) * vin(iii,j,k)
             enddo
#else
             do kk=slo(2),shi(2)
                kkk = k+kk
                do ii=slo(0),shi(0)
                   iii = i + ii
                   vtot = vtot + sten(n)%val(ii,kk) * vin(iii,j,kkk)
                enddo
             enddo
#endif
             vout(i,j,k) = vtot
          endif
       enddo
    
#if BL_SPACEDIM == 3
    else

       do n = 0, Nsten-1
          i = sten(n) % iv(0)
          j = sten(n) % iv(1)
          k = sten(n) % iv(2)

          if (i.ge.lo(0) .and. i.le.hi(0) &
               .and. j.ge.lo(1) .and. j.le.hi(1)&
               .and. k.ge.lo(2) .and. k.le.hi(2) ) then
             
             vtot = 0.d0

             do jj=slo(1),shi(1)
                jjj = j+jj
                do ii=slo(0),shi(0)
                   iii = i+ii
                   vtot = vtot + sten(n)%val(ii,jj) * vin(iii,jjj,k)
                enddo
             enddo
             vout(i,j,k) = vtot
          endif
       enddo
#endif
    endif
    
  end subroutine apply_flux_interp_stencil


end module nbrsTest_nd_module
