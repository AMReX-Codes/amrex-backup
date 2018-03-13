module eb_stencil_types_module

  use amrex_fort_module, only : amrex_real, dim=>bl_spacedim
  implicit none

#if BL_SPACEDIM == 2

  type, bind(c) :: vol_sten
     real(amrex_real) :: val(-1:1,-1:1)
     integer :: iv(0:dim-1)
  end type vol_sten

  type, bind(c) :: face_sten
     real(amrex_real) :: val(-1:1)
     integer :: iv(0:dim-1)
  end type face_sten

  type, bind(c) :: eb_bndry_sten
     real(amrex_real) :: val(-1:1,-1:1)
     real(amrex_real) :: bcval
     integer          :: iv(0:dim-1)
     integer          :: iv_base(0:dim-1)
  end type eb_bndry_sten

#elif BL_SPACEDIM == 3

  type, bind(c) :: vol_sten
     real(amrex_real) :: val(-1:1,-1:1,-1:1)
     integer :: iv(0:dim-1)
  end type vol_sten

  type, bind(c) :: face_sten
     real(amrex_real) :: val(-1:1,-1:1)
     integer :: iv(0:dim-1)
  end type face_sten

  type, bind(c) :: eb_bndry_sten
     real(amrex_real) :: val(-1:1,-1:1, -1:1)
     real(amrex_real) :: bcval
     integer          :: iv(0:dim-1)
     integer          :: iv_base(0:dim-1)
  end type eb_bndry_sten

#endif

  type, bind(c) :: eb_bndry_geom
     real(amrex_real) :: eb_normal(dim)
     real(amrex_real) :: eb_centroid(dim)
     real(amrex_real) :: eb_area
     integer :: iv(0:dim-1)
#if BL_SPACEDIM == 2
     integer :: nbr(-1:1,-1:1)
#elif BL_SPACEDIM == 3
     integer :: nbr(-1:1,-1:1,-1:1)
#endif
  end type eb_bndry_geom

end module eb_stencil_types_module
