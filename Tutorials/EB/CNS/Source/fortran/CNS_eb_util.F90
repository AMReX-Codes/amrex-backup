module CNS_eb_util_module

  use amrex_ebcellflag_module, only : get_neighbor_cells
  use eb_stencil_types_module, only : eb_bndry_geom

  implicit none

contains

  subroutine set_eb_nbr(ebg, Nebg, flag, fglo, fghi) &
       bind(C,name="set_eb_nbr")
    integer,intent(in) :: Nebg, fglo(0:2),fghi(0:2)
    integer,intent(in) :: flag(fglo(0):fghi(0),fglo(1):fghi(1),fglo(2):fghi(2))
    type(eb_bndry_geom),intent(inout) :: ebg(0:Nebg-1)
    integer :: L
    integer nbr(-1:1,-1:1,-1:1)
    do L = 0, Nebg-1
        call get_neighbor_cells(flag(ebg(L)%iv(0),ebg(L)%iv(1),ebg(L)%iv(2)),ebg(L)%nbr)
    enddo
  end subroutine set_eb_nbr

end module CNS_eb_util_module

