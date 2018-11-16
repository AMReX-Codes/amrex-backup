
module my_module

  use iso_c_binding, only : c_int
  implicit none

contains

  ! AMREX_CUDA_FORT_DEVICE = attributes(device)
  AMREX_CUDA_FORT_DEVICE subroutine f(x) bind(c,name='f')
    integer(c_int), intent(inout) :: x
    print *, 'in fortran'
    x = x * 2
  end subroutine f

  ! We would like to have AMREX_CUDA_FORT_HOST_DEVICE that could generate
  ! two versions of the function, one for host and one for device.

end module my_module
