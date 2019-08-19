
module my_module
  use amrex_fort_module, only : amrex_real
  implicit none

contains

  subroutine f5 (a, lo, hi, n) bind(c,name='f5')
    integer, intent(in) :: lo(5), hi(5)
    integer, intent(in), value :: n
    real(amrex_real), intent(inout) :: a(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),lo(4):hi(4),lo(5):hi(5))
    integer :: i1, i2, i3, i4, i5
    print *, "f5: ", lo, hi, n
    do i5 = lo(5), hi(5)
    do i4 = lo(4), hi(4)
    do i3 = lo(3), hi(3)
    do i2 = lo(2), hi(2)
    do i1 = lo(1), hi(1)
       if (a(i1,i2,i3,i4,i5) .ne. 5.) then
          print *, "f5: how did this happend?"
       end if
    end do
    end do
    end do
    end do
    end do
  end subroutine f5

  subroutine f4 (a, lo, hi, n) bind(c,name='f4')
    integer, intent(in) :: lo(4), hi(4)
    integer, intent(in), value :: n
    real(amrex_real), intent(inout) :: a(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),lo(4):hi(4),n)
    integer :: i1, i2, i3, i4, m
    print *, "f4: ", lo, hi, n
    do m = 1, n
    do i4 = lo(4), hi(4)
    do i3 = lo(3), hi(3)
    do i2 = lo(2), hi(2)
    do i1 = lo(1), hi(1)
       if (a(i1,i2,i3,i4,m) .ne. 4.) then
          print *, "f4: how did this happend?"
       end if
    end do
    end do
    end do
    end do
    end do
  end subroutine f4

  subroutine f3 (a, lo, hi) bind(c,name='f3')
    integer, intent(in) :: lo(3), hi(3)
    real(amrex_real), intent(inout) :: a(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
    print *, "f3: ", lo, hi
    a = 1.0_amrex_real;
  end subroutine f3

end module my_module
