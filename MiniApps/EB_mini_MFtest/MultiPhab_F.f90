subroutine init_mf_val(lo, hi, f) bind(C, name = "init_mf_val")

    implicit none

    integer, intent(in) :: lo(3), hi(3)

    double precision, intent(inout) :: f(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))

    double precision, parameter :: small = 1.d-8

    integer          :: i, j, k

    double precision :: val

    val = 0.0
    do k = lo(3)+1, hi(3)-1
       do j = lo(2)+1, hi(2)-1
          do i = lo(1), hi(1)
          f(i,j,k) = 0.0

          if( i==1 .and. j==1) f(i,j,k) = 1.0
          if( i==2 .and. j==1) f(i,j,k) = 2.0
          if( i==3 .and. j==1) f(i,j,k) = 3.0
          if( i==4 .and. j==1) f(i,j,k) = 4.0
          if( i==5 .and. j==1) f(i,j,k) = 5.0

          if( i==1 .and. j==2) f(i,j,k) = 6.0
          if( i==2 .and. j==2) f(i,j,k) = -7.0
          if( i==3 .and. j==2) f(i,j,k) = 8.0
          if( i==4 .and. j==2) f(i,j,k) = 9.0
          if( i==5 .and. j==2) f(i,j,k) = 10.0

          if( i==1 .and. j==3) f(i,j,k) = 11.0
          if( i==2 .and. j==3) f(i,j,k) = -12.0
          if( i==3 .and. j==3) f(i,j,k) = 13.0
          if( i==4 .and. j==3) f(i,j,k) = 14.0
          if( i==5 .and. j==3) f(i,j,k) = 15.0

          if( i==1 .and. j==4) f(i,j,k) = 16.0
          if( i==2 .and. j==4) f(i,j,k) = 17.0
          if( i==3 .and. j==4) f(i,j,k) = -18.0
          if( i==4 .and. j==4) f(i,j,k) = 19.0

          if( i==1 .and. j==5) f(i,j,k) = 20.0
          if( i==2 .and. j==5) f(i,j,k) = 21.0

          enddo
        enddo
    enddo
end subroutine

subroutine set_eb_flag(lo, hi, flag) bind(C, name = "set_eb_flag")

    implicit none

    integer, intent(in) :: lo(3), hi(3)

    integer, intent(inout) :: flag(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))

    integer          :: i, j, k


    do k = lo(3)+1, hi(3)-1
       do j = lo(2)+1, hi(2)-1
          do i = lo(1), hi(1)
          flag(i,j,k) = 0

          if( i==5 .and. j==1) flag(i,j,k) = 1

          if( i==2 .and. j==2) flag(i,j,k) = 2
          if( i==5 .and. j==2) flag(i,j,k) = 1

          if( i==2 .and. j==3) flag(i,j,k) = 2
          if( i==4 .and. j==3) flag(i,j,k) = 1
          if( i==5 .and. j==3) flag(i,j,k) = 1

          if( i==2 .and. j==4) flag(i,j,k) = 1
          if( i==3 .and. j==4) flag(i,j,k) = 2
          if( i==4 .and. j==4) flag(i,j,k) = 1

          if( i==2 .and. j==5) flag(i,j,k) = 1
          if( i==3 .and. j==5) flag(i,j,k) = -1
          if( i==4 .and. j==5) flag(i,j,k) = -1
          if( i==5 .and. j==5) flag(i,j,k) = -1

          enddo
        enddo
    enddo
end subroutine
