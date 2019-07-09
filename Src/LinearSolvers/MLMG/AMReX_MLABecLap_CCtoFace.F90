subroutine cc_to_face(lo, hi, cc_dat, cc_lo, cc_hi, face_dat, face_lo,face_hi, dir) bind(C, name="cc_to_face")
  !  move cell centered data to be face centered
  use amrex_fort_module, only : amrex_real

  implicit none

  integer, intent(in) :: lo(3), hi(3), cc_lo(3), cc_hi(3), face_lo(3), face_hi(3), dir
  real(amrex_real), intent(inout) :: cc_dat(cc_lo(1):cc_hi(1),cc_lo(2):cc_hi(2),cc_lo(3):cc_hi(3))
  real(amrex_real), intent(inout) :: face_dat(face_lo(1):face_hi(1),face_lo(2):face_hi(2),face_lo(3):face_hi(3))

  integer          :: i,j, k


!print *,'lo, hi, cc_lo, cc_hi, face_lo, face_hi=',lo, hi, cc_lo, cc_hi, face_lo, face_hi

!print *,'cc_dat(0,0,0)=',cc_dat(0,0,0)


! We condition based on direction
if (dir .EQ. 0) then
  do i = lo(1),hi(1) + 1
    do j = lo(2),hi(2)
      do k = lo(3),hi(3)
	face_dat(i,j,k) = (-1.d0*cc_dat(i-2,j,k) + 7.d0*cc_dat(i-1,j,k) & 
                             + 7.d0*cc_dat(i,j,k) - 1.d0*cc_dat(i+1,j,k))/(12.d0)

      end do
    end do
  end do
elseif(dir .EQ. 1) then

  do i = lo(1),hi(1)
    do j = lo(2),hi(2) + 1
      do k = lo(3),hi(3)

	face_dat(i,j,k) = (-1.d0*cc_dat(i,j-2,k) + 7.d0*cc_dat(i,j-1,k) & 
                             + 7.d0*cc_dat(i,j,k) - 1.d0*cc_dat(i,j+1,k))/(12.d0)

      end do
    end do
  end do

else
  do i = lo(1),hi(1)
    do j = lo(2),hi(2)
      do k = lo(3),hi(3) + 1

	face_dat(i,j,k) = (-1.d0*cc_dat(i,j,k-2) + 7.d0*cc_dat(i,j,k-1) & 
                             + 7.d0*cc_dat(i,j,k) - 1.d0*cc_dat(i,j,k+1))/(12.d0)

      end do
    end do
  end do


endif






end subroutine cc_to_face
