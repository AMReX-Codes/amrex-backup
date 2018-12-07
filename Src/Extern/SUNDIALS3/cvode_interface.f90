module cvode_interface

  use fcvode_mod
  use fsunmat_dense_mod
  use fsunlinsol_dense_mod
  use fsunlinsol_spgmr_mod

  contains

  integer(c_int) function FCVDense(cvode_mem, N) result(ierr)
    use, intrinsic :: iso_c_binding
    use fnvector_serial_mod
    implicit none 
    type(c_ptr),     value :: cvode_mem
    integer(c_long), value :: N
    type(c_ptr)            :: sunmat_A
    type(c_ptr)            :: sunlinsol_LS
    type(c_ptr)            :: sunvec_y

    sunvec_y = FN_VNewEmpty_Serial(N)
    sunmat_A = FSUNDenseMatrix(N, N)
    sunlinsol_LS = FSUNDenseLinearSolver(sunvec_y, sunmat_A)
    ierr = FCVDlsSetLinearSolver(cvode_mem, sunlinsol_LS, sunmat_A)
  end function FCVDense

  integer(c_int) function FCVIter(cvode_mem, N, nKryBasis) 
    use, intrinsic :: iso_c_binding
    use fnvector_serial_mod
    implicit none 
    type(c_ptr),     value :: cvode_mem
    integer(c_long), value :: N
    integer(c_int), value  :: nKryBasis
    type(c_ptr)            :: sunlinsol_LS
    type(c_ptr)            :: sunvec_y
    integer(c_int)         :: ierr

    sunvec_y = FN_VNewEmpty_Serial(N)
    sunlinsol_LS = FSUNSPGMR(sunvec_y, PREC_NONE, nKryBasis)
    ierr = FCVSpilsSetLinearSolver(cvode_mem, sunlinsol_LS)
  end function FCVIter


end module cvode_interface
