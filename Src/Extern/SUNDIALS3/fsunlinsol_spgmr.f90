! -----------------------------------------------------------------
! Programmer(s): Cody J. Balos @ LLNL
! -----------------------------------------------------------------
! LLNS Copyright Start
! Copyright (c) 2017, Lawrence Livermore National Security
! This work was performed under the auspices of the U.S. Department
! of Energy by Lawrence Livermore National Laboratory in part under
! Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
! Produced at the Lawrence Livermore National Laboratory.
! All rights reserved.
! For details, see the LICENSE file.
! LLNS Copyright End
! -----------------------------------------------------------------
! This file contains a Fortran module for interfacing directly with
! the SUNDIALS sparse matrix using the ISO_C_BINDING module.
! -----------------------------------------------------------------

module fsunlinsol_spgmr_mod

  use, intrinsic :: iso_c_binding, only : c_int
  
  integer(c_int), parameter :: PREC_LEFT    = 1 ! COLAMD
  integer(c_int), parameter :: PREC_RIGHT   = 2
  integer(c_int), parameter :: PREC_NONE    = 0

  !======= Interfaces ========
  interface

    ! =================================================================
    ! Constructors
    ! =================================================================

    type(c_ptr) function FSUNSPGMR(y, pretype, maxl) &
        bind(C,name='SUNSPGMR')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), value    :: y
      integer(c_int), value :: pretype
      integer(c_int), value :: maxl
    end function FSUNSPGMR

    ! =================================================================
    ! Destructors
    ! =================================================================

    subroutine FSUNLinSolFree_SPGMR(LS) &
        bind(C,name='SUNLinSolFree_SPGMR')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), value :: LS
    end subroutine FSUNLinSolFree_SPGMR

    ! =================================================================
    ! Setter/init routines
    ! =================================================================

    !integer(c_int) function FSUNLinSolSetPreconditioner_SPGMR(LS, pdata, pset, psolve) &
    !  bind(C,name='SUNLinSolSetPreconditioner_SPGMR')
    !  use, intrinsic :: iso_c_binding
    !  implicit none
    !  type(c_funptr), value :: pset
    !  type(c_funptr), value :: psolve
    !end function FSUNLinSolSetPreconditioner_SPGMR

    ! =================================================================
    ! Operations
    ! =================================================================

    ! -----------------------------------------------------------------
    ! NOT INTERFACED SUNLinSolGetType_SPGMR
    ! -----------------------------------------------------------------

    integer(c_int) function FSUNLinSolGetType_SPGMR(LS) &
        bind(C,name='SUNLinSolGetType_SPGMR')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), value :: LS
    end function FSUNLinSolGetType_SPGMR

    integer(c_int) function FSUNLinSolInitialize_SPGMR(LS) &
        bind(C,name='SUNLinSolInitialize_SPGMR')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), value :: LS
    end function FSUNLinSolInitialize_SPGMR

    integer(c_int) function FSUNLinSolSetup_SPGMR(LS, A) &
        bind(C,name='SUNLinSolSetup_SPGMR')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), value :: LS
      type(c_ptr), value :: A
    end function FSUNLinSolSetup_SPGMR

    integer(c_int) function FSUNLinSolSolve_SPGMR(LS, A, x, b, tol) &
        bind(C,name='SUNLinSolSolve_SPGMR')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), value    :: LS
      type(c_ptr), value    :: A
      type(c_ptr), value    :: x
      type(c_ptr), value    :: b
      real(c_double), value :: tol
    end function FSUNLinSolSolve_SPGMR

    integer(c_long) function FSUNLinSolLastFlag_SPGMR(LS) &
        bind(C,name='SUNLinSolLastFlag_SPGMR')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), value :: LS
    end function FSUNLinSolLastFlag_SPGMR

    integer(c_int) function FSUNLinSolSpace_SPGMR(LS, lenrwLS, leniwLS) &
        bind(C,name='SUNLinSolSpace_SPGMR')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), value :: LS
      integer(c_long), value :: lenrwLS
      integer(c_long), value :: leniwLS
    end function FSUNLinSolSpace_SPGMR

  end interface

end module fsunlinsol_spgmr_mod
