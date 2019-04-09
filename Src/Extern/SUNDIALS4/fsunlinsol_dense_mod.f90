! This file was automatically generated by SWIG (http://www.swig.org).
! Version 4.0.0
!
! Do not make changes to this file unless you know what you are doing--modify
! the SWIG interface file instead.

! ---------------------------------------------------------------
! Programmer(s): Auto-generated by swig.
! ---------------------------------------------------------------
! SUNDIALS Copyright Start
! Copyright (c) 2002-2019, Lawrence Livermore National Security
! and Southern Methodist University.
! All rights reserved.
!
! See the top-level LICENSE and NOTICE files for details.
!
! SPDX-License-Identifier: BSD-3-Clause
! SUNDIALS Copyright End
! ---------------------------------------------------------------

module fsunlinsol_dense_mod
 use, intrinsic :: ISO_C_BINDING
 use fsunlinsol_mod
 use fsundials_types_mod
 use fnvector_mod
 use fsundials_types_mod
 implicit none
 private

 ! DECLARATION CONSTRUCTS
 public :: FSUNLinSol_Dense
 public :: FSUNDenseLinearSolver
 public :: FSUNLinSolGetType_Dense
 public :: FSUNLinSolInitialize_Dense
 public :: FSUNLinSolSetup_Dense
 public :: FSUNLinSolSolve_Dense
 public :: FSUNLinSolLastFlag_Dense
 public :: FSUNLinSolSpace_Dense
 public :: FSUNLinSolFree_Dense

! WRAPPER DECLARATIONS
interface
function FSUNLinSol_Dense(y, a) &
bind(C, name="SUNLinSol_Dense") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: y
type(C_PTR), value :: a
type(C_PTR) :: fresult
end function

function FSUNDenseLinearSolver(y, a) &
bind(C, name="SUNDenseLinearSolver") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: y
type(C_PTR), value :: a
type(C_PTR) :: fresult
end function

function FSUNLinSolGetType_Dense(s) &
bind(C, name="SUNLinSolGetType_Dense") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: s
integer(C_INT) :: fresult
end function

function FSUNLinSolInitialize_Dense(s) &
bind(C, name="SUNLinSolInitialize_Dense") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: s
integer(C_INT) :: fresult
end function

function FSUNLinSolSetup_Dense(s, a) &
bind(C, name="SUNLinSolSetup_Dense") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: s
type(C_PTR), value :: a
integer(C_INT) :: fresult
end function

function FSUNLinSolSolve_Dense(s, a, x, b, tol) &
bind(C, name="SUNLinSolSolve_Dense") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: s
type(C_PTR), value :: a
type(C_PTR), value :: x
type(C_PTR), value :: b
real(C_DOUBLE), value :: tol
integer(C_INT) :: fresult
end function

function FSUNLinSolLastFlag_Dense(s) &
bind(C, name="SUNLinSolLastFlag_Dense") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: s
integer(C_LONG) :: fresult
end function

function FSUNLinSolSpace_Dense(s, lenrwls, leniwls) &
bind(C, name="SUNLinSolSpace_Dense") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: s
integer(C_LONG) :: lenrwls
integer(C_LONG) :: leniwls
integer(C_INT) :: fresult
end function

function FSUNLinSolFree_Dense(s) &
bind(C, name="SUNLinSolFree_Dense") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: s
integer(C_INT) :: fresult
end function

end interface


end module