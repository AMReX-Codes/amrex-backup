AMREX_HOME ?= ../../../..
KOKKOS_INSTALL ?= $(AMREX_HOME)/Tools/kokkoslib/install


DEBUG	= FALSE
#DEBUG	= TRUE

DIM	= 3

COMP    = gnu 

USE_MPI   = TRUE 
USE_OMP   = FALSE
USE_CUDA  = TRUE
TINY_PROFILE = TRUE

KOKKOS_DEVICES=Cuda
KOKKOS_CUDA_OPTIONS=enable_lambda
KOKKOS_ARCH = "Volta70"

INCLUDE_LOCATIONS += $(KOKKOS_INSTALL)/include
LIBRARY_LOCATIONS += $(KOKKOS_INSTALL)/lib
LIBRARIES += -lkokkos

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ./Make.package
include $(AMREX_HOME)/Src/Base/Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
