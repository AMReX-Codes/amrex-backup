# AMREX_HOME defines the directory in which we will find all the AMReX code.
# If you set AMREX_HOME as an environment variable, this line will be ignored
AMREX_HOME ?= ../../../../../
KOKKOS_INSTALL ?= $(AMREX_HOME)/Tools/kokkoslib/install

DEBUG     = FALSE
USE_MPI   = TRUE 
USE_OMP   = FALSE
USE_CUDA  = TRUE 
#TINY_PROFILE = TRUE
COMP      = gcc
DIM       = 3 

Bpack   := ./Make.package
Blocs   := .

KOKKOS_DEVICES=Cuda
KOKKOS_CUDA_OPTIONS=enable_lambda
KOKKOS_ARCH = "Volta70"

#CXXFLAGS += $(KOKKOS_CPP_DEPENDS)
#LINKFLAGS += $(KOKKOS_LINK_DEPENDS) $(KOKKOS_LDFLAGS)

CPPFLAGS += -DAMREX_KOKKOS
INCLUDE_LOCATIONS += $(KOKKOS_INSTALL)/include
LIBRARY_LOCATIONS += $(KOKKOS_INSTALL)/lib 
LIBRARIES += -lkokkos

#ifneq (,$(findstring Cuda,$(KOKKOS_DEVICES)))
#CXX = ${KOKKOS_PATH}/bin/nvcc_wrapper
#else
#CXX = g++
#endif

#include $(KOKKOS_PATH)/Makefile.kokkos
include ../../Make.HEAT


#$(EXE): $(OBJ) $(KOKKOS_LINK_DEPENDS)
#	$(LINK) $(KOKKOS_LDFLAGS) $(LINKFLAGS) $(EXTRA_PATH) $(OBJ) $(KOKKOS_LIBS) $(LIB) -o $(EXE)

# Compilation rules

#%.o:%.cpp $(KOKKOS_CPP_DEPENDS) $(HEADERS)

