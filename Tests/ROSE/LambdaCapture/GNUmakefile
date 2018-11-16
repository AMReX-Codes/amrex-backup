DEBUG = FALSE
DIM = 3
USE_MPI = FALSE
USE_OMP = FALSE

USE_CUDA = TRUE
COMP = gnu

AMREX_HOME ?= ../../..
include $(AMREX_HOME)/Tools/GNUMake/Make.defs

CEXE_sources += main.cpp
include $(AMREX_HOME)/Src/Base/Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
