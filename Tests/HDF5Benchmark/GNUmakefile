AMREX_HOME ?= ../../

DEBUG	= TRUE
DEBUG	= FALSE

DIM	= 3

COMP    = gnu 

PRECISION = DOUBLE

USE_MPI   = TRUE
USE_OMP   = FALSE

TINY_PROFILE = TRUE

###################################################

EBASE     = main

USE_HDF5 = TRUE

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ./Make.package
include $(AMREX_HOME)/Src/Base/Make.package
include $(AMREX_HOME)/Src/Particle/Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.rules

