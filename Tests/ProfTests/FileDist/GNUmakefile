AMREX_HOME = ../../..

DEBUG    = FALSE
#DEBUG    = TRUE

USE_MPI  = FALSE
#USE_MPI  = TRUE

USE_OMP  = TRUE
USE_OMP  = FALSE

#PROFILE        = TRUE
#COMM_PROFILE   = TRUE
#TRACE_PROFILE  = TRUE
PROFILE        = FALSE
COMM_PROFILE   = FALSE
TRACE_PROFILE  = FALSE


#USE_PROFPARSER = TRUE

COMP      = gnu
DIM       = 2 
PRECISION = DOUBLE

EBASE     = FileDist 

CEXE_sources += FileDist.cpp

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include $(AMREX_HOME)/Src/Base/Make.package
#include $(AMREX_HOME)/Src/Extern/ProfParser/Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.rules

