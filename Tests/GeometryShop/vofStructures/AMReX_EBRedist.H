#ifndef _EBRedist_H
#define _EBRedist_H

#include <AMReX_REAL.H>

#ifndef BL_LANG_FORT

typedef std::array<int,3>  intDIM;

struct NbrSten
{
    amrex::Real val[3][3][3];
    intDIM iv;
};

struct FaceSten
{
#if BL_SPACEDIM == 2
    amrex::Real val[3];
#else
    amrex::Real val[3][3];
#endif
    intDIM iv;
};

#endif

#endif
