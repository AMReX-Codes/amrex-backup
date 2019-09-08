#ifndef AMREX_BASEFAB_UTILITY_H_
#define AMREX_BASEFAB_UTILITY_H_

#include <AMReX_BaseFab.H>

namespace amrex {

template <class Tto, class Tfrom>
AMREX_GPU_HOST_DEVICE
void
cast (BaseFab<Tto>& tofab, BaseFab<Tfrom> const& fromfab,
      Box const& bx, SrcComp scomp, DestComp dcomp, NumComps ncomp) noexcept
{
    ArrayN<Tto> const& d = amrex::shift(tofab.arraynd(), dcomp.i);
    ArrayN<Tfrom const> const& s = amrex::shift(fromfab.const_arraynd(), scomp.i);
    amrex::LoopConcurrent(bx, ncomp.n, [=] (IntVect const& iv, int n) noexcept
    {
        d(iv,n) = static_cast<Tto>(s(iv,n));
    });
}

}

#endif
