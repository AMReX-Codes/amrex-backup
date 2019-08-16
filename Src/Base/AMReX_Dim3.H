#ifndef AMREX_DIM3_H_
#define AMREX_DIM3_H_

#include <iostream>
#include <AMReX_REAL.H>
#include <AMReX_Array.H>

namespace amrex {

struct Dim3 { int x; int y; int z; };
struct XDim3 { Real x; Real y; Real z; };

inline std::ostream& operator<< (std::ostream& os, const Dim3& d) {
    os << '(' << d.x << ',' << d.y << ',' << d.z << ')';
    return os;
}

#if (AMREX_SPACEDIM == 1)
struct DimN { int i0; };
struct XDimN { Real x0; };
inline std::ostream& operator<< (std::ostream& os, const DimN& d) {
    os << '(' << d.i0 <<')';
    return os;
}
#elif (AMREX_SPACEDIM == 2)
struct DimN { int i0; int i1; };
struct XDimN { Real x0; Real x1; };
inline std::ostream& operator<< (std::ostream& os, const DimN& d) {
    os << '(' << d.i0 << ',' << d.i1 <<')';
    return os;
}
#elif (AMREX_SPACEDIM == 3)
struct DimN { int i0; int i1; int i2; };
struct XDimN { Real x0; Real x1; Real x2; };
inline std::ostream& operator<< (std::ostream& os, const DimN& d) {
    os << '(' << d.i0 << ',' << d.i1 << ',' << d.i2 <<')';
    return os;
}
#elif (AMREX_SPACEDIM == 4)
struct DimN { int i0; int i1; int i2; int i3; };
struct XDimN { Real x0; Real x1; Real x2; Real x3; };
inline std::ostream& operator<< (std::ostream& os, const DimN& d) {
    os << '(' << d.i0 << ',' << d.i1 << ',' << d.i2 << ',' << d.i3 <<')';
    return os;
}
#elif (AMREX_SPACEDIM == 5)
struct DimN { int i0; int i1; int i2; int i3; int i4; int i5; };
struct XDimN { Real x0; Real x1; Real x2; Real x3; Real x4; Real x5; };
inline std::ostream& operator<< (std::ostream& os, const DimN& d) {
    os << '(' << d.i0 << ',' << d.i1 << ',' << d.i2
       << ',' << d.i3 << ',' << d.i4 <<')';
    return os;
}
#elif (AMREX_SPACEDIM == 6)
struct DimN { int i0; int i1; int i2; int i3; int i4; int i5; };
struct XDimN { Real x0; Real x1; Real x2; Real x3; Real x4; Real x5; };
inline std::ostream& operator<< (std::ostream& os, const DimN& d) {
    os << '(' << d.i0 << ',' << d.i1 << ',' << d.i2
       << ',' << d.i3 << ',' << d.i4 << ',' << d.i5 <<')';
    return os;
}
#endif

inline XDim3 makeXDim3 (const Array<Real,AMREX_SPACEDIM>& a) noexcept
{
#if (AMREX_SPACEDIM == 1)
    return XDim3{a[0], 0., 0.};
#elif (AMREX_SPACEDIM == 2)
    return XDim3{a[0], a[1], 0.};
#else
    return XDim3{a[0], a[1], a[2]};
#endif
}

}

#endif
