#ifndef AMREX_LOOP_H_
#define AMREX_LOOP_H_

#include <AMReX_Box.H>

namespace amrex {

template <class F>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void Loop (Dim3 lo, Dim3 hi, F&& f) noexcept
{
    for (int k = lo.z; k <= hi.z; ++k) {
    for (int j = lo.y; j <= hi.y; ++j) {
    for (int i = lo.x; i <= hi.x; ++i) {
        f(i,j,k);
    }}}
}

template <class F>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void Loop (Dim3 lo, Dim3 hi, int ncomp, F&& f) noexcept
{
    for (int n = 0; n < ncomp; ++n) {
    for (int k = lo.z; k <= hi.z; ++k) {
    for (int j = lo.y; j <= hi.y; ++j) {
    for (int i = lo.x; i <= hi.x; ++i) {
        f(i,j,k,n);
    }}}}
}

template <class F>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void LoopConcurrent (Dim3 lo, Dim3 hi, F&& f) noexcept
{
    for (int k = lo.z; k <= hi.z; ++k) {
    for (int j = lo.y; j <= hi.y; ++j) {
    AMREX_PRAGMA_SIMD
    for (int i = lo.x; i <= hi.x; ++i) {
        f(i,j,k);
    }}}
}

template <class F>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void LoopConcurrent (Dim3 lo, Dim3 hi, int ncomp, F&& f) noexcept
{
    for (int n = 0; n < ncomp; ++n) {
    for (int k = lo.z; k <= hi.z; ++k) {
    for (int j = lo.y; j <= hi.y; ++j) {
    AMREX_PRAGMA_SIMD
    for (int i = lo.x; i <= hi.x; ++i) {
        f(i,j,k,n);
    }}}}
}

template <class F>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
auto Loop (Box const& bx, F&& f) noexcept
    -> decltype(f(0,0,0))
{
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);
    for (int k = lo.z; k <= hi.z; ++k) {
    for (int j = lo.y; j <= hi.y; ++j) {
    for (int i = lo.x; i <= hi.x; ++i) {
        f(i,j,k);
    }}}
}

template <class F>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
auto Loop (Box const& bx, int ncomp, F&& f) noexcept
    -> decltype(f(0,0,0,0))
{
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);
    for (int n = 0; n < ncomp; ++n) {
    for (int k = lo.z; k <= hi.z; ++k) {
    for (int j = lo.y; j <= hi.y; ++j) {
    for (int i = lo.x; i <= hi.x; ++i) {
        f(i,j,k,n);
    }}}}
}

template <class F>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
auto LoopConcurrent (Box const& bx, F&& f) noexcept
    -> decltype(f(0,0,0))
{
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);
    for (int k = lo.z; k <= hi.z; ++k) {
    for (int j = lo.y; j <= hi.y; ++j) {
    AMREX_PRAGMA_SIMD
    for (int i = lo.x; i <= hi.x; ++i) {
        f(i,j,k);
    }}}
}

template <class F>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
auto LoopConcurrent (Box const& bx, int ncomp, F&& f) noexcept
    ->decltype(f(0,0,0,0))
{
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);
    for (int n = 0; n < ncomp; ++n) {
    for (int k = lo.z; k <= hi.z; ++k) {
    for (int j = lo.y; j <= hi.y; ++j) {
    AMREX_PRAGMA_SIMD
    for (int i = lo.x; i <= hi.x; ++i) {
        f(i,j,k,n);
    }}}}
}

template <class F>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
auto Loop (Box const& bx, F&& f) noexcept
    -> decltype(f(IntVect()))
{
    const auto lo = amrex::lboundn(bx);
    const auto hi = amrex::uboundn(bx);
#if (AMREX_SPACEDIM == 6)
    for (int i5 = lo.i5; i5 <= hi.i5; ++i5) {
#endif
#if (AMREX_SPACEDIM >= 5)
    for (int i4 = lo.i4; i4 <= hi.i4; ++i4) {
#endif
#if (AMREX_SPACEDIM >= 4)
    for (int i3 = lo.i3; i3 <= hi.i3; ++i3) {
#endif
#if (AMREX_SPACEDIM >= 3)
    for (int i2 = lo.i2; i2 <= hi.i2; ++i2) {
#endif
#if (AMREX_SPACEDIM >= 2)
    for (int i1 = lo.i1; i1 <= hi.i1; ++i1) {
#endif
    for (int i0 = lo.i0; i0 <= hi.i0; ++i0) {
        f(IntVect(AMREX_D6_DECL(i0,i1,i2,i3,i4,i5)));
    AMREX_D6_TERM(},},},},},})
}

template <class F>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
auto Loop (Box const& bx, int ncomp, F&& f) noexcept
    -> decltype(f(IntVect(),0))
{
    const auto lo = amrex::lboundn(bx);
    const auto hi = amrex::uboundn(bx);
    for (int n = 0; n < ncomp; ++n) {
#if (AMREX_SPACEDIM == 6)
    for (int i5 = lo.i5; i5 <= hi.i5; ++i5) {
#endif
#if (AMREX_SPACEDIM >= 5)
    for (int i4 = lo.i4; i4 <= hi.i4; ++i4) {
#endif
#if (AMREX_SPACEDIM >= 4)
    for (int i3 = lo.i3; i3 <= hi.i3; ++i3) {
#endif
#if (AMREX_SPACEDIM >= 3)
    for (int i2 = lo.i2; i2 <= hi.i2; ++i2) {
#endif
#if (AMREX_SPACEDIM >= 2)
    for (int i1 = lo.i1; i1 <= hi.i1; ++i1) {
#endif
    for (int i0 = lo.i0; i0 <= hi.i0; ++i0) {
        f(IntVect(AMREX_D6_DECL(i0,i1,i2,i3,i4,i5)),n);
    AMREX_D6_TERM(},},},},},})
    }
}

template <class F>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
auto LoopConcurrent (Box const& bx, F&& f) noexcept
    -> decltype(f(IntVect()))
{
    const auto lo = amrex::lboundn(bx);
    const auto hi = amrex::uboundn(bx);
#if (AMREX_SPACEDIM == 6)
    for (int i5 = lo.i5; i5 <= hi.i5; ++i5) {
#endif
#if (AMREX_SPACEDIM >= 5)
    for (int i4 = lo.i4; i4 <= hi.i4; ++i4) {
#endif
#if (AMREX_SPACEDIM >= 4)
    for (int i3 = lo.i3; i3 <= hi.i3; ++i3) {
#endif
#if (AMREX_SPACEDIM >= 3)
    for (int i2 = lo.i2; i2 <= hi.i2; ++i2) {
#endif
#if (AMREX_SPACEDIM >= 2)
    for (int i1 = lo.i1; i1 <= hi.i1; ++i1) {
#endif
    AMREX_PRAGMA_SIMD
    for (int i0 = lo.i0; i0 <= hi.i0; ++i0) {
        f(IntVect(AMREX_D6_DECL(i0,i1,i2,i3,i4,i5)));
    AMREX_D6_TERM(},},},},},})
}

template <class F>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
auto LoopConcurrent (Box const& bx, int ncomp, F&& f) noexcept
    ->decltype(f(IntVect(),0))
{
    const auto lo = amrex::lboundn(bx);
    const auto hi = amrex::uboundn(bx);
    for (int n = 0; n < ncomp; ++n) {
#if (AMREX_SPACEDIM == 6)
    for (int i5 = lo.i5; i5 <= hi.i5; ++i5) {
#endif
#if (AMREX_SPACEDIM >= 5)
    for (int i4 = lo.i4; i4 <= hi.i4; ++i4) {
#endif
#if (AMREX_SPACEDIM >= 4)
    for (int i3 = lo.i3; i3 <= hi.i3; ++i3) {
#endif
#if (AMREX_SPACEDIM >= 3)
    for (int i2 = lo.i2; i2 <= hi.i2; ++i2) {
#endif
#if (AMREX_SPACEDIM >= 2)
    for (int i1 = lo.i1; i1 <= hi.i1; ++i1) {
#endif
    AMREX_PRAGMA_SIMD
    for (int i0 = lo.i0; i0 <= hi.i0; ++i0) {
        f(IntVect(AMREX_D6_DECL(i0,i1,i2,i3,i4,i5)),n);
    AMREX_D6_TERM(},},},},},})
    }
}

// The functions above are __host__ __device__ functions.  If f is not a
// __host__ __device__ function, we will get warning about calling __host__
// function from a __host__ __device__ function.  This is ugly.  To get rid
// of the warning, we have to use the functions below for those situations.

template <class F>
AMREX_FORCE_INLINE
void LoopOnCpu (Dim3 lo, Dim3 hi, F&& f) noexcept
{
    for (int k = lo.z; k <= hi.z; ++k) {
    for (int j = lo.y; j <= hi.y; ++j) {
    for (int i = lo.x; i <= hi.x; ++i) {
        f(i,j,k);
    }}}
}

template <class F>
AMREX_FORCE_INLINE
void LoopOnCpu (Dim3 lo, Dim3 hi, int ncomp, F&& f) noexcept
{
    for (int n = 0; n < ncomp; ++n) {
    for (int k = lo.z; k <= hi.z; ++k) {
    for (int j = lo.y; j <= hi.y; ++j) {
    for (int i = lo.x; i <= hi.x; ++i) {
        f(i,j,k,n);
    }}}}
}

template <class F>
AMREX_FORCE_INLINE
void LoopConcurrentOnCpu (Dim3 lo, Dim3 hi, F&& f) noexcept
{
    for (int k = lo.z; k <= hi.z; ++k) {
    for (int j = lo.y; j <= hi.y; ++j) {
    AMREX_PRAGMA_SIMD
    for (int i = lo.x; i <= hi.x; ++i) {
        f(i,j,k);
    }}}
}

template <class F>
AMREX_FORCE_INLINE
void LoopConcurrentOnCpu (Dim3 lo, Dim3 hi, int ncomp, F&& f) noexcept
{
    for (int n = 0; n < ncomp; ++n) {
    for (int k = lo.z; k <= hi.z; ++k) {
    for (int j = lo.y; j <= hi.y; ++j) {
    AMREX_PRAGMA_SIMD
    for (int i = lo.x; i <= hi.x; ++i) {
        f(i,j,k,n);
    }}}}
}

template <class F>
AMREX_FORCE_INLINE
auto LoopOnCpu (Box const& bx, F&& f) noexcept
    -> decltype(f(0,0,0))
{
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);
    for (int k = lo.z; k <= hi.z; ++k) {
    for (int j = lo.y; j <= hi.y; ++j) {
    for (int i = lo.x; i <= hi.x; ++i) {
        f(i,j,k);
    }}}
}

template <class F>
AMREX_FORCE_INLINE
auto LoopOnCpu (Box const& bx, int ncomp, F&& f) noexcept
    -> decltype(f(0,0,0,0))
{
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);
    for (int n = 0; n < ncomp; ++n) {
    for (int k = lo.z; k <= hi.z; ++k) {
    for (int j = lo.y; j <= hi.y; ++j) {
    for (int i = lo.x; i <= hi.x; ++i) {
        f(i,j,k,n);
    }}}}
}

template <class F>
AMREX_FORCE_INLINE
auto LoopConcurrentOnCpu (Box const& bx, F&& f) noexcept
    -> decltype(f(0,0,0))
{
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);
    for (int k = lo.z; k <= hi.z; ++k) {
    for (int j = lo.y; j <= hi.y; ++j) {
    AMREX_PRAGMA_SIMD
    for (int i = lo.x; i <= hi.x; ++i) {
        f(i,j,k);
    }}}
}

template <class F>
AMREX_FORCE_INLINE
auto LoopConcurrentOnCpu (Box const& bx, int ncomp, F&& f) noexcept
    -> decltype(f(0,0,0,0))
{
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);
    for (int n = 0; n < ncomp; ++n) {
    for (int k = lo.z; k <= hi.z; ++k) {
    for (int j = lo.y; j <= hi.y; ++j) {
    AMREX_PRAGMA_SIMD
    for (int i = lo.x; i <= hi.x; ++i) {
        f(i,j,k,n);
    }}}}
}

template <class F>
AMREX_FORCE_INLINE
auto LoopOnCpu (Box const& bx, F&& f) noexcept
    -> decltype(f(IntVect()))
{
    const auto lo = amrex::lboundn(bx);
    const auto hi = amrex::uboundn(bx);
#if (AMREX_SPACEDIM == 6)
    for (int i5 = lo.i5; i5 <= hi.i5; ++i5) {
#endif
#if (AMREX_SPACEDIM >= 5)
    for (int i4 = lo.i4; i4 <= hi.i4; ++i4) {
#endif
#if (AMREX_SPACEDIM >= 4)
    for (int i3 = lo.i3; i3 <= hi.i3; ++i3) {
#endif
#if (AMREX_SPACEDIM >= 3)
    for (int i2 = lo.i2; i2 <= hi.i2; ++i2) {
#endif
#if (AMREX_SPACEDIM >= 2)
    for (int i1 = lo.i1; i1 <= hi.i1; ++i1) {
#endif
    for (int i0 = lo.i0; i0 <= hi.i0; ++i0) {
        f(IntVect(AMREX_D6_DECL(i0,i1,i2,i3,i4,i5)));
    AMREX_D6_TERM(},},},},},})
}

template <class F>
AMREX_FORCE_INLINE
auto LoopOnCpu (Box const& bx, int ncomp, F&& f) noexcept
    -> decltype(f(IntVect(),0))
{
    const auto lo = amrex::lboundn(bx);
    const auto hi = amrex::uboundn(bx);
    for (int n = 0; n < ncomp; ++n) {
#if (AMREX_SPACEDIM == 6)
    for (int i5 = lo.i5; i5 <= hi.i5; ++i5) {
#endif
#if (AMREX_SPACEDIM >= 5)
    for (int i4 = lo.i4; i4 <= hi.i4; ++i4) {
#endif
#if (AMREX_SPACEDIM >= 4)
    for (int i3 = lo.i3; i3 <= hi.i3; ++i3) {
#endif
#if (AMREX_SPACEDIM >= 3)
    for (int i2 = lo.i2; i2 <= hi.i2; ++i2) {
#endif
#if (AMREX_SPACEDIM >= 2)
    for (int i1 = lo.i1; i1 <= hi.i1; ++i1) {
#endif
    for (int i0 = lo.i0; i0 <= hi.i0; ++i0) {
        f(IntVect(AMREX_D6_DECL(i0,i1,i2,i3,i4,i5)),n);
    AMREX_D6_TERM(},},},},},})
    }
}

template <class F>
AMREX_FORCE_INLINE
auto LoopConcurrentOnCpu (Box const& bx, F&& f) noexcept
    -> decltype(f(IntVect()))
{
    const auto lo = amrex::lboundn(bx);
    const auto hi = amrex::uboundn(bx);
#if (AMREX_SPACEDIM == 6)
    for (int i5 = lo.i5; i5 <= hi.i5; ++i5) {
#endif
#if (AMREX_SPACEDIM >= 5)
    for (int i4 = lo.i4; i4 <= hi.i4; ++i4) {
#endif
#if (AMREX_SPACEDIM >= 4)
    for (int i3 = lo.i3; i3 <= hi.i3; ++i3) {
#endif
#if (AMREX_SPACEDIM >= 3)
    for (int i2 = lo.i2; i2 <= hi.i2; ++i2) {
#endif
#if (AMREX_SPACEDIM >= 2)
    for (int i1 = lo.i1; i1 <= hi.i1; ++i1) {
#endif
    AMREX_PRAGMA_SIMD
    for (int i0 = lo.i0; i0 <= hi.i0; ++i0) {
        f(IntVect(AMREX_D6_DECL(i0,i1,i2,i3,i4,i5)));
    AMREX_D6_TERM(},},},},},})
}

template <class F>
AMREX_FORCE_INLINE
auto LoopConcurrentOnCpu (Box const& bx, int ncomp, F&& f) noexcept
    -> decltype(f(IntVect(),0))
{
    const auto lo = amrex::lboundn(bx);
    const auto hi = amrex::uboundn(bx);
    for (int n = 0; n < ncomp; ++n) {
#if (AMREX_SPACEDIM == 6)
    for (int i5 = lo.i5; i5 <= hi.i5; ++i5) {
#endif
#if (AMREX_SPACEDIM >= 5)
    for (int i4 = lo.i4; i4 <= hi.i4; ++i4) {
#endif
#if (AMREX_SPACEDIM >= 4)
    for (int i3 = lo.i3; i3 <= hi.i3; ++i3) {
#endif
#if (AMREX_SPACEDIM >= 3)
    for (int i2 = lo.i2; i2 <= hi.i2; ++i2) {
#endif
#if (AMREX_SPACEDIM >= 2)
    for (int i1 = lo.i1; i1 <= hi.i1; ++i1) {
#endif
    AMREX_PRAGMA_SIMD
    for (int i0 = lo.i0; i0 <= hi.i0; ++i0) {
        f(IntVect(AMREX_D6_DECL(i0,i1,i2,i3,i4,i5)),n);
    AMREX_D6_TERM(},},},},},})
    }
}

}

#endif
