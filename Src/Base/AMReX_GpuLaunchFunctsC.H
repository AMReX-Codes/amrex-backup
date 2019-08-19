#ifndef AMREX_GPU_LAUNCH_FUNCTS_C_H_
#define AMREX_GPU_LAUNCH_FUNCTS_C_H_

namespace amrex {

template<typename T, typename L>
AMREX_FORCE_INLINE
void launch (T const& n, L&& f, std::size_t shared_mem_bytes=0) noexcept
{
    f(n);
}

template <typename T, typename L, typename M=amrex::EnableIf_t<std::is_integral<T>::value> >
AMREX_FORCE_INLINE
void For (T n, L&& f, std::size_t shared_mem_bytes=0) noexcept
{
    for (T i = 0; i < n; ++i) {
        f(i);
    }
}

template <typename T, typename L, typename M=amrex::EnableIf_t<std::is_integral<T>::value> >
AMREX_FORCE_INLINE
void ParallelFor (T n, L&& f, std::size_t shared_mem_bytes=0) noexcept
{
    AMREX_PRAGMA_SIMD
    for (T i = 0; i < n; ++i) {
        f(i);
    }
}

template <typename L>
AMREX_FORCE_INLINE
auto For (Box const& box, L&& f, std::size_t shared_mem_bytes=0) noexcept
    -> decltype(f(0,0,0))
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    for (int k = lo.z; k <= hi.z; ++k) {
    for (int j = lo.y; j <= hi.y; ++j) {
    for (int i = lo.x; i <= hi.x; ++i) {
        f(i,j,k);
    }}}
}

template <typename L>
AMREX_FORCE_INLINE
auto ParallelFor (Box const& box, L&& f, std::size_t shared_mem_bytes=0) noexcept
    -> decltype(f(0,0,0))
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    for (int k = lo.z; k <= hi.z; ++k) {
    for (int j = lo.y; j <= hi.y; ++j) {
    AMREX_PRAGMA_SIMD
    for (int i = lo.x; i <= hi.x; ++i) {
        f(i,j,k);
    }}}
}

template <typename T, typename L, typename M=amrex::EnableIf_t<std::is_integral<T>::value> >
AMREX_FORCE_INLINE
auto For (Box const& box, T ncomp, L&& f, std::size_t shared_mem_bytes=0) noexcept
    -> decltype(f(0,0,0,0))
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    for (T n = 0; n < ncomp; ++n) {
        for (int k = lo.z; k <= hi.z; ++k) {
        for (int j = lo.y; j <= hi.y; ++j) {
        for (int i = lo.x; i <= hi.x; ++i) {
            f(i,j,k,n);
        }}}
    }
}

template <typename T, typename L, typename M=amrex::EnableIf_t<std::is_integral<T>::value> >
AMREX_FORCE_INLINE
auto ParallelFor (Box const& box, T ncomp, L&& f, std::size_t shared_mem_bytes=0) noexcept
    -> decltype(f(0,0,0,0))
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    for (T n = 0; n < ncomp; ++n) {
        for (int k = lo.z; k <= hi.z; ++k) {
        for (int j = lo.y; j <= hi.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            f(i,j,k,n);
        }}}
    }
}

template <typename N, typename T, typename L1, typename L2,
          typename M=amrex::EnableIf_t<std::is_integral<N>::value> >
void FabReduce (Box const& box, N ncomp, T const& init_val,
                L1&& f1, L2&& f2, std::size_t shared_mem_bytes=0) noexcept
{
    auto r = init_val;
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    for (N n = 0; n < ncomp; ++n) {
        for (int k = lo.z; k <= hi.z; ++k) {
        for (int j = lo.y; j <= hi.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            f1(i,j,k,n,&r);
        }}}
    }
    f2(r);
}

template <typename T, typename L1, typename L2>
void FabReduce (Box const& box, T const& init_val,
                L1&& f1, L2&& f2, std::size_t shared_mem_bytes=0) noexcept
{
    auto r = init_val;
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    for (int k = lo.z; k <= hi.z; ++k) {
    for (int j = lo.y; j <= hi.y; ++j) {
    AMREX_PRAGMA_SIMD
    for (int i = lo.x; i <= hi.x; ++i) {
        f1(i,j,k,&r);
    }}}
    f2(r);
}

template <typename N, typename T, typename L1, typename L2,
          typename M=amrex::EnableIf_t<std::is_integral<N>::value> >
void VecReduce (N n, T const& init_val,
                L1&& f1, L2&& f2, std::size_t shared_mem_bytes=0) noexcept
{
    auto r = init_val;
    AMREX_PRAGMA_SIMD
    for (N i = 0; i < n; ++i) {
        f1(i,&r);
    }
    f2(r);
}

//
// IntVect versions
//

template <typename L>
AMREX_FORCE_INLINE
auto For (Box const& box, L&& f, std::size_t shared_mem_bytes=0) noexcept
    -> decltype(f(IntVect()))
{
    const auto lo = amrex::lboundn(box);
    const auto hi = amrex::uboundn(box);
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

template <typename L>
AMREX_FORCE_INLINE
auto ParallelFor (Box const& box, L&& f, std::size_t shared_mem_bytes=0) noexcept
    -> decltype(f(IntVect()))
{
    const auto lo = amrex::lboundn(box);
    const auto hi = amrex::uboundn(box);
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

template <typename T, typename L, typename M=amrex::EnableIf_t<std::is_integral<T>::value> >
AMREX_FORCE_INLINE
auto For (Box const& box, T ncomp, L&& f, std::size_t shared_mem_bytes=0) noexcept
    -> decltype(f(IntVect(),0))
{
    const auto lo = amrex::lboundn(box);
    const auto hi = amrex::uboundn(box);
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

template <typename T, typename L, typename M=amrex::EnableIf_t<std::is_integral<T>::value> >
AMREX_FORCE_INLINE
auto ParallelFor (Box const& box, T ncomp, L&& f, std::size_t shared_mem_bytes=0) noexcept
    -> decltype(f(IntVect(),0))
{
    const auto lo = amrex::lboundn(box);
    const auto hi = amrex::uboundn(box);
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
