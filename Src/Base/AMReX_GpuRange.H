#ifndef AMREX_GPU_RANGE_H_
#define AMREX_GPU_RANGE_H_

#include <AMReX_GpuQualifiers.H>
#include <AMReX_GpuControl.H>
#include <AMReX_GpuError.H>
#include <AMReX_Box.H>
#include <AMReX_TypeTraits.H>

namespace amrex {
namespace Gpu {

namespace range_detail {

//! integer version
template <typename T, typename = amrex::EnableIf_t<std::is_integral<T>::value> >
AMREX_GPU_HOST_DEVICE
long size (T const& b) noexcept { return static_cast<long>(b); }

template <typename T, typename = amrex::EnableIf_t<std::is_integral<T>::value> >
AMREX_GPU_HOST_DEVICE
long at (T const& b, long offset) noexcept { return offset; }

//! Box version
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE long size (Box const& b) noexcept
{
#ifdef AMREX_DEVICE_COMPILE
    return b.numPts();
#else
    return 1;
#endif
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE Box at (Box const& b, long offset) noexcept
{
#ifdef AMREX_DEVICE_COMPILE
    auto len = b.lengthnd();
    auto lo = b.loVectnd();
#if (AMREX_SPACEDIM == 1)
    long i0 = offset;
#elif (AMREX_SPACEDIM == 2)
    long i1 = offset / len[0];
    long i0 = offset - i1*len[0];
#elif (AMREX_SPACEDIM == 3)
    long i2 = offset / (len[0]*len[1]);
    long i1 = (offset - i2*(len[0]*len[1])) / len[0];
    long i0 = (offset - i2*(len[0]*len[1])) - i1*len[0];
#elif (AMREX_SPACEDIM == 4)
    long stride2 = len[0]*len[1];
    long stride3 = stride2*len[2];
    long i3 = offset / stride3;
    long i2 = (offset - i3*stride3) / stride2;
    long i1 = (offset - i3*stride3 - i2*stride2) / len[0];
    long i0 =  offset - i3*stride3 - i2*stride2 - i1*len[0];
#elif (AMREX_SPACEDIM == 5)
    long stride2 = len[0]*len[1];
    long stride3 = stride2*len[2];
    long stride4 = stride3*len[3];
    long i4 = offset / stride4;
    long i3 = (offset - i4*stride4) / stride3;
    long i2 = (offset - i4*stride4 - i3*stride3) / stride2;
    long i1 = (offset - i4*stride4 - i3*stride3 - i2*stride2) / len[0];
    long i0 =  offset - i4*stride4 - i3*stride3 - i2*stride2 - i1*len[0];
#else
    long stride2 = len[0]*len[1];
    long stride3 = stride2*len[2];
    long stride4 = stride3*len[3];
    long stride5 = stride4*len[4];
    long i5 = offset / stride5;
    long i4 = (offset - i5*stride5) / stride4;
    long i3 = (offset - i5*stride5 - i4*stride4) / stride3;
    long i2 = (offset - i5*stride5 - i4*stride4 - i3*stride3) / stride2;
    long i1 = (offset - i5*stride5 - i4*stride4 - i3*stride3 - i2*stride2) / len[0];
    long i0 =  offset - i5*stride5 - i4*stride4 - i3*stride3 - i2*stride2 - i1*len[0];
#endif
    IntVect iv{AMREX_D6_DECL(static_cast<int>(i0) + lo[0],
                             static_cast<int>(i1) + lo[1],
                             static_cast<int>(i2) + lo[2],
                             static_cast<int>(i3) + lo[3],
                             static_cast<int>(i4) + lo[4],
                             static_cast<int>(i5) + lo[5])};
    return Box(iv,iv,b.type(),b.dimension());
#else
    return b;
#endif
}

template <typename T>
struct range_impl
{
    AMREX_GPU_HOST_DEVICE
    explicit range_impl (T const& b) noexcept : m_b(b), m_n(range_detail::size(b)) {}

    struct iterator
    {
        AMREX_GPU_HOST_DEVICE
        iterator (T const& b, long i, long s) noexcept : mi_b(b), mi_i(i), mi_s(s)  {}

        AMREX_GPU_HOST_DEVICE
        void operator++ () noexcept { mi_i += mi_s; }

        AMREX_GPU_HOST_DEVICE
        bool operator!= (iterator const& rhs) const noexcept { return mi_i < rhs.mi_i; }

        AMREX_GPU_HOST_DEVICE
        T operator* () const noexcept { return range_detail::at(mi_b,mi_i); }

    private:
        T const& mi_b;
        long mi_i;
        long mi_s;
    };

    AMREX_GPU_HOST_DEVICE
    iterator begin () const noexcept {
#ifdef AMREX_DEVICE_COMPILE
        return iterator(m_b, blockDim.x*blockIdx.x+threadIdx.x, blockDim.x*gridDim.x);
#else
        return iterator(m_b,0,1);
#endif
    }

    AMREX_GPU_HOST_DEVICE
    iterator end () const noexcept { return iterator(m_b,m_n,0); }

private:
    T m_b;
    long m_n;
};
}

template <typename T>
AMREX_GPU_HOST_DEVICE
range_detail::range_impl<T> Range(T const& b) noexcept { return range_detail::range_impl<T>(b); }

}}

#endif
