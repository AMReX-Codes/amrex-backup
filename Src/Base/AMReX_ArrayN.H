#ifndef AMREX_ARRAY_N_H_
#define AMREX_ARRAY_N_H_

#include <AMReX_IntVect.H>

namespace amrex {

    template <typename T>
    struct ArrayN
    {
        T* AMREX_RESTRICT p;
        LDimN stride;
        DimN begin;
        DimN end;
        int ncomp;

        AMREX_GPU_HOST_DEVICE
        constexpr ArrayN (T* a_p, DimN const& a_begin, DimN const& a_end, int a_ncomp) noexcept
            : p(a_p),
#if (AMREX_SPACEDIM == 1)
              stride{ (long)(a_end.i0-a_begin.i0) },
#elif (AMREX_SPACEDIM == 2)
              stride{ (long)(a_end.i0-a_begin.i0)
                     *(long)(a_end.i1-a_begin.i1) ,
                      (long)(a_end.i0-a_begin.i0) },
#elif (AMREX_SPACEDIM == 3)
              stride{ (long)(a_end.i0-a_begin.i0)
                     *(long)(a_end.i1-a_begin.i1)
                     *(long)(a_end.i2-a_begin.i2) ,
                      (long)(a_end.i0-a_begin.i0) ,
                      (long)(a_end.i0-a_begin.i0)
                     *(long)(a_end.i1-a_begin.i1) },
#elif (AMREX_SPACEDIM == 4)
              stride{ (long)(a_end.i0-a_begin.i0)
                     *(long)(a_end.i1-a_begin.i1)
                     *(long)(a_end.i2-a_begin.i2)
                     *(long)(a_end.i3-a_begin.i3) ,
                      (long)(a_end.i0-a_begin.i0) ,
                      (long)(a_end.i0-a_begin.i0)
                     *(long)(a_end.i1-a_begin.i1) ,
                      (long)(a_end.i0-a_begin.i0)
                     *(long)(a_end.i1-a_begin.i1)
                     *(long)(a_end.i2-a_begin.i2) },
#elif (AMREX_SPACEDIM == 5)
              stride{ (long)(a_end.i0-a_begin.i0)
                     *(long)(a_end.i1-a_begin.i1)
                     *(long)(a_end.i2-a_begin.i2)
                     *(long)(a_end.i3-a_begin.i3)
                     *(long)(a_end.i4-a_begin.i4) ,
                      (long)(a_end.i0-a_begin.i0) ,
                      (long)(a_end.i0-a_begin.i0)
                     *(long)(a_end.i1-a_begin.i1) ,
                      (long)(a_end.i0-a_begin.i0)
                     *(long)(a_end.i1-a_begin.i1)
                     *(long)(a_end.i2-a_begin.i2) ,
                      (long)(a_end.i0-a_begin.i0)
                     *(long)(a_end.i1-a_begin.i1)
                     *(long)(a_end.i2-a_begin.i2)
                     *(long)(a_end.i3-a_begin.i3) },
#else
              stride{ (long)(a_end.i0-a_begin.i0)
                     *(long)(a_end.i1-a_begin.i1)
                     *(long)(a_end.i2-a_begin.i2)
                     *(long)(a_end.i3-a_begin.i3)
                     *(long)(a_end.i4-a_begin.i4)
                     *(long)(a_end.i5-a_begin.i5),
                      (long)(a_end.i0-a_begin.i0) ,
                      (long)(a_end.i0-a_begin.i0)
                     *(long)(a_end.i1-a_begin.i1) ,
                      (long)(a_end.i0-a_begin.i0)
                     *(long)(a_end.i1-a_begin.i1)
                     *(long)(a_end.i2-a_begin.i2) ,
                      (long)(a_end.i0-a_begin.i0)
                     *(long)(a_end.i1-a_begin.i1)
                     *(long)(a_end.i2-a_begin.i2)
                     *(long)(a_end.i3-a_begin.i3) ,
                      (long)(a_end.i0-a_begin.i0)
                     *(long)(a_end.i1-a_begin.i1)
                     *(long)(a_end.i2-a_begin.i2)
                     *(long)(a_end.i3-a_begin.i3)
                     *(long)(a_end.i4-a_begin.i4) },
#endif
              begin(a_begin),
              end(a_end),
              ncomp(a_ncomp)
            {}

        template <class U,
                  class = typename std::enable_if
                  <std::is_same<typename std::remove_const<T>::type,
                                typename std::remove_const<U>::type>::value>::type >
        AMREX_GPU_HOST_DEVICE
        constexpr ArrayN (ArrayN<U> const& rhs, int start_comp) noexcept
            : p((T*)(rhs.p+start_comp*rhs.stride.i0)),
              stride(rhs.stride),
              begin(rhs.begin),
              end(rhs.end),
              ncomp(rhs.ncomp-start_comp)
        {}

        template <class U=T, class = typename std::enable_if<!std::is_void<U>::value>::type >
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        U& operator() (IntVect const& iv) const noexcept {
#if defined(AMREX_DEBUG) || defined(AMREX_BOUND_CHECK)
            index_assert(iv,0);
#endif
#if (AMREX_SPACEDIM == 1)
            return p[iv[0]-begin.i0];
#elif (AMREX_SPACEDIM == 2)
            return p[(iv[0]-begin.i0)+(iv[1]-begin.i1)*stride.i1];
#elif (AMREX_SPACEDIM == 3)
            return p[(iv[0]-begin.i0)+(iv[1]-begin.i1)*stride.i1+(iv[2]-begin.i2)*stride.i2];
#elif (AMREX_SPACEDIM == 4)
            return p[(iv[0]-begin.i0)+(iv[1]-begin.i1)*stride.i1+(iv[2]-begin.i2)*stride.i2
                     +(iv[3]-begin.i3)*stride.i3];
#elif (AMREX_SPACEDIM == 5)
            return p[(iv[0]-begin.i0)+(iv[1]-begin.i1)*stride.i1+(iv[2]-begin.i2)*stride.i2
                     +(iv[3]-begin.i3)*stride.i3+(iv[4]-begin.i4)*stride.i4];
#elif (AMREX_SPACEDIM == 6)
            return p[(iv[0]-begin.i0)+(iv[1]-begin.i1)*stride.i1+(iv[2]-begin.i2)*stride.i2
                     +(iv[3]-begin.i3)*stride.i3+(iv[4]-begin.i4)*stride.i4+(iv[5]-begin.i5)*stride.i5];
#endif
        }

        template <class U=T, class = typename std::enable_if<!std::is_void<U>::value>::type >
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        U& operator() (IntVect const& iv, int n) const noexcept {
#if defined(AMREX_DEBUG) || defined(AMREX_BOUND_CHECK)
            index_assert(iv,n);
#endif
#if (AMREX_SPACEDIM == 1)
            return p[iv[0]-begin.i0 + n*stride.i0];
#elif (AMREX_SPACEDIM == 2)
            return p[(iv[0]-begin.i0)+(iv[1]-begin.i1)*stride.i1 + n*stride.i0];
#elif (AMREX_SPACEDIM == 3)
            return p[(iv[0]-begin.i0)+(iv[1]-begin.i1)*stride.i1+(iv[2]-begin.i2)*stride.i2
                     + n*stride.i0];
#elif (AMREX_SPACEDIM == 4)
            return p[(iv[0]-begin.i0)+(iv[1]-begin.i1)*stride.i1+(iv[2]-begin.i2)*stride.i2
                     +(iv[3]-begin.i3)*stride.i3 + n*stride.i0];
#elif (AMREX_SPACEDIM == 5)
            return p[(iv[0]-begin.i0)+(iv[1]-begin.i1)*stride.i1+(iv[2]-begin.i2)*stride.i2
                     +(iv[3]-begin.i3)*stride.i3+(iv[4]-begin.i4)*stride.i4 + n*stride.i0];
#elif (AMREX_SPACEDIM == 6)
            return p[(iv[0]-begin.i0)+(iv[1]-begin.i1)*stride.i1+(iv[2]-begin.i2)*stride.i2
                     +(iv[3]-begin.i3)*stride.i3+(iv[4]-begin.i4)*stride.i4+(iv[5]-begin.i5)*stride.i5
                     + n*stride.i0];
#endif
        }

        template <class U=T, class = typename std::enable_if<!std::is_void<U>::value>::type >
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        T* ptr (IntVect const& iv) const noexcept {
#if defined(AMREX_DEBUG) || defined(AMREX_BOUND_CHECK)
            index_assert(iv,0);
#endif
#if (AMREX_SPACEDIM == 1)
            return p + (iv[0]-begin.i0);
#elif (AMREX_SPACEDIM == 2)
            return p + ((iv[0]-begin.i0)+(iv[1]-begin.i1)*stride.i1);
#elif (AMREX_SPACEDIM == 3)
            return p + ((iv[0]-begin.i0)+(iv[1]-begin.i1)*stride.i1+(iv[2]-begin.i2)*stride.i2);
#elif (AMREX_SPACEDIM == 4)
            return p + ((iv[0]-begin.i0)+(iv[1]-begin.i1)*stride.i1+(iv[2]-begin.i2)*stride.i2
                        +(iv[3]-begin.i3)*stride.i3);
#elif (AMREX_SPACEDIM == 5)
            return p + ((iv[0]-begin.i0)+(iv[1]-begin.i1)*stride.i1+(iv[2]-begin.i2)*stride.i2
                     +(iv[3]-begin.i3)*stride.i3+(iv[4]-begin.i4)*stride.i4);
#elif (AMREX_SPACEDIM == 6)
            return p + ((iv[0]-begin.i0)+(iv[1]-begin.i1)*stride.i1+(iv[2]-begin.i2)*stride.i2
                     +(iv[3]-begin.i3)*stride.i3+(iv[4]-begin.i4)*stride.i4+(iv[5]-begin.i5)*stride.i5);
#endif
        }

        template <class U=T, class = typename std::enable_if<!std::is_void<U>::value>::type >
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        T* ptr (IntVect const& iv, int n) const noexcept {
#if defined(AMREX_DEBUG) || defined(AMREX_BOUND_CHECK)
            index_assert(iv,n);
#endif
#if (AMREX_SPACEDIM == 1)
            return p + (iv[0]-begin.i0 + n*stride.i0);
#elif (AMREX_SPACEDIM == 2)
            return p + ((iv[0]-begin.i0)+(iv[1]-begin.i1)*stride.i1 + n*stride.i0);
#elif (AMREX_SPACEDIM == 3)
            return p + ((iv[0]-begin.i0)+(iv[1]-begin.i1)*stride.i1+(iv[2]-begin.i2)*stride.i2
                     + n*stride.i0);
#elif (AMREX_SPACEDIM == 4)
            return p + ((iv[0]-begin.i0)+(iv[1]-begin.i1)*stride.i1+(iv[2]-begin.i2)*stride.i2
                     +(iv[3]-begin.i3)*stride.i3 + n*stride.i0);
#elif (AMREX_SPACEDIM == 5)
            return p + ((iv[0]-begin.i0)+(iv[1]-begin.i1)*stride.i1+(iv[2]-begin.i2)*stride.i2
                     +(iv[3]-begin.i3)*stride.i3+(iv[4]-begin.i4)*stride.i4 + n*stride.i0);
#elif (AMREX_SPACEDIM == 6)
            return p + ((iv[0]-begin.i0)+(iv[1]-begin.i1)*stride.i1+(iv[2]-begin.i2)*stride.i2
                     +(iv[3]-begin.i3)*stride.i3+(iv[4]-begin.i4)*stride.i4+(iv[5]-begin.i5)*stride.i5
                     + n*stride.i0);
#endif
        }

        template <class U=T, class = typename std::enable_if<!std::is_void<U>::value>::type >
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        U& operator() (DimN const& iv) const noexcept {
#if defined(AMREX_DEBUG) || defined(AMREX_BOUND_CHECK)
            index_assert(iv,0);
#endif
#if (AMREX_SPACEDIM == 1)
            return p[iv.i0-begin.i0];
#elif (AMREX_SPACEDIM == 2)
            return p[(iv.i0-begin.i0)+(iv.i1-begin.i1)*stride.i1];
#elif (AMREX_SPACEDIM == 3)
            return p[(iv.i0-begin.i0)+(iv.i1-begin.i1)*stride.i1+(iv.i2-begin.i2)*stride.i2];
#elif (AMREX_SPACEDIM == 4)
            return p[(iv.i0-begin.i0)+(iv.i1-begin.i1)*stride.i1+(iv.i2-begin.i2)*stride.i2
                     +(iv.i3-begin.i3)*stride.i3];
#elif (AMREX_SPACEDIM == 5)
            return p[(iv.i0-begin.i0)+(iv.i1-begin.i1)*stride.i1+(iv.i2-begin.i2)*stride.i2
                     +(iv.i3-begin.i3)*stride.i3+(iv.i4-begin.i4)*stride.i4];
#elif (AMREX_SPACEDIM == 6)
            return p[(iv.i0-begin.i0)+(iv.i1-begin.i1)*stride.i1+(iv.i2-begin.i2)*stride.i2
                     +(iv.i3-begin.i3)*stride.i3+(iv.i4-begin.i4)*stride.i4+(iv.i5-begin.i5)*stride.i5];
#endif
        }

        template <class U=T, class = typename std::enable_if<!std::is_void<U>::value>::type >
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        U& operator() (DimN const& iv, int n) const noexcept {
#if defined(AMREX_DEBUG) || defined(AMREX_BOUND_CHECK)
            index_assert(iv,n);
#endif
#if (AMREX_SPACEDIM == 1)
            return p[iv.i0-begin.i0 + n*stride.i0];
#elif (AMREX_SPACEDIM == 2)
            return p[(iv.i0-begin.i0)+(iv.i1-begin.i1)*stride.i1 + n*stride.i0];
#elif (AMREX_SPACEDIM == 3)
            return p[(iv.i0-begin.i0)+(iv.i1-begin.i1)*stride.i1+(iv.i2-begin.i2)*stride.i2
                     + n*stride.i0];
#elif (AMREX_SPACEDIM == 4)
            return p[(iv.i0-begin.i0)+(iv.i1-begin.i1)*stride.i1+(iv.i2-begin.i2)*stride.i2
                     +(iv.i3-begin.i3)*stride.i3 + n*stride.i0];
#elif (AMREX_SPACEDIM == 5)
            return p[(iv.i0-begin.i0)+(iv.i1-begin.i1)*stride.i1+(iv.i2-begin.i2)*stride.i2
                     +(iv.i3-begin.i3)*stride.i3+(iv.i4-begin.i4)*stride.i4 + n*stride.i0];
#elif (AMREX_SPACEDIM == 6)
            return p[(iv.i0-begin.i0)+(iv.i1-begin.i1)*stride.i1+(iv.i2-begin.i2)*stride.i2
                     +(iv.i3-begin.i3)*stride.i3+(iv.i4-begin.i4)*stride.i4+(iv.i5-begin.i5)*stride.i5
                     + n*stride.i0];
#endif
        }

        template <class U=T, class = typename std::enable_if<!std::is_void<U>::value>::type >
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        T* ptr (DimN const& iv) const noexcept {
#if defined(AMREX_DEBUG) || defined(AMREX_BOUND_CHECK)
            index_assert(iv,0);
#endif
#if (AMREX_SPACEDIM == 1)
            return p + (iv.i0-begin.i0);
#elif (AMREX_SPACEDIM == 2)
            return p + ((iv.i0-begin.i0)+(iv.i1-begin.i1)*stride.i1);
#elif (AMREX_SPACEDIM == 3)
            return p + ((iv.i0-begin.i0)+(iv.i1-begin.i1)*stride.i1+(iv.i2-begin.i2)*stride.i2);
#elif (AMREX_SPACEDIM == 4)
            return p + ((iv.i0-begin.i0)+(iv.i1-begin.i1)*stride.i1+(iv.i2-begin.i2)*stride.i2
                        +(iv.i3-begin.i3)*stride.i3);
#elif (AMREX_SPACEDIM == 5)
            return p + ((iv.i0-begin.i0)+(iv.i1-begin.i1)*stride.i1+(iv.i2-begin.i2)*stride.i2
                     +(iv.i3-begin.i3)*stride.i3+(iv.i4-begin.i4)*stride.i4);
#elif (AMREX_SPACEDIM == 6)
            return p + ((iv.i0-begin.i0)+(iv.i1-begin.i1)*stride.i1+(iv.i2-begin.i2)*stride.i2
                     +(iv.i3-begin.i3)*stride.i3+(iv.i4-begin.i4)*stride.i4+(iv.i5-begin.i5)*stride.i5);
#endif
        }

        template <class U=T, class = typename std::enable_if<!std::is_void<U>::value>::type >
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        T* ptr (DimN const& iv, int n) const noexcept {
#if defined(AMREX_DEBUG) || defined(AMREX_BOUND_CHECK)
            index_assert(iv,n);
#endif
#if (AMREX_SPACEDIM == 1)
            return p + (iv.i0-begin.i0 + n*stride.i0);
#elif (AMREX_SPACEDIM == 2)
            return p + ((iv.i0-begin.i0)+(iv.i1-begin.i1)*stride.i1 + n*stride.i0);
#elif (AMREX_SPACEDIM == 3)
            return p + ((iv.i0-begin.i0)+(iv.i1-begin.i1)*stride.i1+(iv.i2-begin.i2)*stride.i2
                     + n*stride.i0);
#elif (AMREX_SPACEDIM == 4)
            return p + ((iv.i0-begin.i0)+(iv.i1-begin.i1)*stride.i1+(iv.i2-begin.i2)*stride.i2
                     +(iv.i3-begin.i3)*stride.i3 + n*stride.i0);
#elif (AMREX_SPACEDIM == 5)
            return p + ((iv.i0-begin.i0)+(iv.i1-begin.i1)*stride.i1+(iv.i2-begin.i2)*stride.i2
                     +(iv.i3-begin.i3)*stride.i3+(iv.i4-begin.i4)*stride.i4 + n*stride.i0);
#elif (AMREX_SPACEDIM == 6)
            return p + ((iv.i0-begin.i0)+(iv.i1-begin.i1)*stride.i1+(iv.i2-begin.i2)*stride.i2
                     +(iv.i3-begin.i3)*stride.i3+(iv.i4-begin.i4)*stride.i4+(iv.i5-begin.i5)*stride.i5
                     + n*stride.i0);
#endif
        }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        T* dataPtr () const noexcept {
            return this->p;
        }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        std::size_t size () const noexcept {
            return this->stride.i0 * this->ncomp;
        }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        int nComp () const noexcept { return ncomp; }

#if defined(AMREX_DEBUG) || defined(AMREX_BOUND_CHECK)
        AMREX_GPU_HOST_DEVICE inline
        void index_assert (DimN const& iv, int n) const
        {
#if (AMREX_SPACEDIM == 1)
            if (iv.i0 < begin.i0 || iv.i0 >= end.i0 ||
                n < 0 || n >= ncomp) {
#ifdef AMREX_DEVICE_COMPILE
                std::printf(" (%d,%d) is out of bound (%d:%d,0:%d)\n",
                            iv.i0, n,
                            begin.i0, end.i0-1,
                            ncomp-1);
                amrex::Abort();
#else
                std::stringstream ss;
                ss << " (" << iv.i0 << "," <<  n
                   << ") is out of bound ("
                   << begin.i0 << ":" << end.i0-1 << ","
                   << "0:" << ncomp-1 << ")";
                amrex::Abort(ss.str());
#endif
            }
#elif (AMREX_SPACEDIM == 2)
            if (iv.i0 < begin.i0 || iv.i0 >= end.i0 ||
                iv.i1 < begin.i1 || iv.i1 >= end.i1 ||
                n < 0 || n >= ncomp) {
#ifdef AMREX_DEVICE_COMPILE
                std::printf(" (%d,%d,%d) is out of bound (%d:%d,%d:%d,0:%d)\n",
                            iv.i0, iv.i1, n,
                            begin.i0, end.i0-1,
                            begin.i1, end.i1-1,
                            ncomp-1);
                amrex::Abort();
#else
                std::stringstream ss;
                ss << " (" << iv.i0 << "," << iv.i1 << "," <<  n
                   << ") is out of bound ("
                   << begin.i0 << ":" << end.i0-1 << ","
                   << begin.i1 << ":" << end.i1-1 << ","
                   << "0:" << ncomp-1 << ")";
                amrex::Abort(ss.str());
#endif
            }
#elif (AMREX_SPACEDIM == 3)
            if (iv.i0 < begin.i0 || iv.i0 >= end.i0 ||
                iv.i1 < begin.i1 || iv.i1 >= end.i1 ||
                iv.i2 < begin.i2 || iv.i2 >= end.i2 ||
                n < 0 || n >= ncomp) {
#ifdef AMREX_DEVICE_COMPILE
                std::printf(" (%d,%d,%d,%d) is out of bound (%d:%d,%d:%d,%d:%d,0:%d)\n",
                            iv.i0, iv.i1, iv.i2, n,
                            begin.i0, end.i0-1,
                            begin.i1, end.i1-1,
                            begin.i2, end.i2-1,
                            ncomp-1);
                amrex::Abort();
#else
                std::stringstream ss;
                ss << " (" << iv.i0 << "," << iv.i1 << "," << iv.i2 << "," <<  n
                   << ") is out of bound ("
                   << begin.i0 << ":" << end.i0-1 << ","
                   << begin.i1 << ":" << end.i1-1 << ","
                   << begin.i2 << ":" << end.i2-1 << ","
                   << "0:" << ncomp-1 << ")";
                amrex::Abort(ss.str());
#endif
            }
#elif (AMREX_SPACEDIM == 4)
            if (iv.i0 < begin.i0 || iv.i0 >= end.i0 ||
                iv.i1 < begin.i1 || iv.i1 >= end.i1 ||
                iv.i2 < begin.i2 || iv.i2 >= end.i2 ||
                iv.i3 < begin.i3 || iv.i3 >= end.i3 ||
                n < 0 || n >= ncomp) {
#ifdef AMREX_DEVICE_COMPILE
                std::printf(" (%d,%d,%d,%d,%d) is out of bound (%d:%d,%d:%d,%d:%d,%d:%d,0:%d)\n",
                            iv.i0, iv.i1, iv.i2, iv.i3, n,
                            begin.i0, end.i0-1,
                            begin.i1, end.i1-1,
                            begin.i2, end.i2-1,
                            begin.i3, end.i3-1,
                            ncomp-1);
                amrex::Abort();
#else
                std::stringstream ss;
                ss << " (" << iv.i0 << "," << iv.i1 << "," << iv.i2 << "," << iv.i3 << "," <<  n
                   << ") is out of bound ("
                   << begin.i0 << ":" << end.i0-1 << ","
                   << begin.i1 << ":" << end.i1-1 << ","
                   << begin.i2 << ":" << end.i2-1 << ","
                   << begin.i3 << ":" << end.i3-1 << ","
                   << "0:" << ncomp-1 << ")";
                amrex::Abort(ss.str());
#endif
            }
#elif (AMREX_SPACEDIM == 5)
            if (iv.i0 < begin.i0 || iv.i0 >= end.i0 ||
                iv.i1 < begin.i1 || iv.i1 >= end.i1 ||
                iv.i2 < begin.i2 || iv.i2 >= end.i2 ||
                iv.i3 < begin.i3 || iv.i3 >= end.i3 ||
                iv.i4 < begin.i4 || iv.i4 >= end.i4 ||
                n < 0 || n >= ncomp) {
#ifdef AMREX_DEVICE_COMPILE
                std::printf(" (%d,%d,%d,%d,%d,%d) is out of bound (%d:%d,%d:%d,%d:%d,%d:%d,%d:%d,0:%d)\n",
                            iv.i0, iv.i1, iv.i2, iv.i3, iv.i4, n,
                            begin.i0, end.i0-1,
                            begin.i1, end.i1-1,
                            begin.i2, end.i2-1,
                            begin.i3, end.i3-1,
                            begin.i4, end.i4-1,
                            ncomp-1);
                amrex::Abort();
#else
                std::stringstream ss;
                ss << " (" << iv.i0 << "," << iv.i1 << "," << iv.i2 << "," << iv.i3 << ","
                   << iv.i4 << ", " << n << ") is out of bound ("
                   << begin.i0 << ":" << end.i0-1 << ","
                   << begin.i1 << ":" << end.i1-1 << ","
                   << begin.i2 << ":" << end.i2-1 << ","
                   << begin.i3 << ":" << end.i3-1 << ","
                   << begin.i4 << ":" << end.i4-1 << ","
                   << "0:" << ncomp-1 << ")";
                amrex::Abort(ss.str());
#endif
            }
#elif (AMREX_SPACEDIM == 6)
            if (iv.i0 < begin.i0 || iv.i0 >= end.i0 ||
                iv.i1 < begin.i1 || iv.i1 >= end.i1 ||
                iv.i2 < begin.i2 || iv.i2 >= end.i2 ||
                iv.i3 < begin.i3 || iv.i3 >= end.i3 ||
                iv.i4 < begin.i4 || iv.i4 >= end.i4 ||
                iv.i5 < begin.i5 || iv.i5 >= end.i5 ||
                n < 0 || n >= ncomp) {
#ifdef AMREX_DEVICE_COMPILE
                std::printf(" (%d,%d,%d,%d,%d,%d,%d) is out of bound (%d:%d,%d:%d,%d:%d,%d:%d,%d:%d,%d:%d,0:%d)\n",
                            iv.i0, iv.i1, iv.i2, iv.i3, iv.i4, iv.i5, n,
                            begin.i0, end.i0-1,
                            begin.i1, end.i1-1,
                            begin.i2, end.i2-1,
                            begin.i3, end.i3-1,
                            begin.i4, end.i4-1,
                            begin.i5, end.i5-1,
                            ncomp-1);
                amrex::Abort();
#else
                std::stringstream ss;
                ss << " (" << iv.i0 << "," << iv.i1 << "," << iv.i2 << "," << iv.i3 << ","
                   << iv.i4 << ", " << iv.i5 << ", " << n << ") is out of bound ("
                   << begin.i0 << ":" << end.i0-1 << ","
                   << begin.i1 << ":" << end.i1-1 << ","
                   << begin.i2 << ":" << end.i2-1 << ","
                   << begin.i3 << ":" << end.i3-1 << ","
                   << begin.i4 << ":" << end.i4-1 << ","
                   << begin.i5 << ":" << end.i5-1 << ","
                   << "0:" << ncomp-1 << ")";
                amrex::Abort(ss.str());
#endif
            }
#endif
        }

        AMREX_GPU_HOST_DEVICE inline
        void index_assert (IntVect const& iv, int n) const
        {
            index_assert(DimN{AMREX_D6_DECL(iv[0],iv[1],iv[2],iv[3],iv[4],iv[5])},n);
        }
#endif
    };

    template <class T>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    DimN lboundn (ArrayN<T> const& a) noexcept
    {
        return a.begin;
    }

    template <class T>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    DimN uboundn (ArrayN<T> const& a) noexcept
    {
        return DimN{AMREX_D6_DECL(a.end.i0-1,a.end.i1-1,a.end.i2-1,a.end.i3-1,a.end.i4-1,a.end.i5-1)};
    }

    template <class T>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    DimN lengthn (ArrayN<T> const& a) noexcept
    {
        return DimN{AMREX_D6_DECL(a.end.i0-a.begin.i0,a.end.i1-a.begin.i1,a.end.i2-a.begin.i2,
                                  a.end.i3-a.begin.i3,a.end.i4-a.begin.i4,a.end.i5-a.begin.i5)};
    }

    template <typename T>
    AMREX_GPU_HOST_DEVICE
    ArrayN<T> shift (ArrayN<T> const& src, int comp) noexcept
    {
        ArrayN<T> r = src;
        r.p += comp * src.stride.i0;
        r.ncomp -= comp;
        return r;
    }

    template <typename T>
    AMREX_GPU_HOST_DEVICE
    ArrayN<T> shift (ArrayN<T> const& src, IntVect const& off, int comp) noexcept
    {
        ArrayN<T> r = src;
        r.p += comp * src.stride.i0;
        r.ncomp -= comp;
        AMREX_D6_TERM(r.begin.i0 += off[0];, r.begin.i1 += off[1];, r.begin.i2 += off[2];,
                      r.begin.i3 += off[3];, r.begin.i4 += off[4];, r.begin.i5 += off[5];);
        AMREX_D6_TERM(  r.end.i0 += off[0];,   r.end.i1 += off[1];,   r.end.i2 += off[2];,
                        r.end.i3 += off[3];,   r.end.i4 += off[4];,   r.end.i5 += off[5];);
        return r;
    }

    template <typename T>
    std::ostream& operator<< (std::ostream& os, const ArrayN<T>& a) {
        os << "((" << lboundn(a) << ',' << uboundn(a) << ")," << a.ncomp << ')';
        return os;
    }
}

#endif
