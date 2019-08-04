#ifndef AMREX_SCAN_H_
#define AMREX_SCAN_H_

#include <AMReX_Gpu.H>
#include <AMReX_Arena.H>
#include <cstdint>
#include <type_traits>

namespace amrex {
namespace Scan {

#ifdef AMREX_USE_GPU

namespace detail {

template <typename T>
struct STVA
{
    char status;
    T value;
};

template <typename T, bool SINGLE_WORD> struct BlockStatus {};

template <typename T>
struct BlockStatus<T, true>
{
    template<typename U>
    union Data {
        STVA<U> s;
        uint64_t i;
        void operator=(Data<U> const&) = delete;
        void operator=(Data<U> &&) = delete;
    };
    Data<T> d;

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    void write (char a_status, T a_value) {
        volatile uint64_t tmp;
        reinterpret_cast<STVA<T> volatile&>(tmp).status = a_status;
        reinterpret_cast<STVA<T> volatile&>(tmp).value  = a_value;
        d.i = tmp;
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    T get_aggregate() const { return d.s.value; }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    STVA<T> read () volatile {
        volatile uint64_t tmp = d.i;
        return {reinterpret_cast<STVA<T> volatile&>(tmp).status,
                reinterpret_cast<STVA<T> volatile&>(tmp).value };
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    void set_status (char a_status) { d.s.status = a_status; }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    STVA<T> wait () volatile {
        STVA<T> r;
        do {
            __threadfence_block();
            r = read();
        } while (r.status == 'x');
        return r;
    }
};

template <typename T>
struct BlockStatus<T, false>
{
    T aggregate;
    T inclusive;
    char status;

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    void write (char a_status, T a_value) {
        if (a_status == 'a') {
            aggregate = a_value;
        } else {
            inclusive = a_value;
        }
        __threadfence();
        status = a_status;
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    T get_aggregate() const { return aggregate; }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    STVA<T> read () volatile {
        if (status == 'x') {
            return {'x', 0};
        } else if (status == 'a') {
            return {'a', aggregate};
        } else {
            return {'p', inclusive};
        }
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    void set_status (char a_status) { status = a_status; }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    STVA<T> wait () volatile {
        STVA<T> r;
        do {
            r = read();
            __threadfence();
        } while (r.status == 'x');
        return r;
    }
};

template <int N, typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void StripedToBlocked (T* items, T* temp)
{
    #pragma unroll
    for (int i = 0; i < N; ++i) {
        temp[i*blockDim.x + threadIdx.x] = items[i];
    }
    __syncthreads();
    #pragma unroll
    for (int i = 0; i < N; ++i) {
        items[i] = temp[threadIdx.x*N + i];
    }
}

template <int N, typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void BlockedToStriped (T* items, T* temp)
{
    #pragma unroll
    for (int i = 0; i < N; ++i) {
        temp[threadIdx.x*N + i] = items[i];
    }
    __syncthreads();
    #pragma unroll
    for (int i = 0; i < N; ++i) {
        items[i] = temp[i*blockDim.x+threadIdx.x];
    }
}

}

enum class Type { inclusive, exclusive };

template <typename T, typename FIN, typename FOUT>
T PrefixSum (int n, FIN && fin, FOUT && fout, Type type)
{
    if (n <= 0) return 0;
    constexpr int nwarps = 4;
    constexpr int nthreads = nwarps*Gpu::Device::warp_size;
    constexpr int nitems = 12;
    constexpr int nelms_per_block = nthreads * nitems;
    int nblocks = (n + nelms_per_block - 1) / nelms_per_block;
    std::size_t sm = sizeof(T) * (nthreads*nitems + nwarps+1) + sizeof(int);
    auto stream = Gpu::gpuStream();

    typedef typename std::conditional<sizeof(detail::STVA<T>) <= 8,
        detail::BlockStatus<T,true>, detail::BlockStatus<T,false> >::type BlockStatusT;

    // temporary memory
    std::size_t nbytes_blockstatus = Arena::align(sizeof(BlockStatusT)*nblocks);
    std::size_t nbytes_blockid = Arena::align(sizeof(unsigned int));
    std::size_t nbytes_totalsum = Arena::align(sizeof(T));
    auto dp = (char*)(The_Device_Arena()->alloc(  nbytes_blockstatus
                                                + nbytes_blockid
                                                + nbytes_totalsum));
    BlockStatusT* AMREX_RESTRICT block_status_p = (BlockStatusT*)dp;
    unsigned int* AMREX_RESTRICT virtual_block_id_p = (unsigned int*)(dp + nbytes_blockstatus);
    T* AMREX_RESTRICT totalsum_p = (T*)(dp + nbytes_blockstatus + nbytes_blockid);

    amrex::ParallelFor(nblocks, [=] AMREX_GPU_DEVICE (int i) noexcept {
        BlockStatusT& block_status = block_status_p[i];
        block_status.set_status('x');
        if (i == 0) {
            *virtual_block_id_p = 0;
            *totalsum_p = 0;
        }
    });

    amrex::launch_global<<<nblocks, nthreads, sm, stream>>>(
    [=] AMREX_GPU_DEVICE () noexcept
    {
        int lane = threadIdx.x % Gpu::Device::warp_size;
        int warp = threadIdx.x / Gpu::Device::warp_size;

        amrex::Gpu::SharedMemory<T> gsm;
        T* shared = gsm.dataPtr();
        T* shared2 = shared + blockDim.x*nitems;

        // First of all, get block virtual id.  We must do this to
        // avoid deadlock because CUDA may launch blocks in any order.
        // Anywhere in this function, we should not use blockIdx.
        int virtual_block_id = 0;
        if (gridDim.x > 1) {
            int& virtual_block_id_shared = *((int*)(shared+blockDim.x*nitems+nwarps+1));
            if (threadIdx.x == 0) {
                unsigned int bid = Gpu::Atomic::Inc(virtual_block_id_p, gridDim.x);
                virtual_block_id_shared = bid;
            }
            __syncthreads();
            virtual_block_id = virtual_block_id_shared;
        }

        // Each block processes [ibegin,iend).
        int ibegin = nelms_per_block * virtual_block_id;
        int iend = amrex::min(ibegin+nelms_per_block, n);
        BlockStatusT& block_status = block_status_p[virtual_block_id];

        //
        // This is based on "Single-pass Parallel Prefix Scan with
        // Decoupled Look-back" by D. Merrill & M. Garland.
        //

        // Each block is responsible for threadDim.x*nitems elements.
        // Each threads is responsible for nitems elements.
        T items[nitems];
        for (int i = 0; i < nitems; ++i) {
            int offset = ibegin + i*blockDim.x + threadIdx.x;
            if (offset < iend) {
                items[i] = fin(offset);
            } else {
                items[i] = 0;
            }
        }
        detail::StripedToBlocked<nitems>(items, shared);

        // Inclusive sum within a thread
        #pragma unroll
        for (int i = 1; i < nitems; ++i) {
            items[i] += items[i-1];
        }

        // Inclusive sum within a warp
        T x = items[nitems-1];
        for (int i = 1; i <= Gpu::Device::warp_size; i *= 2) {
            T s = __shfl_up_sync(0xffffffff, x, i);
            if (lane >= i) x += s;
        }

        T sum_prev_threads = x - items[nitems-1];

        // Store the sum of warp in shared memory.
        if (lane == Gpu::Device::warp_size-1) {
            shared2[warp] = x;
        }
        __syncthreads();

        T sum_prev_warps = 0;
        if (warp > 0) {
            sum_prev_warps = shared2[0];
            for (int i = 1; i < warp; ++i) {
                sum_prev_warps += shared2[i];
            }
        }
        bool last_thread_in_block = threadIdx.x == blockDim.x-1;
        // inclusive_sum_block is only valid for the last thread in the block
        T inclusive_sum_block;
        if (last_thread_in_block) {
            inclusive_sum_block = items[nitems-1] + sum_prev_warps + sum_prev_threads;
            block_status.write((virtual_block_id == 0) ? 'p' : 'a', inclusive_sum_block);
        }

        T sum_prev_blocks = 0;

        if (virtual_block_id > 0) {
            if (warp == nwarps-1) {
                T exclusive_prefix = 0;
                BlockStatusT volatile* pbs = block_status_p;
                for (int iblock0 = virtual_block_id-1; iblock0 >= 0; iblock0 -= Gpu::Device::warp_size)
                {
                    int iblock = iblock0-lane;
                    detail::STVA<T> stva{'p', 0};
                    if (iblock >= 0) {
                        stva = pbs[iblock].wait();
                    }

                    T x = stva.value;

                    unsigned const status_bf = __ballot_sync(0xffffffff, stva.status == 'p');
                    bool stop_lookback = status_bf & 0x1u;
                    if (stop_lookback == false) {
                        if (status_bf != 0) {
                            T y = x;
                            if (lane > 0) x = 0;
                            unsigned int bit_mask = 0x1u;
                            for (int i = 1; i < Gpu::Device::warp_size; ++i) {
                                bit_mask <<= 1;
                                if (i == lane) x = y;
                                if (status_bf & bit_mask) {
                                    stop_lookback = true;
                                    break;
                                }
                            }
                        }

                        for (int i = Gpu::Device::warp_size/2; i > 0; i /= 2) {
                            x += __shfl_down_sync(0xffffffff, x, i);
                        }
                    }

                    if (lane == 0) { exclusive_prefix += x; }
                    if (stop_lookback) break;
                }

                if (lane == 0) {
                    block_status.write('p', block_status.get_aggregate() + exclusive_prefix);
                    shared2[nwarps] = exclusive_prefix;
                }
            }

            __syncthreads();

            sum_prev_blocks = shared2[nwarps];
        }

        T s = sum_prev_warps + sum_prev_threads;
        if (type == Type::inclusive) {
            #pragma unroll
            for (int i = 0; i < nitems; ++i) {
                items[i] += s;
            }
        } else {
            #pragma unroll
            for (int i = nitems-1; i > 0; --i) {
                items[i] = items[i-1] + s;
            }
            items[0] = s;
        }

        detail::BlockedToStriped<nitems>(items, shared);
        for (int i = 0; i < nitems; ++i) {
            int offset = ibegin + i*blockDim.x + threadIdx.x;
            if (offset >= iend) break;
            fout(offset, items[i]+sum_prev_blocks);
        }

        if (virtual_block_id == gridDim.x-1 && last_thread_in_block) {
            *totalsum_p = sum_prev_blocks + inclusive_sum_block;
        }
    });

    T totalsum;
    Gpu::dtoh_memcpy(&totalsum, totalsum_p, sizeof(T));

    The_Device_Arena()->free(dp);

    AMREX_GPU_ERROR_CHECK();

    return totalsum;
}

// The return value is the total sum.
template <typename N, typename T, typename M=amrex::EnableIf_t<std::is_integral<N>::value> >
T InclusiveSum (N n, T const* in, T * out)
{
    AMREX_ALWAYS_ASSERT(static_cast<long>(n) < static_cast<long>(std::numeric_limits<int>::max()));
    return PrefixSum<T>(n,
                 [=] AMREX_GPU_DEVICE (int i) -> T { return in[i]; },
                 [=] AMREX_GPU_DEVICE (int i, T const& x) { out[i] = x; },
                 Type::inclusive);
}

// The return value is the total sum.
template <typename N, typename T, typename M=amrex::EnableIf_t<std::is_integral<N>::value> >
T ExclusiveSum (N n, T const* in, T * out)
{
    AMREX_ALWAYS_ASSERT(static_cast<long>(n) < static_cast<long>(std::numeric_limits<int>::max()));
    return PrefixSum<T>(n,
                 [=] AMREX_GPU_DEVICE (int i) -> T { return in[i]; },
                 [=] AMREX_GPU_DEVICE (int i, T const& x) { out[i] = x; },
                 Type::exclusive);
}

#endif

}}

#endif
