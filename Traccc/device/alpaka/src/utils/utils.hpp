/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include <alpaka/alpaka.hpp>

namespace traccc::alpaka {

using Dim = ::alpaka::DimInt<1>;
using Idx = uint32_t;
using WorkDiv = ::alpaka::WorkDivMembers<Dim, Idx>;

// Get alpaka accelerator - based on alpaka/examples/ExampleDefaultAcc.hpp
#if defined(ALPAKA_ACC_GPU_CUDA_ENABLED)
using Acc = ::alpaka::AccGpuCudaRt<Dim, Idx>;
#elif defined(ALPAKA_ACC_GPU_HIP_ENABLED)
using Acc = ::alpaka::AccGpuHipRt<Dim, Idx>;
#elif defined(ALPAKA_ACC_SYCL_ENABLED)
#if defined(ALPAKA_SYCL_ONEAPI_CPU)
using Acc = ::alpaka::AccCpuSycl<Dim, Idx>;
#elif defined(ALPAKA_SYCL_ONEAPI_FPGA)
using Acc = ::alpaka::AccFpgaSyclIntel<Dim, Idx>;
#elif defined(ALPAKA_SYCL_ONEAPI_GPU)
using Acc = ::alpaka::AccGpuSyclIntel<Dim, Idx>;
#endif
#elif defined(ALPAKA_ACC_CPU_B_SEQ_T_THREADS_ENABLED)
using Acc = ::alpaka::AccCpuThreads<Dim, Idx>;
#else
#error "No supported backend selected." //we definitely want to fail the build if no matching accelerator is found
#endif

// Set the queue type to blocking for CPU based accelerators, due to the
// synchronous nature of vecmem there.
#if defined(ALPAKA_ACC_CPU_B_SEQ_T_THREADS_ENABLED)
using QueueType = ::alpaka::Blocking;
#else
using QueueType = ::alpaka::NonBlocking;
#endif

using Host = ::alpaka::DevCpu;
using Queue = ::alpaka::Queue<Acc, QueueType>;

template <typename TAcc>
consteval Idx getWarpSize() {
    if constexpr (::alpaka::accMatchesTags<TAcc, ::alpaka::TagGpuCudaRt,
                                           ::alpaka::TagGpuSyclIntel>) {
        return 32;
    }
    if constexpr (::alpaka::accMatchesTags<TAcc, ::alpaka::TagGpuHipRt>) {
        return 64;
    } else {
        return 4;
    }
}

template <typename TAcc>
inline WorkDiv makeWorkDiv(Idx blocks, Idx threadsOrElements) {
    const Idx blocksPerGrid = std::max(Idx{1}, blocks);
    if constexpr (::alpaka::isMultiThreadAcc<TAcc>) {
        const Idx threadsPerBlock(threadsOrElements);
        const Idx elementsPerThread = Idx{1};
        return WorkDiv{blocksPerGrid, threadsPerBlock, elementsPerThread};
    } else {
        const Idx threadsPerBlock = Idx{1};
        const Idx elementsPerThread(threadsOrElements);
        return WorkDiv{blocksPerGrid, threadsPerBlock, elementsPerThread};
    }
}

}  // namespace traccc::alpaka
