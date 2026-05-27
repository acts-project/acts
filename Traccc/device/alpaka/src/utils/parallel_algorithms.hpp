/**
 * traccc library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "utils.hpp"

// Project include(s).
#include "traccc/utils/memory_resource.hpp"

// Thrust include(s).
#if !defined(ALPAKA_ACC_SYCL_ENABLED)
#include <thrust/binary_search.h>
#include <thrust/copy.h>
#include <thrust/execution_policy.h>
#include <thrust/fill.h>
#include <thrust/scan.h>
#include <thrust/sort.h>
#include <thrust/unique.h>
#endif

// OneDPL include.
//
// This is left to a separate file to turn off warnings from oneDPL.
#if defined(ALPAKA_ACC_SYCL_ENABLED)
#include "oneDPL.hpp"
#endif

namespace traccc::alpaka::details {

inline auto getExecutionPolicy([[maybe_unused]] Queue &q,
                               [[maybe_unused]] const memory_resource &mr) {
#if defined(ALPAKA_ACC_GPU_CUDA_ENABLED)
    auto stream = ::alpaka::getNativeHandle(q);
    return thrust::cuda::par_nosync(
               std::pmr::polymorphic_allocator<std::byte>(&(mr.main)))
        .on(stream);
#elif defined(ALPAKA_ACC_GPU_HIP_ENABLED)
    auto stream = ::alpaka::getNativeHandle(q);
    return thrust::hip_rocprim::par_nosync(
               std::pmr::polymorphic_allocator<std::byte>(&(mr.main)))
        .on(stream);
#elif defined(ALPAKA_ACC_SYCL_ENABLED)
    auto queue = ::alpaka::getNativeHandle(q);
    return oneapi::dpl::execution::device_policy{queue};
#else
    return thrust::host;
#endif
}

template <typename RandomAccessIterator, typename Compare>
void sort(Queue &q, const memory_resource mr, RandomAccessIterator first,
          RandomAccessIterator last, Compare comp) {
    auto execPolicy = getExecutionPolicy(q, mr);

#if defined(ALPAKA_ACC_SYCL_ENABLED)
    oneapi::dpl::sort(execPolicy, first, last, comp);
#else
    thrust::sort(execPolicy, first, last, comp);
#endif
}

template <typename RandomAccessIterator1, typename RandomAccessIterator2,
          typename Compare>
void sort_by_key(Queue &q, const memory_resource &mr,
                 RandomAccessIterator1 keys_first,
                 RandomAccessIterator1 keys_last,
                 RandomAccessIterator2 values_first, Compare comp) {
    auto execPolicy = getExecutionPolicy(q, mr);

#if defined(ALPAKA_ACC_SYCL_ENABLED)
    oneapi::dpl::sort_by_key(execPolicy, keys_first, keys_last, values_first,
                             comp);
#else
    thrust::sort_by_key(execPolicy, keys_first, keys_last, values_first, comp);
#endif
}

template <typename RandomAccessIterator1, typename RandomAccessIterator2>
void sort_by_key(Queue &q, const memory_resource &mr,
                 RandomAccessIterator1 keys_first,
                 RandomAccessIterator1 keys_last,
                 RandomAccessIterator2 values_first) {
    auto execPolicy = getExecutionPolicy(q, mr);

#if defined(ALPAKA_ACC_SYCL_ENABLED)
    oneapi::dpl::sort_by_key(execPolicy, keys_first, keys_last, values_first);
#else
    thrust::sort_by_key(execPolicy, keys_first, keys_last, values_first);
#endif
}

template <typename ForwardIt1, typename ForwardIt2, typename OutputIt,
          typename Compare>
void upper_bound(Queue &q, const memory_resource &mr, ForwardIt1 first1,
                 ForwardIt1 last1, ForwardIt2 first2, ForwardIt2 last2,
                 OutputIt d_first, Compare comp) {

    auto execPolicy = getExecutionPolicy(q, mr);
#if defined(ALPAKA_ACC_SYCL_ENABLED)
    oneapi::dpl::upper_bound(execPolicy, first1, last1, first2, last2, d_first,
                             comp);
#else
    thrust::upper_bound(execPolicy, first1, last1, first2, last2, d_first,
                        comp);
#endif
}

template <typename InputIt, typename OutputIt, typename Compare>
OutputIt unique_copy(Queue &q, const memory_resource &mr, InputIt first,
                     InputIt last, OutputIt d_first, Compare comp) {
    auto execPolicy = getExecutionPolicy(q, mr);

#if defined(ALPAKA_ACC_SYCL_ENABLED)
    return oneapi::dpl::unique_copy(execPolicy, first, last, d_first, comp);
#else
    return thrust::unique_copy(execPolicy, first, last, d_first, comp);
#endif
}

template <typename InputIterator, typename UnaryFunction>
void for_each(Queue &q, const memory_resource &mr, InputIterator first,
              InputIterator last, UnaryFunction f) {

    auto execPolicy = getExecutionPolicy(q, mr);
#if defined(ALPAKA_ACC_SYCL_ENABLED)
    oneapi::dpl::for_each(execPolicy, first, last, f);
#else
    thrust::for_each(execPolicy, first, last, f);
#endif
}

}  // namespace traccc::alpaka::details
