/**
 * traccc library, part of the ACTS project (R&D line)
 *
 * (c) 2024-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "../utils/cuda_error_handling.hpp"
#include "../utils/global_index.hpp"
#include "../utils/utils.hpp"
#include "traccc/cuda/utils/stream.hpp"

// VecMem include(s).
#include <vecmem/memory/memory_resource.hpp>
#include <vecmem/memory/unique_ptr.hpp>
#include <vecmem/utils/copy.hpp>

// CUDA include
#include <cuda_runtime.h>

// System include
#include <concepts>
#include <utility>

namespace traccc::cuda {
namespace kernels {
template <typename CONTAINER, std::semiregular R, typename VIEW>
    requires std::regular_invocable<R,
                                    decltype(std::declval<CONTAINER>().at(0)),
                                    decltype(std::declval<CONTAINER>().at(0))>
__global__ void is_ordered_on_kernel(R relation, VIEW _in, bool* out) {

    const device::global_index_t tid = details::global_index1();

    const CONTAINER in(_in);

    if (tid > 0 && tid < in.size()) {
        if (!relation(in.at(tid - 1), in.at(tid))) {
            *out = false;
        }
    }
}
}  // namespace kernels

/**
 * @brief Sanity check that a given container is ordered on a given relation.
 *
 * For a container $v$ to be ordered on a relation $R$, it must be the case that
 * for all indices $i$ and $j$, if $i < j$, then $R(i, j)$.
 *
 * @note This function runs in O(n) time.
 *
 * @note Although functions like `std::sort` requires the relation to be strict
 * weak order, this function is more lax in its requirements. Rather, the
 * relation should be a total preorder, i.e. a non-strict weak order.
 *
 * @note For any strict weak order $R$, `is_ordered_on(sort(R, v))` is true.
 *
 * @tparam CONTAINER The type of the (device) container.
 * @tparam R The type of relation $R$, a callable which returns a bool if the
 * first argument can be immediately before the second type.
 * @tparam VIEW The type of the view for the container.
 * @param relation A relation object of type `R`.
 * @param mr A memory resource used for allocating intermediate memory.
 * @param view The container which to check for ordering.
 * @return true If the container is ordered on `R`.
 * @return false Otherwise.
 */
template <typename CONTAINER, std::semiregular R, typename VIEW>
    requires std::regular_invocable<R,
                                    decltype(std::declval<CONTAINER>().at(0)),
                                    decltype(std::declval<CONTAINER>().at(0))>
bool is_ordered_on(R&& relation, vecmem::memory_resource& mr,
                   const vecmem::copy& copy, stream& stream, const VIEW& view) {

    // This should never be a performance-critical step, so we can keep the
    // block size fixed.
    constexpr int block_size = 512;

    cudaStream_t cuda_stream = details::get_stream(stream);

    // Grab the number of elements in our container.
    const typename VIEW::size_type n = copy.get_size(view);

    // Exit early for empty containers.
    if (n == 0) {
        return true;
    }

    // Initialize the output boolean.
    vecmem::unique_alloc_ptr<bool> out = vecmem::make_unique_alloc<bool>(mr);
    bool initial_out = true;
    TRACCC_CUDA_ERROR_CHECK(
        cudaMemcpyAsync(out.get(), &initial_out, sizeof(bool),
                        cudaMemcpyHostToDevice, cuda_stream));

    // Launch the kernel which will write its result to the `out` boolean.
    kernels::is_ordered_on_kernel<CONTAINER>
        <<<(n + block_size - 1) / block_size, block_size, 0, cuda_stream>>>(
            relation, view, out.get());

    TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());

    // Copy the output to host, then return it.
    bool host_out;

    TRACCC_CUDA_ERROR_CHECK(cudaMemcpyAsync(&host_out, out.get(), sizeof(bool),
                                            cudaMemcpyDeviceToHost,
                                            cuda_stream));

    stream.synchronize();

    return host_out;
}
}  // namespace traccc::cuda
