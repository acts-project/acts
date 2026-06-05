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
#include <vecmem/containers/data/vector_buffer.hpp>
#include <vecmem/containers/device_vector.hpp>
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

/// Kernel used in implementing @c traccc::cuda::is_contiguous_on
template <typename CONTAINER, std::semiregular P, typename VIEW,
          std::equality_comparable S>
    requires std::regular_invocable<P,
                                    decltype(std::declval<CONTAINER>().at(0))>
__global__ void is_contiguous_on_compress_adjacent(
    P projection, VIEW _in, vecmem::data::vector_view<S> out_view) {

    const device::global_index_t tid = details::global_index1();

    const CONTAINER in(_in);
    vecmem::device_vector<S> out(out_view);

    if (tid > 0 && tid < in.size()) {
        S v1 = projection(in.at(tid - 1));
        S v2 = projection(in.at(tid));

        if (v1 != v2) {
            out.push_back(v2);
        }
    } else if (tid == 0) {
        out.push_back(projection(in.at(tid)));
    }
}

/// Kernel used in implementing @c traccc::cuda::is_contiguous_on
template <std::equality_comparable T>
__global__ void is_contiguous_on_all_unique(
    vecmem::data::vector_view<T> in_view, bool* out) {

    const device::global_index_t tid_x = threadIdx.x + blockIdx.x * blockDim.x;
    const device::global_index_t tid_y = threadIdx.y + blockIdx.y * blockDim.y;

    const vecmem::device_vector<T> in(in_view);

    if (tid_x < in.size() && tid_y < in.size() && tid_x != tid_y &&
        in.at(tid_x) == in.at(tid_y)) {
        *out = false;
    }
}
}  // namespace kernels

/**
 * @brief Sanity check that a given container is contiguous on a given
 *        projection.
 *
 * For a container $v$ to be contiguous on a projection $\pi$, it must be the
 * case that for all indices $i$ and $j$, if $v_i = v_j$, then all indices $k$
 * between $i$ and $j$, $v_i = v_j = v_k$.
 *
 * @note This function runs in O(n^2) time.
 *
 * @tparam CONTAINER The type of the (device) container.
 * @tparam P The type of projection $\pi$, a callable which returns some
 * comparable type.
 * @tparam VIEW The type of the view for the container.
 * @param projection A projection object of type `P`.
 * @param mr A memory resource used for allocating intermediate memory.
 * @param view The container which to check for contiguity.
 * @return true If the container is contiguous on `P`.
 * @return false Otherwise.
 */
template <typename CONTAINER, std::semiregular P, typename VIEW>
    requires std::regular_invocable<P,
                                    decltype(std::declval<CONTAINER>().at(0))>
bool is_contiguous_on(P&& projection, vecmem::memory_resource& mr,
                      const vecmem::copy& copy, stream& stream,
                      const VIEW& view) {

    // This should never be a performance-critical step, so we can keep the
    // block size fixed.
    constexpr int block_size = 512;
    constexpr int block_size_2d = 32;

    cudaStream_t cuda_stream = details::get_stream(stream);

    // Grab the number of elements in our container.
    const typename VIEW::size_type n = copy.get_size(view);

    // Exit early for empty containers.
    if (n == 0) {
        return true;
    }

    // Get the output type of the projection.
    using projection_t =
        std::invoke_result_t<P, decltype(std::declval<CONTAINER>().at(0))>;

    // Allocate memory for intermediate values and outputs, then set them up.
    vecmem::data::vector_buffer<projection_t> iout(
        n, mr, vecmem::data::buffer_type::resizable);
    copy.setup(iout)->ignore();
    vecmem::unique_alloc_ptr<bool> out = vecmem::make_unique_alloc<bool>(mr);

    bool initial_out = true;

    TRACCC_CUDA_ERROR_CHECK(
        cudaMemcpyAsync(out.get(), &initial_out, sizeof(bool),
                        cudaMemcpyHostToDevice, cuda_stream));

    // Launch the first kernel, which will squash consecutive equal elements
    // into one element.
    kernels::is_contiguous_on_compress_adjacent<CONTAINER>
        <<<(n + block_size - 1) / block_size, block_size, 0, cuda_stream>>>(
            projection, view, iout);

    TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());

    // Launch the second kernel, which will check if the values are unique.
    uint32_t grid_size_rd =
        (copy.get_size(iout) + block_size_2d - 1) / block_size_2d;
    dim3 all_unique_grid_size(grid_size_rd, grid_size_rd);
    dim3 all_unique_block_size(block_size_2d, block_size_2d);

    kernels::is_contiguous_on_all_unique<<<
        all_unique_grid_size, all_unique_block_size, 0, cuda_stream>>>(
        iout, out.get());

    TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());

    // Get the result from the device and return it.
    bool host_out;

    TRACCC_CUDA_ERROR_CHECK(cudaMemcpyAsync(&host_out, out.get(), sizeof(bool),
                                            cudaMemcpyDeviceToHost,
                                            cuda_stream));

    stream.synchronize();

    return host_out;
}
}  // namespace traccc::cuda
