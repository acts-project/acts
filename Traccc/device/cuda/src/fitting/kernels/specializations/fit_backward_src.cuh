/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include "../../../utils/global_index.hpp"
#include "../fit_backward.hpp"
#include "traccc/fitting/device/fit_backward.hpp"

namespace traccc::cuda {
namespace kernels {
template <typename fitter_t>
__global__ __launch_bounds__(128) void fit_backward(
    const fitting_config cfg, const device::fit_payload<fitter_t> payload) {
    device::fit_backward<fitter_t>(details::global_index1(), cfg, payload);
}
}  // namespace kernels

template <typename fitter_t>
void fit_backward(const dim3& grid_size, const dim3& block_size,
                  std::size_t shared_mem_size, const cudaStream_t& stream,
                  const fitting_config& cfg,
                  const device::fit_payload<fitter_t>& payload) {
    kernels::fit_backward<fitter_t>
        <<<grid_size, block_size, shared_mem_size, stream>>>(cfg, payload);
}

}  // namespace traccc::cuda
