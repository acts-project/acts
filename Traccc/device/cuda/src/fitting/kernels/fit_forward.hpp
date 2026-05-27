/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include <cuda_runtime.h>

#include "traccc/fitting/device/fit.hpp"
#include "traccc/fitting/fitting_config.hpp"

namespace traccc::cuda {

template <typename fitter_t>
void fit_forward(const dim3& grid_size, const dim3& block_size,
                 std::size_t shared_mem_size, const cudaStream_t& stream,
                 const fitting_config& cfg,
                 const device::fit_payload<fitter_t>& payload);

}  // namespace traccc::cuda
