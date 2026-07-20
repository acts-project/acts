/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/fitting/device/fit_payload.hpp"
#include "traccc/fitting/fitting_config.hpp"

// CUDA include(s).
#include <cuda_runtime.h>

namespace traccc::cuda {

template <typename fitter_t>
void fit_forward(const dim3& grid_size, const dim3& block_size,
                 std::size_t shared_mem_size, const cudaStream_t& stream,
                 const fitting_config& cfg, const device::fit_payload& payload,
                 const device::fit_tpayload<
                     typename fitter_t::detector_type::const_view_type,
                     typename fitter_t::bfield_type,
                     typename fitter_t::surface_type>* tpayload);

}  // namespace traccc::cuda
