/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "../../../utils/global_index.hpp"
#include "../fit_backward.hpp"

// Project include(s).
#include "traccc/fitting/device/fit_backward.hpp"

namespace traccc::cuda {
namespace kernels {

template <typename fitter_t>
__global__ __launch_bounds__(128) void fit_backward(
    const fitting_config cfg, const device::fit_payload payload,
    const device::fit_tpayload<
        typename fitter_t::detector_type::const_view_type,
        typename fitter_t::bfield_type, typename fitter_t::surface_type>*
        tpayload) {

    device::fit_backward<fitter_t>(details::global_index1(), cfg, payload,
                                   *tpayload);
}

}  // namespace kernels

template <typename fitter_t>
void fit_backward(const dim3& grid_size, const dim3& block_size,
                  std::size_t shared_mem_size, const cudaStream_t& stream,
                  const fitting_config& cfg, const device::fit_payload& payload,
                  const device::fit_tpayload<
                      typename fitter_t::detector_type::const_view_type,
                      typename fitter_t::bfield_type,
                      typename fitter_t::surface_type>* tpayload) {

    kernels::fit_backward<fitter_t>
        <<<grid_size, block_size, shared_mem_size, stream>>>(cfg, payload,
                                                             tpayload);
}

}  // namespace traccc::cuda
