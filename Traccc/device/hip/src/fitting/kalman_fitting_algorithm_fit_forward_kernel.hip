/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "../utils/global_index.hpp"
#include "../utils/hip_error_handling.hpp"
#include "../utils/magnetic_field_types.hpp"
#include "../utils/utils.hpp"
#include "traccc/hip/fitting/kalman_fitting_algorithm.hpp"

// Project include(s).
#include "traccc/fitting/details/kalman_fitting_types.hpp"
#include "traccc/fitting/device/fit_forward.hpp"
#include "traccc/geometry/detector.hpp"
#include "traccc/utils/detector_buffer_bfield_visitor.hpp"

// HIP include(s).
#include <hip/hip_runtime.h>

namespace traccc::hip {
namespace kernels {

template <typename fitter_t>
__global__ __launch_bounds__(128) void fit_forward(
    const fitting_config cfg, const device::fit_payload payload,
    const device::fit_tpayload<
        typename fitter_t::detector_type::const_view_type,
        typename fitter_t::bfield_type, typename fitter_t::surface_type>*
        tpayload) {
  device::fit_forward<fitter_t>(details::global_index1(), cfg, payload,
                                *tpayload);
}

}  // namespace kernels

void kalman_fitting_algorithm::fit_forward_kernel(
    const fitting_config& config, const fit_payload& payload) const {
  return detector_buffer_magnetic_field_visitor<detector_type_list,
                                                hip::bfield_type_list<scalar>>(
      payload.detector, payload.field,
      [&]<typename detector_traits_t, typename bfield_view_t>(
          const typename detector_traits_t::view&, const bfield_view_t&) {
        // Get the number of tracks.
        const unsigned int n_tracks = payload.payload.tracks.tracks.capacity();
        assert(n_tracks == copy().get_size(payload.payload.tracks.tracks));

        // Launch parameters for the kernel.
        const unsigned int nThreads = warp_size() * 4;
        const unsigned int nBlocks = (n_tracks + nThreads - 1) / nThreads;

        // Fitter type to use.
        using fitter_t =
            traccc::details::kalman_fitter_t<typename detector_traits_t::device,
                                             bfield_view_t>;

        // Run the track fitting
        hipLaunchKernelGGL(kernels::fit_forward<fitter_t>, nBlocks, nThreads, 0,
                           details::get_stream(stream()), config,
                           payload.payload, payload.get_tpayload<fitter_t>());
        TRACCC_HIP_ERROR_CHECK(hipGetLastError());
      });
}

}  // namespace traccc::hip
