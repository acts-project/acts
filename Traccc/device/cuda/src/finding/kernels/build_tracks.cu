/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "../../utils/global_index.hpp"
#include "build_tracks.cuh"

// Project include(s).
#include "traccc/edm/track_parameters.hpp"
#include "traccc/finding/candidate_link.hpp"
#include "traccc/finding/device/build_tracks.hpp"
#include "traccc/finding/finding_config.hpp"

namespace traccc::cuda::kernels {

__global__ void build_tracks(
    const __grid_constant__ bool run_mbf,
    const __grid_constant__ measurement_selector::config calib_cfg,
    const __grid_constant__ device::build_tracks_payload payload) {

    device::build_tracks(details::global_index1(), run_mbf, calib_cfg, payload);
}
}  // namespace traccc::cuda::kernels
