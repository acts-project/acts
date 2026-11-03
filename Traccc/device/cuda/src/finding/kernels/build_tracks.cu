// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

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
