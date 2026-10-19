// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "traccc/finding/device/build_tracks.hpp"
#include "traccc/finding/finding_config.hpp"

namespace traccc::cuda::kernels {

__global__ void build_tracks(
    const __grid_constant__ bool run_mbf,
    const __grid_constant__ measurement_selector::config calib_cfg,
    const __grid_constant__ device::build_tracks_payload payload);
}
