// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "traccc/definitions/primitives.hpp"
#include "traccc/finding/device/gather_best_tips_per_measurement.hpp"

namespace traccc::cuda::kernels {

__global__ void gather_best_tips_per_measurement(
    const __grid_constant__
        device::gather_best_tips_per_measurement_payload<default_algebra>
            payload);

}  // namespace traccc::cuda::kernels
