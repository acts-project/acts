// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Local include(s).
#include "../../utils/barrier.hpp"
#include "../../utils/global_index.hpp"
#include "gather_best_tips_per_measurement.cuh"

namespace traccc::cuda::kernels {

__global__ void gather_best_tips_per_measurement(
    const __grid_constant__
        device::gather_best_tips_per_measurement_payload<default_algebra>
            payload) {
  device::gather_best_tips_per_measurement(details::global_index1(), barrier{},
                                           payload);
}

}  // namespace traccc::cuda::kernels
