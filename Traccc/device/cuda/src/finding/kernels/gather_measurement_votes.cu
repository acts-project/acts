// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Local include(s).
#include "../../utils/global_index.hpp"
#include "gather_measurement_votes.cuh"

namespace traccc::cuda::kernels {

__global__ void gather_measurement_votes(
    const __grid_constant__ device::gather_measurement_votes_payload payload) {
  device::gather_measurement_votes(details::global_index1(), payload);
}

}  // namespace traccc::cuda::kernels
