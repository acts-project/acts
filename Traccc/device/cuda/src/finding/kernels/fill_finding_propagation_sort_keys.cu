// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Local include(s).
#include "../../utils/global_index.hpp"
#include "fill_finding_propagation_sort_keys.cuh"

// Project include(s).
#include "traccc/finding/device/fill_finding_propagation_sort_keys.hpp"

namespace traccc::cuda::kernels {

__global__ void fill_finding_propagation_sort_keys(
    const __grid_constant__ device::fill_finding_propagation_sort_keys_payload
        payload) {
  device::fill_finding_propagation_sort_keys(details::global_index1(), payload);
}

}  // namespace traccc::cuda::kernels
