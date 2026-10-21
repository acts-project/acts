// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "../../utils/global_index.hpp"
#include "remove_duplicates.cuh"
#include "traccc/finding/device/remove_duplicates.hpp"

namespace traccc::cuda::kernels {

__global__ void remove_duplicates(
    const __grid_constant__ finding_config cfg,
    const __grid_constant__ device::remove_duplicates_payload payload) {
  device::remove_duplicates(details::global_index1(), cfg, payload);
}

}  // namespace traccc::cuda::kernels
