// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "traccc/ambiguity_resolution/device/sort_updated_tracks.hpp"

namespace traccc::cuda::kernels {

__global__ void sort_updated_tracks(
    device::sort_updated_tracks_payload payload);
}
