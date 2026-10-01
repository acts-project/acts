// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "traccc/edm/device/sort_key.hpp"
#include "traccc/edm/track_collection.hpp"

// CUDA include(s).
#include <cuda_runtime.h>

// VecMem include(s).
#include <vecmem/containers/data/vector_view.hpp>

namespace traccc::cuda {

/// Function calling a kernel for @c traccc::device::fill_fitting_sort_keys
void fill_fitting_sort_keys(
    const dim3& grid_size, const dim3& block_size, cudaStream_t stream,
    edm::track_collection<default_algebra>::const_view track_candidates_view,
    vecmem::data::vector_view<device::sort_key> keys_view,
    vecmem::data::vector_view<unsigned int> ids_view);

}  // namespace traccc::cuda
