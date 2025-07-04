// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <cstdint>
#include <utility>

#include <cuda_runtime_api.h>

namespace Acts::detail {

/// Function to perform the basic junction reomval algorithm on a graph with
/// scored edges: In case a node has more than one in/out edge, only the
/// edge with the largest score is kept.
/// NOTE: The function expects pointers on the device
/// NOTE: The function returns a pointer to device memory. The caller is
/// responsible for freeing the memory with cudaFreeAsync(ptr, stream).
/// TODO: Use some type of RAII type in the future
std::pair<std::int64_t *, std::size_t> junctionRemovalCuda(
    std::size_t nEdges, std::size_t nNodes, const float *scores,
    const std::int64_t *srcNodes, const std::int64_t *dstNodes,
    cudaStream_t stream);

}  // namespace Acts::detail
