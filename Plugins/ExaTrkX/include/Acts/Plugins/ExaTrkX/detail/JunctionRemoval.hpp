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

std::pair<std::int64_t *, std::size_t> junctionRemovalCuda(
    std::size_t nEdges, std::size_t nNodes, const float *scores,
    const std::int64_t *srcNodes, const std::int64_t *dstNodes,
    cudaStream_t stream);

}
