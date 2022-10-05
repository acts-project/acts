// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <torch/script.h>
#include <torch/torch.h>

namespace Acts {

// Edge building using FRNN and CUDA
torch::Tensor buildEdges(at::Tensor& embedFeatures, int64_t numSpacepoints,
                         int dim, float rVal, int kVal,
                         bool flipDirections = false);

// This function is kept for debugging. The main loop is not parallelized to not
// require additional dependencies. However, this can be easily achieved with
// TBB or the parallel STL if necessary
torch::Tensor buildEdgesBruteForce(at::Tensor& embedFeatures,
                                   int64_t numSpacepoints, int dim, float rVal,
                                   int kVal);

}  // namespace Acts
