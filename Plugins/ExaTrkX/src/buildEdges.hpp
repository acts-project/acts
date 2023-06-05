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

}  // namespace Acts
