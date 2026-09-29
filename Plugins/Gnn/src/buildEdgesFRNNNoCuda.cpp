// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/Gnn/detail/buildEdges.hpp"

#include <torch/torch.h>

#include "CudaUnavailable.hpp"

torch::Tensor ActsPlugins::detail::buildEdgesFRNN(
    torch::Tensor & /*embedFeatures*/, float /*rVal*/, int /*kVal*/,
    bool /*flipDirections*/) {
  throwNoCudaSupport("run FRNN edge building");
}
