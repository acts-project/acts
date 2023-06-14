// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <cstdint>

namespace at {
class Tensor;
}

namespace Acts {

/// Edge building using FRNN and CUDA. If CUDA is not available,
/// a brute-force CPU method is used
/// TODO implement better CPU method (K-D-Tree, ...)
/// TODO make parameters concise (the tensor should have the numSpacepoints, dim
/// info)
///
/// @param embedFeatures Tensor of shape (n_nodes, embedding_dim)
/// @param numSpacepoints number of spacepoints
/// @param dim embedding embedding dim
/// @param rVal radius for NN search
/// @param kVal max number of neighbours in NN search
/// @param flipDirections if we want to randomly flip directions of the
/// edges after the edge building
at::Tensor buildEdges(at::Tensor& embedFeatures, int64_t numSpacepoints,
                      int dim, float rVal, int kVal,
                      bool flipDirections = false);

}  // namespace Acts
