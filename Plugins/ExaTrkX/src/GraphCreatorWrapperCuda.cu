// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/ExaTrkX/detail/GraphCreatorWrapper.hpp"

#include <CUDA_graph_creator>
#include <TTree_hits>
#include <graph>

#include "oldApiHelper.hpp"

namespace Acts::detail {

GraphCreatorWrapperCuda::GraphCreatorWrapperCuda(const std::string &path,
                                                 int device, int blocks) {
  m_graphCreator = std::make_unique<CUDA_graph_creator<float>>(
      blocks, device, path, 10,
      std::pair<float, float>{0.f, std::numeric_limits<float>::max()});
}

GraphCreatorWrapperCuda::~GraphCreatorWrapperCuda() {}

std::pair<at::Tensor, at::Tensor> GraphCreatorWrapperCuda::build(
    const std::vector<float> &features,
    const std::vector<std::uint64_t> &moduleIds, const Acts::Logger &logger) {
  auto builder = [&](auto &hits, bool print) {
    CUDA_graph_creator<float>::graph_building_stats stats;
    return m_graphCreator->build_impl(hits, stats, print);
  };
  return oldApiBuild(features, moduleIds, logger, builder, 1000.f, 3.14159f,
                     1000.f);
}

}  // namespace Acts::detail
