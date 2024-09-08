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

namespace Acts::detail {

GraphCreatorWrapperCuda::GraphCreatorWrapperCuda(const std::string &path,
                                                 int device, int blocks) {
  m_graphCreator = std::make_unique<CUDA_graph_creator<float>>(
      blocks, device, path, 10,
      std::pair<float, float>{0.f, std::numeric_limits<float>::max()});
}

GraphCreatorWrapperCuda::~GraphCreatorWrapperCuda() {}

graph<float> GraphCreatorWrapperCuda::build(TTree_hits<float> &hits) {
  CUDA_graph_creator<float>::graph_building_stats stats;
  return m_graphCreator->build_impl(hits, stats, false);
}

}  // namespace Acts::detail
