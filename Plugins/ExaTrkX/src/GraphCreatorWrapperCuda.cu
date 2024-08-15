// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/ExaTrkX/detail/GraphCreaterWrapper.hpp"

#include <TTree_hits>
#include <graph>
#include <CUDA_graph_creator>

namespace Acts::detail {

GraphCreaterWrapperCuda::GraphCreaterWrapperCuda(const std::string &path, int device) {
  m_graphCreator = std::make_unique<CUDA_graph_creator<float>>(
      device, path, 10,
      std::pair<float, float>{0.f, std::numeric_limits<float>::max()});
}

graph<float> GraphCreaterWrapperCuda::build(TTree_hits<float> &hits) {
  CUDA_graph_creator<float>::graph_building_stats stats;ś
  return m_graphCreator->build_impl(hits, stats);
}

}  // namespace Acts::detail

