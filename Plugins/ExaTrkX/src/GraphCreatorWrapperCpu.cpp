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
#include <graph_creator>

namespace Acts::detail {

GraphCreaterWrapperCpu::GraphCreaterWrapperCpu(const std::string &path) {
  m_graphCreator = std::make_unique<graph_creator<float>>(
      path, 10,
      std::pair<float, float>{0.f, std::numeric_limits<float>::max()});
}

graph<float> GraphCreaterWrapperCpu::build(TTree_hits<float> &hits) {
  TTree_particles<float> particlesTree;
  std::string eventId = "no-id";

  auto [graph, _] =
      m_graphCreator->build_impl(hitsTree, particlesTree, eventId, print);

  return graph;
}

}  // namespace Acts::detail
