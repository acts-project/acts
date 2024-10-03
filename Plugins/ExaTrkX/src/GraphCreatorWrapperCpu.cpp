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

#include <TTree_hits>
#include <graph>
#include <graph_creator>

namespace Acts::detail {

GraphCreatorWrapperCpu::GraphCreatorWrapperCpu(const std::string &path) {
  m_graphCreator = std::make_unique<graph_creator<float>>(
      path, 10,
      std::pair<float, float>{0.f, std::numeric_limits<float>::max()});
}

GraphCreatorWrapperCpu::~GraphCreatorWrapperCpu() {}

graph<float> GraphCreatorWrapperCpu::build(TTree_hits<float> &hits,
                                           bool print) {
  TTree_particles<float> particles;
  std::string eventId = "no-id";

  auto [graph, _] = m_graphCreator->build_impl(hits, particles, eventId, print);

  return graph;
}

}  // namespace Acts::detail
