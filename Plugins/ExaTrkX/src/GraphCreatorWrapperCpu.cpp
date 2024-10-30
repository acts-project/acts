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

#include "oldApiHelper.hpp"

namespace Acts::detail {

GraphCreatorWrapperCpu::GraphCreatorWrapperCpu(const std::string &path) {
#ifdef OLD_MMG_API
  m_graphCreator = std::make_unique<graph_creator<float>>(
      path, 10,
      std::pair<float, float>{0.f, std::numeric_limits<float>::max()});
#endif
}

GraphCreatorWrapperCpu::~GraphCreatorWrapperCpu() {}

std::pair<at::Tensor, at::Tensor> GraphCreatorWrapperCpu::build(
    const std::vector<float> &features,
    const std::vector<std::uint64_t> &moduleIds, const Acts::Logger &logger) {
#ifdef OLD_MMG_API
  auto builder = [&](auto &hits, bool print) {
    TTree_particles<float> particles;
    std::string eventId = "no-id";

    auto [graph, _] =
        m_graphCreator->build_impl(hits, particles, eventId, print);

    return graph;
  };

  // TODO hardocded
  return oldApiBuild(features, moduleIds, logger, builder, 1000.f, 3.14159f,
                     1000.f);
#else
  return {};
#endif
}

}  // namespace Acts::detail
