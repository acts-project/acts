// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Traccc/ActsSpToTracccAlg.hpp"

#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"

#include <stdexcept>

namespace ActsExamples {

ActsSpToTracccAlg::ActsSpToTracccAlg(const Config& cfg,
                                     std::unique_ptr<const Acts::Logger> logger)
    : IAlgorithm("ActsSpToTracccAlg", std::move(logger)), m_cfg(cfg) {
  m_inputSpacePoints.initialize(m_cfg.inputSpacePoints);
  m_outputTracccSpacepoints.initialize(m_cfg.outputTracccSpacepoints);
}

ProcessCode ActsSpToTracccAlg::execute(const AlgorithmContext& ctx) const {
  const auto& actsSpacePoints = m_inputSpacePoints(ctx);

  traccc::edm::spacepoint_collection::host tracccSpacepoints{m_mr};
  std::size_t nConverted = 0;
  std::size_t nSkipped = 0;

  for (const auto& sp : actsSpacePoints) {
    const auto sourceLinks = sp.sourceLinks();
    if (sourceLinks.empty()) {
      ++nSkipped;
      continue;
    }

    // Get first source link — pixel or first strip of pair
    const auto& sl1 = sourceLinks[0].get<IndexSourceLink>();
    const std::size_t tracccIdx1 = sl1.index();

    tracccSpacepoints.push_back({});
    auto tsp = tracccSpacepoints.at(tracccSpacepoints.size() - 1);
    tsp.measurement_index_1() = tracccIdx1;
    tsp.global()[0] = sp.x();
    tsp.global()[1] = sp.y();
    tsp.global()[2] = sp.z();
    tsp.z_variance() = sp.varianceZ();
    tsp.radius_variance() = sp.varianceR();

    ++nConverted;
  }

  ACTS_INFO("ActsSpToTracccAlg: converted "
            << nConverted << " / " << actsSpacePoints.size() << " spacepoints"
            << " (skipped " << nSkipped << ")");

  m_outputTracccSpacepoints(ctx, std::move(tracccSpacepoints));
  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
