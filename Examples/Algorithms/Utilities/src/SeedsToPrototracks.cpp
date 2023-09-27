// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Utilities/SeedsToPrototracks.hpp"

#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsExamples/EventData/SimSeed.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/Utilities/EventDataTransforms.hpp"

namespace ActsExamples {

SeedsToPrototracks::SeedsToPrototracks(Config cfg, Acts::Logging::Level lvl)
    : IAlgorithm("TrajectoriesToPrototracks", lvl), m_cfg(std::move(cfg)) {
  m_inputSeeds.initialize(m_cfg.inputSeeds);
  m_outputProtoTracks.initialize(m_cfg.outputProtoTracks);
}

ProcessCode SeedsToPrototracks::execute(const AlgorithmContext& ctx) const {
  const auto seeds = m_inputSeeds(ctx);

  auto tracks = seedsToPrototracks(seeds);

  m_outputProtoTracks(ctx, std::move(tracks));

  return ProcessCode::SUCCESS;
}
}  // namespace ActsExamples
