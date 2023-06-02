// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Utilities/PrototracksToSeeds.hpp"

#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsExamples/EventData/SimSeed.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/Utilities/EventDataTransforms.hpp"

#include <algorithm>

namespace ActsExamples {

PrototracksToSeeds::PrototracksToSeeds(Config cfg, Acts::Logging::Level lvl)
    : IAlgorithm("PrototracksToSeeds", lvl), m_cfg(std::move(cfg)) {
  m_outputSeeds.initialize(m_cfg.outputSeeds);
  m_outputProtoTracks.initialize(m_cfg.outputProtoTracks);
  m_inputProtoTracks.initialize(m_cfg.inputProtoTracks);
  m_inputSpacePoints.initialize(m_cfg.inputSpacePoints);
}

ProcessCode PrototracksToSeeds::execute(const AlgorithmContext& ctx) const {
  auto prototracks = m_inputProtoTracks(ctx);

  prototracks.erase(std::remove_if(prototracks.begin(), prototracks.end(),
                                   [](const auto& t) { return t.size() < 3; }),
                    prototracks.end());

  SimSeedContainer seeds;
  seeds.reserve(prototracks.size());

  const ProtoTrackToSeed prototrackToSeed{m_inputSpacePoints(ctx)};
  std::transform(prototracks.begin(), prototracks.end(),
                 std::back_inserter(seeds), prototrackToSeed);

  m_outputSeeds(ctx, std::move(seeds));
  m_outputProtoTracks(ctx, std::move(prototracks));

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
