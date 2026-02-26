// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Utilities/ProtoTracksToSeeds.hpp"

#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsExamples/EventData/SimSeed.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/Utilities/EventDataTransforms.hpp"

#include <algorithm>

namespace ActsExamples {

ProtoTracksToSeeds::ProtoTracksToSeeds(
    Config cfg, std::unique_ptr<const Acts::Logger> logger)
    : IAlgorithm("ProtoTracksToSeeds", std::move(logger)),
      m_cfg(std::move(cfg)) {
  m_outputSeeds.initialize(m_cfg.outputSeeds);
  m_outputProtoTracks.initialize(m_cfg.outputProtoTracks);
  m_inputProtoTracks.initialize(m_cfg.inputProtoTracks);
  m_inputSpacePoints.initialize(m_cfg.inputSpacePoints);
}

ProcessCode ProtoTracksToSeeds::execute(const AlgorithmContext& ctx) const {
  auto protoTracks = m_inputProtoTracks(ctx);

  const auto nBefore = protoTracks.size();
  protoTracks.erase(std::remove_if(protoTracks.begin(), protoTracks.end(),
                                   [](const auto& t) { return t.size() < 3; }),
                    protoTracks.end());
  ACTS_DEBUG("Discarded " << protoTracks.size() - nBefore
                          << " proto tracks with less then 3 hits");

  SimSeedContainer seeds;
  seeds.reserve(protoTracks.size());

  const auto& sps = m_inputSpacePoints(ctx);
  std::transform(protoTracks.begin(), protoTracks.end(),
                 std::back_inserter(seeds),
                 [&](const auto& pt) { return protoTrackToSeed(pt, sps); });

  m_outputSeeds(ctx, std::move(seeds));
  m_outputProtoTracks(ctx, std::move(protoTracks));

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
