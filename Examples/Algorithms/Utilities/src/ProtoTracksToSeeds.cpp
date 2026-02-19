// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Utilities/ProtoTracksToSeeds.hpp"

#include "ActsExamples/Utilities/EventDataTransforms.hpp"

#include <algorithm>
#include <cstddef>

namespace ActsExamples {

ProtoTracksToSeeds::ProtoTracksToSeeds(Config cfg, Acts::Logging::Level lvl)
    : IAlgorithm("PrototracksToSeeds", lvl), m_cfg(std::move(cfg)) {
  m_inputSpacePoints.initialize(m_cfg.inputSpacePoints);
  m_inputProtoTracks.initialize(m_cfg.inputProtoTracks);
  m_outputSeeds.initialize(m_cfg.outputSeeds);
  m_outputProtoTracks.initialize(m_cfg.outputProtoTracks);
}

ProcessCode ProtoTracksToSeeds::execute(const AlgorithmContext& ctx) const {
  const SpacePointContainer& spacePoints = m_inputSpacePoints(ctx);
  ProtoTrackContainer protoTracks = m_inputProtoTracks(ctx);

  const std::size_t nBefore = protoTracks.size();
  protoTracks.erase(
      std::remove_if(protoTracks.begin(), protoTracks.end(),
                     [](const ProtoTrack& pt) { return pt.size() < 3; }),
      protoTracks.end());

  if (protoTracks.size() < nBefore) {
    ACTS_DEBUG("Discarded " << nBefore - protoTracks.size()
                            << " proto tracks with less then 3 hits");
  }

  SeedContainer seeds;
  seeds.reserve(protoTracks.size());

  for (const ProtoTrack& pt : protoTracks) {
    protoTrackToSeed(pt, spacePoints, seeds);
  }

  m_outputSeeds(ctx, std::move(seeds));
  m_outputProtoTracks(ctx, std::move(protoTracks));

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
