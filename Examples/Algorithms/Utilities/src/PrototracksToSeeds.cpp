// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Utilities/PrototracksToSeeds.hpp"

#include "ActsExamples/Utilities/EventDataTransforms.hpp"

#include <algorithm>
#include <cstddef>

namespace ActsExamples {

PrototracksToSeeds::PrototracksToSeeds(Config cfg, Acts::Logging::Level lvl)
    : IAlgorithm("PrototracksToSeeds", lvl), m_cfg(std::move(cfg)) {
  m_inputSpacePoints.initialize(m_cfg.inputSpacePoints);
  m_inputProtoTracks.initialize(m_cfg.inputProtoTracks);
  m_outputSeeds.initialize(m_cfg.outputSeeds);
  m_outputProtoTracks.initialize(m_cfg.outputProtoTracks);
}

ProcessCode PrototracksToSeeds::execute(const AlgorithmContext& ctx) const {
  const SpacePointContainer& spacePoints = m_inputSpacePoints(ctx);
  ProtoTrackContainer prototracks = m_inputProtoTracks(ctx);

  const std::size_t nBefore = prototracks.size();
  prototracks.erase(
      std::remove_if(prototracks.begin(), prototracks.end(),
                     [](const ProtoTrack& pt) { return pt.size() < 3; }),
      prototracks.end());

  if (prototracks.size() < nBefore) {
    ACTS_DEBUG("Discarded " << nBefore - prototracks.size()
                            << " prototracks with less then 3 hits");
  }

  SeedContainer seeds;
  seeds.reserve(prototracks.size());

  for (const ProtoTrack& pt : prototracks) {
    prototrackToSeed(pt, spacePoints, seeds);
  }

  m_outputSeeds(ctx, std::move(seeds));
  m_outputProtoTracks(ctx, std::move(prototracks));

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
