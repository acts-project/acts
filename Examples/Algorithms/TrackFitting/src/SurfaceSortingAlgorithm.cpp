// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFitting/SurfaceSortingAlgorithm.hpp"

#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsFatras/EventData/Hit.hpp"

#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

namespace ActsExamples {

SurfaceSortingAlgorithm::SurfaceSortingAlgorithm(Config cfg,
                                                 Acts::Logging::Level level)
    : IAlgorithm("SurfaceSortingAlgorithm", level), m_cfg(std::move(cfg)) {
  if (m_cfg.inputProtoTracks.empty()) {
    throw std::invalid_argument("Missing input proto track collection");
  }
  if (m_cfg.inputSimHits.empty()) {
    throw std::invalid_argument("Missing input simulated hits collection");
  }
  if (m_cfg.inputMeasurementSimHitsMap.empty()) {
    throw std::invalid_argument("Missing input measurement sim hits map");
  }
  if (m_cfg.outputProtoTracks.empty()) {
    throw std::invalid_argument("Missing output proto track collection");
  }

  m_inputProtoTracks.initialize(m_cfg.inputProtoTracks);
  m_inputSimHits.initialize(m_cfg.inputSimHits);
  m_inputMeasurementSimHitsMap.initialize(m_cfg.inputMeasurementSimHitsMap);
  m_outputProtoTracks.initialize(m_cfg.outputProtoTracks);
}

ProcessCode SurfaceSortingAlgorithm::execute(
    const AlgorithmContext& ctx) const {
  const auto& protoTracks = m_inputProtoTracks(ctx);
  const auto& simHits = m_inputSimHits(ctx);
  const auto& simHitsMap = m_inputMeasurementSimHitsMap(ctx);

  ProtoTrackContainer sortedTracks;
  sortedTracks.reserve(protoTracks.size());
  std::map<double, Index> trackHitList;

  for (std::size_t itrack = 0; itrack < protoTracks.size(); ++itrack) {
    const auto& protoTrack = protoTracks[itrack];

    ProtoTrack sortedProtoTrack;
    sortedProtoTrack.reserve(protoTrack.size());
    trackHitList.clear();

    if (protoTrack.empty()) {
      continue;
    }

    for (const auto hit : protoTrack) {
      const auto simHitIndex = simHitsMap.find(hit)->second;
      auto simHit = simHits.nth(simHitIndex);
      auto simHitTime = simHit->time();
      trackHitList.insert(std::make_pair(simHitTime, hit));
    }

    /// Map will now be sorted by truth hit time
    for (auto const& [time, hit] : trackHitList) {
      sortedProtoTrack.emplace_back(hit);
    }

    sortedTracks.emplace_back(std::move(sortedProtoTrack));
  }

  m_outputProtoTracks(ctx, std::move(sortedTracks));

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
