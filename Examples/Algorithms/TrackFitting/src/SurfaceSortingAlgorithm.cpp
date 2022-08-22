// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFitting/SurfaceSortingAlgorithm.hpp"

#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"

ActsExamples::SurfaceSortingAlgorithm::SurfaceSortingAlgorithm(
    Config cfg, Acts::Logging::Level level)
    : ActsExamples::BareAlgorithm("SurfaceSortingAlgorithm", level),
      m_cfg(std::move(cfg)) {
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
}

ActsExamples::ProcessCode ActsExamples::SurfaceSortingAlgorithm::execute(
    const ActsExamples::AlgorithmContext& ctx) const {
  using HitSimHitsMap = IndexMultimap<Index>;

  const auto& protoTracks =
      ctx.eventStore.get<ProtoTrackContainer>(m_cfg.inputProtoTracks);
  const auto& simHits = ctx.eventStore.get<SimHitContainer>(m_cfg.inputSimHits);
  const auto& simHitsMap =
      ctx.eventStore.get<HitSimHitsMap>(m_cfg.inputMeasurementSimHitsMap);

  ProtoTrackContainer sortedTracks;
  sortedTracks.reserve(protoTracks.size());
  TrackHitList trackHitList;

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

  ctx.eventStore.add(m_cfg.outputProtoTracks, std::move(sortedTracks));

  return ActsExamples::ProcessCode::SUCCESS;
}
