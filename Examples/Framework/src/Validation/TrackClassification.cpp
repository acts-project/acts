// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Validation/TrackClassification.hpp"

#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/EventData/Trajectories.hpp"
#include "ActsExamples/Utilities/Range.hpp"

#include <algorithm>

namespace ActsExamples {

namespace {

/// Increase the hit count for the given particle id by one.
inline void increaseHitCount(std::vector<ParticleHitCount>& particleHitCounts,
                             ActsFatras::Barcode particleId) {
  // linear search since there is no ordering
  auto it = std::ranges::find_if(particleHitCounts, [=](const auto& phc) {
    return (phc.particleId == particleId);
  });
  // either increase count if we saw the particle before or add it
  if (it != particleHitCounts.end()) {
    it->hitCount += 1u;
  } else {
    particleHitCounts.push_back({particleId, 1u});
  }
}

/// Sort hit counts by decreasing values, i.e. majority particle comes first.
inline void sortHitCount(std::vector<ParticleHitCount>& particleHitCounts) {
  std::ranges::sort(particleHitCounts, std::greater{},
                    [](const auto& p) { return p.hitCount; });
}

}  // namespace

void identifyContributingParticles(
    const IndexMultimap<ActsFatras::Barcode>& hitParticlesMap,
    const ProtoTrack& protoTrack,
    std::vector<ParticleHitCount>& particleHitCounts) {
  particleHitCounts.clear();

  for (auto hitIndex : protoTrack) {
    // register all particles that generated this hit
    for (const auto& [_, value] :
         makeRange(hitParticlesMap.equal_range(hitIndex))) {
      increaseHitCount(particleHitCounts, value);
    }
  }
  sortHitCount(particleHitCounts);
}

void identifyContributingParticles(
    const IndexMultimap<ActsFatras::Barcode>& hitParticlesMap,
    const Trajectories& trajectories, std::size_t tip,
    std::vector<ParticleHitCount>& particleHitCounts) {
  particleHitCounts.clear();

  if (!trajectories.hasTrajectory(tip)) {
    return;
  }

  trajectories.multiTrajectory().visitBackwards(tip, [&](const auto& state) {
    // no truth info with non-measurement state
    if (!state.typeFlags().isMeasurement()) {
      return true;
    }
    // skip outliers
    if (state.typeFlags().isOutlier()) {
      return true;
    }
    // register all particles that generated this hit
    IndexSourceLink sl =
        state.getUncalibratedSourceLink().template get<IndexSourceLink>();
    auto hitIndex = sl.index();
    for (const auto& [_, value] :
         makeRange(hitParticlesMap.equal_range(hitIndex))) {
      increaseHitCount(particleHitCounts, value);
    }
    return true;
  });
  sortHitCount(particleHitCounts);
}

void identifyContributingParticles(
    const IndexMultimap<ActsFatras::Barcode>& hitParticlesMap,
    const ConstTrackContainer::ConstTrackProxy& track,
    std::vector<ParticleHitCount>& particleHitCounts) {
  particleHitCounts.clear();

  for (const auto& state : track.trackStatesReversed()) {
    // no truth info with non-measurement state
    if (!state.typeFlags().isMeasurement()) {
      continue;
    }
    // skip outliers
    if (state.typeFlags().isOutlier()) {
      continue;
    }
    // register all particles that generated this hit
    IndexSourceLink sl =
        state.getUncalibratedSourceLink().template get<IndexSourceLink>();
    auto hitIndex = sl.index();
    for (const auto& [_, value] :
         makeRange(hitParticlesMap.equal_range(hitIndex))) {
      increaseHitCount(particleHitCounts, value);
    }
  }
  sortHitCount(particleHitCounts);
}

}  // namespace ActsExamples
