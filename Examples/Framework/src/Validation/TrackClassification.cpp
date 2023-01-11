// This file is part of the Acts project.
//
// Copyright (C) 2019-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Validation/TrackClassification.hpp"

#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/Utilities/Range.hpp"

#include <algorithm>

namespace {

/// Increase the hit count for the given particle id by one.
inline void increaseHitCount(
    std::vector<ActsExamples::ParticleHitCount>& particleHitCounts,
    ActsFatras::Barcode particleId) {
  // linear search since there is no ordering
  auto it = std::find_if(particleHitCounts.begin(), particleHitCounts.end(),
                         [=](const ActsExamples::ParticleHitCount& phc) {
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
inline void sortHitCount(
    std::vector<ActsExamples::ParticleHitCount>& particleHitCounts) {
  std::sort(particleHitCounts.begin(), particleHitCounts.end(),
            [](const ActsExamples::ParticleHitCount& lhs,
               const ActsExamples::ParticleHitCount& rhs) {
              return (lhs.hitCount > rhs.hitCount);
            });
}

}  // namespace

void ActsExamples::identifyContributingParticles(
    const IndexMultimap<ActsFatras::Barcode>& hitParticlesMap,
    const ProtoTrack& protoTrack,
    std::vector<ActsExamples::ParticleHitCount>& particleHitCounts) {
  particleHitCounts.clear();

  for (auto hitIndex : protoTrack) {
    // register all particles that generated this hit
    for (auto hitParticle : makeRange(hitParticlesMap.equal_range(hitIndex))) {
      increaseHitCount(particleHitCounts, hitParticle.second);
    }
  }
  sortHitCount(particleHitCounts);
}

void ActsExamples::identifyContributingParticles(
    const IndexMultimap<ActsFatras::Barcode>& hitParticlesMap,
    const Trajectories& trajectories, size_t tip,
    std::vector<ParticleHitCount>& particleHitCounts) {
  particleHitCounts.clear();

  if (not trajectories.hasTrajectory(tip)) {
    return;
  }

  trajectories.multiTrajectory().visitBackwards(tip, [&](const auto& state) {
    // no truth info with non-measurement state
    if (not state.typeFlags().test(Acts::TrackStateFlag::MeasurementFlag)) {
      return true;
    }
    // register all particles that generated this hit
    const IndexSourceLink& sl =
        state.uncalibratedSourceLink().template get<IndexSourceLink>();
    auto hitIndex = sl.index();
    for (auto hitParticle : makeRange(hitParticlesMap.equal_range(hitIndex))) {
      increaseHitCount(particleHitCounts, hitParticle.second);
    }
    return true;
  });
  sortHitCount(particleHitCounts);
}
