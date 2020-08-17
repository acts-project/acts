// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/EventData/SimMultiTrajectory.hpp"

std::vector<ActsExamples::ParticleHitCount>
ActsExamples::SimMultiTrajectory::identifyMajorityParticle(
    const size_t& entryIndex) const {
  std::vector<ActsExamples::ParticleHitCount> particleHitCount;
  particleHitCount.reserve(10);
  if (not m_trackTips.empty()) {
    if (not hasTrajectory(entryIndex)) {
      return particleHitCount;
    }
    m_multiTrajectory.visitBackwards(entryIndex, [&](const auto& state) {
      // No truth info with non-measurement state
      if (not state.typeFlags().test(Acts::TrackStateFlag::MeasurementFlag)) {
        return true;
      }
      // Find the truth particle associated with this state
      const auto particleId = state.uncalibrated().truthHit().particleId();
      // Find if the particle already exists
      auto it = std::find_if(particleHitCount.begin(), particleHitCount.end(),
                             [=](const ActsExamples::ParticleHitCount& phc) {
                               return phc.particleId == particleId;
                             });

      // Either increase count if we saw the particle before or add it
      if (it != particleHitCount.end()) {
        it->hitCount += 1;
      } else {
        particleHitCount.push_back({particleId, 1u});
      }
      return true;
    });
  }
  if (not particleHitCount.empty()) {
    // sort by hit count, i.e. majority particle first
    std::sort(particleHitCount.begin(), particleHitCount.end(),
              [](const ActsExamples::ParticleHitCount& lhs,
                 const ActsExamples::ParticleHitCount& rhs) {
                return lhs.hitCount > rhs.hitCount;
              });
  }

  return particleHitCount;
}
