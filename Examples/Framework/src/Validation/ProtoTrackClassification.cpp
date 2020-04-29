// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/Validation/ProtoTrackClassification.hpp"

#include <algorithm>

#include "ACTFW/Utilities/Range.hpp"

void FW::identifyContributingParticles(
    const IndexMultimap<ActsFatras::Barcode>& hitParticlesMap,
    const ProtoTrack& protoTrack,
    std::vector<FW::ParticleHitCount>& particleHitCount) {
  particleHitCount.clear();

  for (auto hitIndex : protoTrack) {
    // find all particles that generate this hit
    for (auto hitParticle : makeRange(hitParticlesMap.equal_range(hitIndex))) {
      auto particleId = hitParticle.second;
      // search for existing particle in the existing hit counts
      auto isSameParticle = [=](const ParticleHitCount& phc) {
        return (phc.particleId == particleId);
      };
      auto it = std::find_if(particleHitCount.begin(), particleHitCount.end(),
                             isSameParticle);
      // either increase count if we saw the particle before or add it
      if (it != particleHitCount.end()) {
        it->hitCount += 1;
      } else {
        particleHitCount.push_back({particleId, 1u});
      }
    }
  }

  // sort by hit count, i.e. majority particle first
  auto compareHitCount = [](const ParticleHitCount& lhs,
                            const ParticleHitCount& rhs) {
    return lhs.hitCount < rhs.hitCount;
  };
  std::sort(particleHitCount.begin(), particleHitCount.end(), compareHitCount);
}
