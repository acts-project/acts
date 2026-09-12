// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsFatras/Synthetic/SeedingTruth.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace ActsFatras {

Synthetic::EventSummary Synthetic::summarize(const Event& event,
                                             const float ptThreshold) {
  EventSummary summary;
  summary.spacePoints = event.spacePoints.size();
  for (const GeneratedParticle& particle : event.particles) {
    if (particle.primary()) {
      ++summary.primaries;
      summary.primaryHits += particle.numHits;
      if (particle.pt >= ptThreshold && particle.numHits >= 3) {
        ++summary.seedablePrimaries;
      }
    } else {
      ++summary.secondaries;
      summary.secondaryHits += particle.numHits;
    }
  }
  return summary;
}

Synthetic::SeedingSummary Synthetic::evaluateSeeds(
    const Event& event, const Acts::SeedContainer& seeds,
    const float ptThreshold, const std::size_t minTrueSpacePoints) {
  const auto particleColumn =
      event.spacePoints.column<std::uint32_t>("particleId");

  SeedingSummary summary;
  summary.seeds = seeds.size();
  std::vector<bool> matched(event.particles.size(), false);
  for (const auto seed : seeds) {
    const auto indices = seed.spacePointIndices();
    if (indices.empty() || indices.size() < minTrueSpacePoints) {
      continue;
    }
    // The particle contributing the most space points. Scanning the seed
    // against itself beats a map, a seed holding a handful of them.
    std::uint32_t particle = 0;
    std::size_t best = 0;
    for (const Acts::SpacePointIndex index : indices) {
      const std::uint32_t candidate =
          event.spacePoints[index].extra(particleColumn);
      const auto count = static_cast<std::size_t>(
          std::ranges::count_if(indices, [&](Acts::SpacePointIndex other) {
            return event.spacePoints[other].extra(particleColumn) == candidate;
          }));
      if (count > best) {
        best = count;
        particle = candidate;
      }
    }
    if (best < minTrueSpacePoints) {
      continue;
    }
    if (const GeneratedParticle& info = event.particles[particle];
        !info.primary() || info.pt < ptThreshold) {
      continue;
    }
    ++summary.trueSeeds;
    matched[particle] = true;
  }
  summary.matchedPrimaries =
      static_cast<std::size_t>(std::ranges::count(matched, true));
  return summary;
}

}  // namespace ActsFatras
