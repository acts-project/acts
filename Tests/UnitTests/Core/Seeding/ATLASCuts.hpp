// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Seeding/CandidatesForMiddleSp.hpp"
#include "Acts/Seeding/IExperimentCuts.hpp"

namespace Acts {

class ATLASCuts : public IExperimentCuts {
 public:
  /// Returns seed weight bonus/malus depending on detector considerations.
  /// @param bottom bottom space point of the current seed
  /// @param middle middle space point of the current seed
  /// @param top top space point of the current seed
  /// @return seed weight to be added to the seed's weight
  float seedWeight(const ConstInternalSpacePointProxy& bottom,
                   const ConstInternalSpacePointProxy& middle,
                   const ConstInternalSpacePointProxy& top) const override;

  /// @param weight the current seed weight
  /// @param bottom bottom space point of the current seed
  /// @param middle middle space point of the current seed
  /// @param top top space point of the current seed
  /// @return true if the seed should be kept, false if the seed should be
  /// discarded
  bool singleSeedCut(
      float weight, const ConstInternalSpacePointProxy& bottom,
      const ConstInternalSpacePointProxy& /*middle*/,
      const ConstInternalSpacePointProxy& /*top*/) const override;

  /// @param seedCandidates contains collection of seed candidates created for one middle
  /// space point in a std::tuple format
  /// @return vector of seeds that pass the cut
  std::vector<TripletCandidate> cutPerMiddleSP(
      std::vector<TripletCandidate> seedCandidates,
      const InternalSpacePointContainer& spacePoints) const override;
};

float ATLASCuts::seedWeight(const ConstInternalSpacePointProxy& bottom,
                            const ConstInternalSpacePointProxy& /*middle*/,
                            const ConstInternalSpacePointProxy& top) const {
  float weight = 0;
  if (bottom.radius() > 150) {
    weight = 400;
  }
  if (top.radius() < 150) {
    weight = 200;
  }
  return weight;
}

bool ATLASCuts::singleSeedCut(float weight,
                              const ConstInternalSpacePointProxy& b,
                              const ConstInternalSpacePointProxy& /*m*/,
                              const ConstInternalSpacePointProxy& /*t*/) const {
  return !(b.radius() > 150. && weight < 380.);
}

std::vector<TripletCandidate> ATLASCuts::cutPerMiddleSP(
    std::vector<TripletCandidate> seedCandidates,
    const InternalSpacePointContainer& spacePoints) const {
  std::vector<TripletCandidate> newSeedsVector;
  if (seedCandidates.size() <= 1) {
    return seedCandidates;
  }

  newSeedsVector.push_back(seedCandidates[0]);
  std::size_t itLength = std::min(seedCandidates.size(), std::size_t{5});
  // don't cut first element
  for (std::size_t i(1); i < itLength; i++) {
    float weight = seedCandidates[i].weight;
    const auto& bottom = spacePoints.at(seedCandidates[i].bottom);
    if (weight > 200. || bottom.radius() > 43.) {
      newSeedsVector.push_back(seedCandidates[i]);
    }
  }
  return newSeedsVector;
}

}  // namespace Acts
