// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Seeding/IExperimentCuts.hpp"

#include <algorithm>

namespace Acts {
class ATLASCuts : public IExperimentCuts {
 public:
  /// Returns seed weight bonus/malus depending on detector considerations.
  /// @param bottom bottom space point of the current seed
  /// @param middle middle space point of the current seed
  /// @param top top space point of the current seed
  /// @return seed weight to be added to the seed's weight
  float seedWeight(const Acts::SpacePoint& bottom,
                   const Acts::SpacePoint& middle,
                   const Acts::SpacePoint& top) const;
  /// @param weight the current seed weight
  /// @param bottom bottom space point of the current seed
  /// @param middle middle space point of the current seed
  /// @param top top space point of the current seed
  /// @return true if the seed should be kept, false if the seed should be
  /// discarded
  bool singleSeedCut(float weight, const Acts::SpacePoint& bottom,
                     const Acts::SpacePoint&, const Acts::SpacePoint&) const;

  /// @param seeds contains pairs of weight and seed created for one middle
  /// space
  /// point
  /// @return vector of seeds that pass the cut
  std::vector<std::pair<float, std::unique_ptr<const InternalSeed>>>
  cutPerMiddleSP(
      std::vector<std::pair<float, std::unique_ptr<const InternalSeed>>> seeds)
      const;
};

float ATLASCuts::seedWeight(const Acts::SpacePoint& bottom,
                            const Acts::SpacePoint&,
                            const Acts::SpacePoint& top) const {
  float weight = 0;
  if (bottom.radius() > 150) {
    weight = 400;
  }
  if (top.radius() < 150) {
    weight = 200;
  }
  return weight;
}

bool ATLASCuts::singleSeedCut(float weight, const Acts::SpacePoint& b,
                              const Acts::SpacePoint&,
                              const Acts::SpacePoint&) const {
  return !(b.radius() > 150. && weight < 380.);
}

std::vector<std::pair<float, std::unique_ptr<const InternalSeed>>>
ATLASCuts::cutPerMiddleSP(
    std::vector<std::pair<float, std::unique_ptr<const InternalSeed>>> seeds)
    const {
  std::vector<std::pair<float, std::unique_ptr<const InternalSeed>>>
      newSeedsVector;
  if (seeds.size() > 1) {
    newSeedsVector.push_back(std::move(seeds[0]));
    size_t itLength = std::min(seeds.size(), size_t(5));
    // don't cut first element
    for (size_t i = 1; i < itLength; i++) {
      if (seeds[i].first > 200. || seeds[i].second->sp[0]->radius() > 43.) {
        newSeedsVector.push_back(std::move(seeds[i]));
      }
    }
    return newSeedsVector;
  }
  return seeds;
}
}  // namespace Acts
