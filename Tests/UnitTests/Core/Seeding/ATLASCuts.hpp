// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <algorithm>

#include "Acts/Seeding/IExperimentCuts.hpp"

namespace Acts {
template <typename SpacePoint>
class ATLASCuts : public IExperimentCuts<SpacePoint> {
 public:
  /// Returns seed weight bonus/malus depending on detector considerations.
  /// @param bottom bottom space point of the current seed
  /// @param middle middle space point of the current seed
  /// @param top top space point of the current seed
  /// @return seed weight to be added to the seed's weight
  float seedWeight(const InternalSpacePoint<SpacePoint>& bottom,
                   const InternalSpacePoint<SpacePoint>& middle,
                   const InternalSpacePoint<SpacePoint>& top) const;
  /// @param weight the current seed weight
  /// @param bottom bottom space point of the current seed
  /// @param middle middle space point of the current seed
  /// @param top top space point of the current seed
  /// @return true if the seed should be kept, false if the seed should be
  /// discarded
  bool singleSeedCut(float weight, const InternalSpacePoint<SpacePoint>& bottom,
                     const InternalSpacePoint<SpacePoint>&,
                     const InternalSpacePoint<SpacePoint>&) const;

  /// @param seeds contains pairs of weight and seed created for one middle
  /// space
  /// point
  /// @return vector of seeds that pass the cut
  std::vector<std::pair<float, std::unique_ptr<const InternalSeed<SpacePoint>>>>
  cutPerMiddleSP(
      std::vector<
          std::pair<float, std::unique_ptr<const InternalSeed<SpacePoint>>>>
          seeds) const;
};

template <typename SpacePoint>
float ATLASCuts<SpacePoint>::seedWeight(
    const InternalSpacePoint<SpacePoint>& bottom,
    const InternalSpacePoint<SpacePoint>&,
    const InternalSpacePoint<SpacePoint>& top) const {
  float weight = 0;
  if (bottom.radius() > 150) {
    weight = 400;
  }
  if (top.radius() < 150) {
    weight = 200;
  }
  return weight;
}

template <typename SpacePoint>
bool ATLASCuts<SpacePoint>::singleSeedCut(
    float weight, const InternalSpacePoint<SpacePoint>& b,
    const InternalSpacePoint<SpacePoint>&,
    const InternalSpacePoint<SpacePoint>&) const {
  return !(b.radius() > 150. && weight < 380.);
}

template <typename SpacePoint>
std::vector<std::pair<float, std::unique_ptr<const InternalSeed<SpacePoint>>>>
ATLASCuts<SpacePoint>::cutPerMiddleSP(
    std::vector<
        std::pair<float, std::unique_ptr<const InternalSeed<SpacePoint>>>>
        seeds) const {
  std::vector<std::pair<float, std::unique_ptr<const InternalSeed<SpacePoint>>>>
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
