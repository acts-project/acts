// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Seeding/IExperimentCuts.hpp"

#include <algorithm>

namespace Acts {
template <typename SpacePoint>
class GenericDetectorCuts : public IExperimentCuts<SpacePoint> {
 public:
  struct Config {
    int maxSeedSize = 5;
    float cutRadius = 200.;
    float cutWeight = 380.;
    float weight_outer = 400.;
    float weight_inner = 200.;
    float keepWeight = 200.;
    float minRadius = 43.;
  };

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
                     const InternalSpacePoint<SpacePoint>& middle,
                     const InternalSpacePoint<SpacePoint>& top) const;

  /// @param seeds contains pairs of weight and seed created for one middle
  /// space point
  /// @return vector of seeds that pass the cut
  std::vector<std::pair<float, std::unique_ptr<const InternalSeed<SpacePoint>>>>
  cutPerMiddleSP(
      std::vector<
          std::pair<float, std::unique_ptr<const InternalSeed<SpacePoint>>>>
          seeds) const;

 private:
  Config m_cfg;
};

template <typename SpacePoint>
float GenericDetectorCuts<SpacePoint>::seedWeight(
    const InternalSpacePoint<SpacePoint>& bottom,
    const InternalSpacePoint<SpacePoint>&,
    const InternalSpacePoint<SpacePoint>& top) const {
  float weight = 0;
  if (bottom.radius() > m_cfg.cutRadius) {
    weight = m_cfg.weight_outer;
  }
  if (top.radius() < m_cfg.cutRadius) {
    weight = m_cfg.weight_inner;
  }
  return weight;
}

template <typename SpacePoint>
bool GenericDetectorCuts<SpacePoint>::singleSeedCut(
    float weight, const InternalSpacePoint<SpacePoint>& b,
    const InternalSpacePoint<SpacePoint>&,
    const InternalSpacePoint<SpacePoint>&) const {
  return !(b.radius() > m_cfg.cutRadius && weight < m_cfg.cutWeight);
}

template <typename SpacePoint>
std::vector<std::pair<float, std::unique_ptr<const InternalSeed<SpacePoint>>>>
GenericDetectorCuts<SpacePoint>::cutPerMiddleSP(
    std::vector<
        std::pair<float, std::unique_ptr<const InternalSeed<SpacePoint>>>>
        seeds) const {
  std::vector<std::pair<float, std::unique_ptr<const InternalSeed<SpacePoint>>>>
      newSeedsVector;
  if (seeds.size() > 1) {
    newSeedsVector.push_back(std::move(seeds[0]));
    size_t itLength = std::min(seeds.size(), size_t(m_cfg.maxSeedSize));
    // don't cut first element
    for (size_t i = 1; i < itLength; i++) {
      if (seeds[i].first > m_cfg.keepWeight ||
          seeds[i].second->sp[0]->radius() > m_cfg.minRadius) {
        newSeedsVector.push_back(std::move(seeds[i]));
      }
    }
    return newSeedsVector;
  }
  return seeds;
}
}  // namespace Acts
