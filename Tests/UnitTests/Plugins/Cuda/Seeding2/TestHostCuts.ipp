// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

template <typename SpacePoint>
float TestHostCuts<SpacePoint>::seedWeight(const SpacePoint& bottom,
                                           const SpacePoint&,
                                           const SpacePoint& top) const {
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
bool TestHostCuts<SpacePoint>::singleSeedCut(float weight, const SpacePoint& b,
                                             const SpacePoint&,
                                             const SpacePoint&) const {
  return !(b.radius() > 150. && weight < 380.);
}

template <typename SpacePoint>
std::vector<typename Acts::CandidatesForMiddleSp<const SpacePoint>::value_type>
TestHostCuts<SpacePoint>::cutPerMiddleSP(
    std::vector<
        typename Acts::CandidatesForMiddleSp<const SpacePoint>::value_type>
        seedCandidates) const {
  std::vector<
      typename Acts::CandidatesForMiddleSp<const SpacePoint>::value_type>
      newSeedsVector;
  if (seedCandidates.size() <= 1) {
    return seedCandidates;
  }

  newSeedsVector.push_back(std::move(seedCandidates[0]));
  std::size_t itLength = std::min(seedCandidates.size(), std::size_t(5));
  // don't cut first element
  for (std::size_t i(1); i < itLength; i++) {
    float weight = seedCandidates[i].weight;
    const auto& bottom = seedCandidates[i].bottom;
    if (weight > 200. or bottom->radius() > 43.) {
      newSeedsVector.push_back(std::move(seedCandidates[i]));
    }
  }
  return newSeedsVector;
}
