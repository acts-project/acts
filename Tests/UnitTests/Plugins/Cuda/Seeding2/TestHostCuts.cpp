// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// Local include(s).
#include "TestHostCuts.hpp"

float TestHostCuts::seedWeight(
    const Acts::InternalSpacePoint<TestSpacePoint>& bottom,
    const Acts::InternalSpacePoint<TestSpacePoint>&,
    const Acts::InternalSpacePoint<TestSpacePoint>& top) const {
  float weight = 0;
  if (bottom.radius() > 150) {
    weight = 400;
  }
  if (top.radius() < 150) {
    weight = 200;
  }
  return weight;
}

bool TestHostCuts::singleSeedCut(
    float weight, const Acts::InternalSpacePoint<TestSpacePoint>& b,
    const Acts::InternalSpacePoint<TestSpacePoint>&,
    const Acts::InternalSpacePoint<TestSpacePoint>&) const {
  return !(b.radius() > 150. && weight < 380.);
}

std::vector<
    std::pair<float, std::unique_ptr<const Acts::InternalSeed<TestSpacePoint>>>>
TestHostCuts::cutPerMiddleSP(
    std::vector<std::pair<
        float, std::unique_ptr<const Acts::InternalSeed<TestSpacePoint>>>>
        seeds) const {
  std::vector<std::pair<
      float, std::unique_ptr<const Acts::InternalSeed<TestSpacePoint>>>>
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
