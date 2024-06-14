// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Surfaces/detail/MergeHelper.hpp"

namespace Acts::detail {

std::pair<ActsScalar, ActsScalar> mergedPhiSector(
    ActsScalar hlPhi1, ActsScalar avgPhi1, ActsScalar hlPhi2,
    ActsScalar avgPhi2, const Logger& logger, ActsScalar tolerance) {
  using namespace Acts::UnitLiterals;
  ActsScalar minPhi1 = detail::radian_sym(-hlPhi1 + avgPhi1);
  ActsScalar maxPhi1 = detail::radian_sym(hlPhi1 + avgPhi1);

  ACTS_VERBOSE("one: [" << minPhi1 / 1_degree << ", " << maxPhi1 / 1_degree
                        << "] ~> " << avgPhi1 / 1_degree << " +- "
                        << hlPhi1 / 1_degree);

  ActsScalar maxPhi2 = detail::radian_sym(hlPhi2 + avgPhi2);
  ActsScalar minPhi2 = detail::radian_sym(-hlPhi2 + avgPhi2);

  ACTS_VERBOSE("two: [" << minPhi2 / 1_degree << ", " << maxPhi2 / 1_degree
                        << "] ~> " << avgPhi2 / 1_degree << " +- "
                        << hlPhi2 / 1_degree);

  ACTS_VERBOSE("Checking for CCW or CW ordering");
  auto same = [tolerance](ActsScalar a, ActsScalar b) {
    return std::abs(a - b) < tolerance;
  };

  ActsScalar newMaxPhi{}, newMinPhi{};
  ActsScalar newHlPhi = hlPhi1 + hlPhi2;

  if (same(minPhi1, maxPhi2)) {
    ACTS_VERBOSE("-> CCW ordering: one is 'left' of two");

    newMinPhi = minPhi2;
    newMaxPhi = maxPhi1;

  } else if (same(maxPhi1, minPhi2)) {
    ACTS_VERBOSE("-> CW ordering: this is 'right' of other");
    newMinPhi = minPhi1;
    newMaxPhi = maxPhi2;

  } else {
    ACTS_ERROR("Phi ranges are incompatible");
    throw std::invalid_argument("Phi ranges are incompatible");
  }

  ActsScalar newAvgPhi = detail::radian_sym(newMinPhi + newHlPhi);

  ACTS_VERBOSE("merged: [" << newMinPhi / 1_degree << ", "
                           << newMaxPhi / 1_degree << "] ~> "
                           << newAvgPhi / 1_degree << " +- "
                           << newHlPhi / 1_degree);

  return {newHlPhi, newAvgPhi};
}

}  // namespace Acts::detail
