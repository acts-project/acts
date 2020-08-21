// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

// Local include(s).
#include "TestSpacePoint.hpp"

// Acts include(s).
#include "Acts/Seeding/IExperimentCuts.hpp"

/// Custom selection cuts for the test, used on the host
class TestHostCuts : public Acts::IExperimentCuts<TestSpacePoint> {
 public:
  /// Returns seed weight bonus/malus depending on detector considerations.
  /// @param bottom bottom space point of the current seed
  /// @param middle middle space point of the current seed
  /// @param top top space point of the current seed
  /// @return seed weight to be added to the seed's weight
  float seedWeight(const Acts::InternalSpacePoint<TestSpacePoint>& bottom,
                   const Acts::InternalSpacePoint<TestSpacePoint>& middle,
                   const Acts::InternalSpacePoint<TestSpacePoint>& top) const;

  /// @param weight the current seed weight
  /// @param bottom bottom space point of the current seed
  /// @param middle middle space point of the current seed
  /// @param top top space point of the current seed
  /// @return true if the seed should be kept, false if the seed should be
  /// discarded
  bool singleSeedCut(float weight,
                     const Acts::InternalSpacePoint<TestSpacePoint>& bottom,
                     const Acts::InternalSpacePoint<TestSpacePoint>&,
                     const Acts::InternalSpacePoint<TestSpacePoint>&) const;

  /// @param seeds contains pairs of weight and seed created for one middle
  /// space
  /// point
  /// @return vector of seeds that pass the cut
  std::vector<std::pair<
      float, std::unique_ptr<const Acts::InternalSeed<TestSpacePoint>>>>
  cutPerMiddleSP(
      std::vector<std::pair<
          float, std::unique_ptr<const Acts::InternalSeed<TestSpacePoint>>>>
          seeds) const;

};  // struct TestHostCuts
