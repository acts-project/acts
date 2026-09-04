// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Units.hpp"
#include "Acts/Seeding/detail/SpacePointGridPhiBinning.hpp"

#include <stdexcept>

using namespace Acts::UnitLiterals;

namespace Acts::Test {

namespace {

// A physically sensible baseline: the minimum-pT helix radius is well above
// rMax / 2, so the phi bin count is driven by the helix deflection.
detail::SpacePointGridPhiBinningConfig makeConfig() {
  detail::SpacePointGridPhiBinningConfig config;
  config.minPt = 500_MeV;
  config.bFieldInZ = 2_T;
  config.rMax = 200_mm;
  config.deltaRMax = 100_mm;
  config.impactMax = 10_mm;
  config.phiBinDeflectionCoverage = 1;
  config.maxPhiBins = 10000;
  return config;
}

}  // namespace

BOOST_AUTO_TEST_SUITE(SpacePointGridPhiBinningTests)

BOOST_AUTO_TEST_CASE(ZeroFieldUsesMaxPhiBins) {
  auto config = makeConfig();
  config.bFieldInZ = 0;
  config.maxPhiBins = 42;
  BOOST_CHECK_EQUAL(detail::computeSpacePointGridPhiBins(config), 42);
}

BOOST_AUTO_TEST_CASE(BaselineGivesFinitePositiveBinCount) {
  const auto config = makeConfig();
  const int phiBins = detail::computeSpacePointGridPhiBins(config);
  BOOST_CHECK_GT(phiBins, 1);
  BOOST_CHECK_LE(phiBins, config.maxPhiBins);
}

BOOST_AUTO_TEST_CASE(MaxPhiBinsCapsResult) {
  auto config = makeConfig();
  const int uncapped = detail::computeSpacePointGridPhiBins(config);
  config.maxPhiBins = uncapped - 1;
  BOOST_CHECK_EQUAL(detail::computeSpacePointGridPhiBins(config),
                    config.maxPhiBins);
}

BOOST_AUTO_TEST_CASE(DeflectionCoverageIncreasesBinCount) {
  auto config = makeConfig();
  const int singleCoverage = detail::computeSpacePointGridPhiBins(config);
  config.phiBinDeflectionCoverage = 2;
  const int doubleCoverage = detail::computeSpacePointGridPhiBins(config);
  BOOST_CHECK_GT(doubleCoverage, singleCoverage);
}

// Regression for the robustness fix of PR 5696: an impact parameter at or
// beyond rMax - deltaRMax must fall back to a coarse binning instead of
// producing a NaN and an invalid (zero) bin count.
BOOST_AUTO_TEST_CASE(LargeImpactParameterStaysFinite) {
  auto config = makeConfig();
  config.impactMax = config.rMax;
  const int phiBins = detail::computeSpacePointGridPhiBins(config);
  BOOST_CHECK_GT(phiBins, 0);
  BOOST_CHECK_LE(phiBins, config.maxPhiBins);
}

BOOST_AUTO_TEST_CASE(TooSmallHelixRadiusThrows) {
  auto config = makeConfig();
  // minHelixRadius = minPt / bFieldInZ falls below rMax / 2
  config.minPt = 10_MeV;
  BOOST_CHECK_THROW(detail::computeSpacePointGridPhiBins(config),
                    std::domain_error);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Acts::Test
