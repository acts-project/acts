// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "ActsFatras/Utilities/GammaDistribution.hpp"

#include <limits>
#include <random>

BOOST_AUTO_TEST_SUITE(FatrasGammaDistribution)

// Test if the standard-library gamma distribution gives the same results as
// underflow-safe GammaDistribution
BOOST_AUTO_TEST_CASE(GammaDistributionSequence) {
  // the same random number generators
  std::mt19937 rnd_std{30059};
  std::mt19937 rnd{30059};

  double alpha = 3.0;

  for (int i = 0; i < 20; ++i) {
    alpha /= 2.0;

    std::gamma_distribution<double> gDist_std(alpha, 1.);
    ActsFatras::GammaDistribution gDist(alpha, 1.);

    const auto u_std = gDist_std(rnd_std);
    const auto u = gDist(rnd);

    CHECK_CLOSE_OR_SMALL(u, u_std, 1e-6, std::numeric_limits<double>::min());
  }
}

BOOST_AUTO_TEST_SUITE_END()
