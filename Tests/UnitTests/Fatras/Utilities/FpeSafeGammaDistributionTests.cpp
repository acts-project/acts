// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "ActsFatras/Utilities/detail/FpeSafeGammaDistribution.hpp"

#include <limits>
#include <random>

BOOST_AUTO_TEST_SUITE(FatrasFpeSafeGammaDistribution)

// Test if the standard-library gamma distribution gives the same results as
// underflow-safe  FpeSafeGammaDistribution
BOOST_AUTO_TEST_CASE(FpeSafeGammaDistributionSequence) {
  std::mt19937 rnd{30059};

  const int num = 20;
#ifdef __GLIBCXX__
  // results from std::gamma_distribution<double>, in libstdc++
  double results[num] = {1.4631785,
                         4.3811862e-01,
                         2.1447861e-01,
                         1.1303697e-01,
                         3.5990972e-06,
                         6.6019759e-09,
                         2.2688236e-15,
                         2.3389734e-11,
                         3.1843972e-02,
                         2.8021629e-85,
                         3.1397155e-185,
                         3.1850069e-284,
                         0.,
                         2.2658164e-53,
                         0.,
                         0.,
                         0.,
                         0.,
                         0.,
                         0.};
#elif defined _LIBCPP_VERSION
  // results from std::gamma_distribution<double>, in libc++
  double results[num] = {3.3655543,
                         6.7167479e-01,
                         1.3436782,
                         8.9406670e-04,
                         7.4754048e-05,
                         9.8973516e-20,
                         2.9430952e-08,
                         8.4810760e-07,
                         7.2823453e-41,
                         0.,
                         0.,
                         3.0934012e-14,
                         0.,
                         0.,
                         0.,
                         0.,
                         0.,
                         0.,
                         0.,
                         0.};
#else
  // unknown library, will fail
  double results[num];
  for (int i = 0; i < num; ++i) {
    results[i] = -1.0;
  }
#endif

  double alpha = 3.0;
  for (int i = 0; i < num; ++i) {
    alpha /= 2.0;
    ActsFatras::detail::FpeSafeGammaDistribution gDist(alpha, 1.);
    const auto u = gDist(rnd);
    CHECK_CLOSE_OR_SMALL(u, results[i], 1e-6,
                         std::numeric_limits<double>::min());
  }
}

BOOST_AUTO_TEST_SUITE_END()
