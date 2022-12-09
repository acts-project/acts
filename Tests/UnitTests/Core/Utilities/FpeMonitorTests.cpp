// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Utilities/FpeMonitor.hpp"

#include <optional>

namespace utf = boost::unit_test;

namespace Acts {

namespace Test {

// BOOST_AUTO_TEST_CASE(FpeMonitorTest) {
// FpeMonitor mon;
// volatile const double x = -1;
// printf("y = %f\n", sqrt(x));
// }

BOOST_AUTO_TEST_CASE(FpeMonitorTestDivByZero) {
  // {
  // FpeMonitor mon;
  // volatile double z = 0;
  // volatile double x = 1 / z;
  // std::cout << "x: " << x << std::endl;
  // }

  {
    FpeMonitor mon;
    volatile float v = 0;
    volatile double w = std::numeric_limits<double>::max();
    v = 2 * w;
    std::cout << v << std::endl;
  }
}  // namespace Test

}  // namespace Test

}  // namespace Acts
