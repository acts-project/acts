// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/**
 * @file MeasurementTests.cpp
 */

// Boost include(s)
#define BOOST_TEST_MODULE Measurement Tests
#include <boost/test/included/unit_test.hpp>

// ACTS include(s)
#include "ACTS/EventData/Measurement.hpp"
#include "ACTS/Surfaces/CylinderSurface.hpp"
#include "ACTS/Utilities/ParameterDefinitions.hpp"

namespace Acts {
namespace Test {
  template <ParID_t... params>
  using Measurement_t = Measurement<unsigned long int, params...>;

  /**
   * @brief Unit test for creation of Measurement object
   */
  BOOST_AUTO_TEST_CASE(measurement_initialization)
  {
    CylinderSurface   cylinder(nullptr, 3, 10);
    ActsSymMatrixD<2> cov;
    cov << 0.04, 0, 0, 0.1;
    Measurement_t<ParDef::eLOC_0, ParDef::eLOC_1> m(
        cylinder, 0, std::move(cov), -0.1, 0.45);
  }
}  // end of namespace Test
}  // end of namespace Acts
