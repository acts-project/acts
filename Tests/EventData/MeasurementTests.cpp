// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// Boost include(s)
#define BOOST_TEST_MODULE Measurement Tests
#include <boost/test/included/unit_test.hpp>

// Acts include(s)
#include "Acts/EventData/Measurement.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Utilities/ParameterDefinitions.hpp"

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
}  // namespace Test
}  // namespace Acts
