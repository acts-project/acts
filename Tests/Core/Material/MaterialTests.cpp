// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// clang-format off
#define BOOST_TEST_MODULE Material Tests
#include <boost/test/included/unit_test.hpp>
// clang-format on

#include <climits>

#include "Acts/Material/Material.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Units.hpp"

namespace Acts {

namespace Test {

  // the maximum tolerance is half the accuracy
  float elMaxTolerance = 0.5 / float(UCHAR_MAX);

  namespace au = Acts::units;

  // first test correct boolean behavior
  BOOST_AUTO_TEST_CASE(Material_boolean_test)
  {
    Material vacuum;
    BOOST_CHECK_EQUAL(bool(vacuum), false);

    Material something(1., 2., 3., 4., 5);
    BOOST_CHECK_EQUAL(bool(something), true);
  }

  // now test thge construction and units
  BOOST_AUTO_TEST_CASE(Material_construction_and_units)
  {

    // density at room temperature
    float X0  = 9.370 * au::_cm;
    float L0  = 46.52 * au::_cm;
    float Z   = 14.;
    float A   = 28.0855;
    float rho = 2.329 * au::_g / (au::_cm * au::_cm * au::_cm);

    Material silicon(X0, L0, A, Z, rho);
    CHECK_CLOSE_REL(silicon.X0(), 93.70 * au::_mm, 0.001);
    CHECK_CLOSE_REL(silicon.L0(), 465.2 * au::_mm, 0.001);
    CHECK_CLOSE_REL(silicon.Z(), 14., 0.001);
    CHECK_CLOSE_REL(silicon.A(), 28.0855, 0.001);
    CHECK_CLOSE_REL(silicon.rho(),
                    0.002329 * au::_g / (au::_mm * au::_mm * au::_mm),
                    0.001);
    CHECK_CLOSE_REL(silicon.zOverAtimesRho(), 14. / 28.0855 * 0.002329, 0.0001);

    Material copiedSilicon(silicon);
    BOOST_CHECK_EQUAL(silicon, copiedSilicon);

    Material moveCopiedSilicon(std::move(copiedSilicon));
    BOOST_CHECK_EQUAL(silicon, moveCopiedSilicon);

    Material assignedSilicon = silicon;
    BOOST_CHECK_EQUAL(silicon, assignedSilicon);

    Material moveAssignedSilicon = std::move(assignedSilicon);
    BOOST_CHECK_EQUAL(silicon, moveAssignedSilicon);
  }
}
}
