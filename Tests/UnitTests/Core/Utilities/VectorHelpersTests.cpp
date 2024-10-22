// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/VectorHelpers.hpp"

#include <numbers>

using Acts::VectorHelpers::eta;
using Acts::VectorHelpers::theta;

BOOST_AUTO_TEST_SUITE(AngleHelpers)

BOOST_AUTO_TEST_CASE(EtaFromVector) {
  CHECK_CLOSE_ABS(0.0, eta(Acts::Vector3{1, 0, 0}), 1e-6);
}

BOOST_AUTO_TEST_CASE(ThetaFromVector) {
  CHECK_CLOSE_ABS(std::numbers::pi / 2, theta(Acts::Vector3{1, 0, 0}), 1e-6);
}

BOOST_AUTO_TEST_SUITE_END()
