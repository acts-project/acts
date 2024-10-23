// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/AngleHelpers.hpp"

#include <numbers>

using Acts::AngleHelpers::etaFromTheta;
using Acts::AngleHelpers::thetaFromEta;

BOOST_AUTO_TEST_SUITE(AngleHelpers)

BOOST_AUTO_TEST_CASE(EtaThetaConversion) {
  CHECK_CLOSE_ABS(0.0, etaFromTheta(std::numbers::pi / 2), 1e-6);

  CHECK_CLOSE_ABS(1.0, etaFromTheta(thetaFromEta(1.0)), 1e-6);

  CHECK_CLOSE_ABS(1.0, thetaFromEta(etaFromTheta(1.0)), 1e-6);
}

BOOST_AUTO_TEST_SUITE_END()
