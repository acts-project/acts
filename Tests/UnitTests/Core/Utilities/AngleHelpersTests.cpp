// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/AngleHelpers.hpp"

#include <cmath>
#include <numbers>

namespace bd = boost::unit_test::data;

using Acts::AngleHelpers::etaFromTheta;
using Acts::AngleHelpers::thetaFromEta;

BOOST_AUTO_TEST_SUITE(AngleHelpers)

BOOST_AUTO_TEST_CASE(EtaThetaConversion) {
  CHECK_CLOSE_ABS(0.0, etaFromTheta(std::numbers::pi / 2), 1e-6);

  CHECK_CLOSE_ABS(1.0, etaFromTheta(thetaFromEta(1.0)), 1e-6);

  CHECK_CLOSE_ABS(1.0, thetaFromEta(etaFromTheta(1.0)), 1e-6);
}

BOOST_DATA_TEST_CASE(EtaFromThetaRobustness, bd::xrange(1, 100, 1), exponent) {
  {
    float thetaRight = std::pow(10.0f, -exponent);
    float etaRight = etaFromTheta(thetaRight);
    BOOST_CHECK(!std::isnan(etaRight));

    float thetaLeft = std::numbers::pi_v<float> - thetaRight;
    float etaLeft = etaFromTheta(thetaLeft);
    BOOST_CHECK(!std::isnan(etaLeft));
  }

  {
    double thetaRight = std::pow(10.0, -exponent);
    double etaRight = etaFromTheta(thetaRight);
    BOOST_CHECK(!std::isnan(etaRight));

    double thetaLeft = std::numbers::pi - thetaRight;
    double etaLeft = etaFromTheta(thetaLeft);
    BOOST_CHECK(!std::isnan(etaLeft));
  }
}

BOOST_DATA_TEST_CASE(ThetaFromEtaRobustness, bd::xrange(1.0, 100.0, 1.0),
                     etaRight) {
  {
    float thetaRight = thetaFromEta(etaRight);
    BOOST_CHECK(!std::isnan(thetaRight));

    float etaLeft = -etaRight;
    float thetaLeft = thetaFromEta(etaLeft);
    BOOST_CHECK(!std::isnan(thetaLeft));
  }

  {
    double thetaRight = thetaFromEta(etaRight);
    BOOST_CHECK(!std::isnan(thetaRight));

    double etaLeft = -etaRight;
    double thetaLeft = thetaFromEta(etaLeft);
    BOOST_CHECK(!std::isnan(thetaLeft));
  }
}

BOOST_AUTO_TEST_SUITE_END()
