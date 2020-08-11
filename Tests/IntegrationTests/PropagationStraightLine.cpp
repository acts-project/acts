// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"

#include <limits>

#include "PropagationDatasets.hpp"
#include "PropagationTests.hpp"

namespace {

namespace ds = ActsTests::PropagationDatasets;
using namespace Acts::UnitLiterals;

constexpr auto epsPos = 1_um;
constexpr auto epsDir = 0.125_mrad;
constexpr auto epsMom = 1_eV;
constexpr bool showDebug = false;

const Acts::GeometryContext geoCtx;
const Acts::MagneticFieldContext magCtx;
const Acts::StraightLineStepper stepper;
const Acts::Propagator<Acts::StraightLineStepper> propagator(stepper);

}  // namespace

BOOST_AUTO_TEST_SUITE(PropagationStraightLine)

BOOST_DATA_TEST_CASE(
    ForwardBackward,
    ds::phi* ds::theta* ds::absMomentum* ds::chargeNonZero* ds::pathLength, phi,
    theta, p, q, s) {
  // phi is ill-defined in forward/backward tracks
  const auto phiFixed = ((0 < theta) and (theta < M_PI)) ? phi : 0.0;
  const auto initial = makeParametersCurvilinear(phiFixed, theta, p, q);
  runForwardBackwardTest(propagator, geoCtx, magCtx, initial, s, epsPos, epsDir,
                         epsMom, showDebug);
}

BOOST_AUTO_TEST_SUITE_END()
