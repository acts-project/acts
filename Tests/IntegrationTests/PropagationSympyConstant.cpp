// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/RiddersPropagator.hpp"
#include "Acts/Propagator/SympyStepper.hpp"

#include "PropagationDatasets.hpp"
#include "PropagationTests.hpp"

namespace {

namespace ds = ActsTests::PropagationDatasets;

using namespace Acts;
using namespace UnitLiterals;

using MagneticField = ConstantBField;
using Stepper = SympyStepper;
using TestPropagator = Propagator<Stepper>;
using RiddersPropagator = RiddersPropagator<TestPropagator>;

// absolute parameter tolerances for position, direction, and absolute momentum
constexpr auto epsPos = 1_um;
constexpr auto epsDir = 0.125_mrad;
constexpr auto epsMom = 1_eV;

const auto geoCtx = GeometryContext::dangerouslyDefaultConstruct();
const MagneticFieldContext magCtx;

inline TestPropagator makePropagator(double bz) {
  auto magField = std::make_shared<MagneticField>(Vector3(0.0, 0.0, bz));
  Stepper stepper(std::move(magField));
  return TestPropagator(std::move(stepper));
}

}  // namespace

BOOST_AUTO_TEST_SUITE(PropagationSympyConstant)

// check that the propagation is reversible and self-consistent

BOOST_DATA_TEST_CASE(ForwardBackward,
                     ds::phi* ds::thetaWithoutBeam* ds::absMomentum*
                         ds::chargeNonZero* ds::pathLength* ds::magneticField,
                     phi, theta, p, q, s, bz) {
  runForwardBackwardTest(makePropagator(bz), geoCtx, magCtx,
                         makeParametersCurvilinear(phi, theta, p, q), s, epsPos,
                         epsDir, epsMom);
}

// TODO implement jacobian/covariance tests. right now covariance comparisons
// with RiddersPropagator fail because of numerical imprecisions

BOOST_AUTO_TEST_SUITE_END()
