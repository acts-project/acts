// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/RiddersPropagator.hpp"
#include "Acts/Propagator/SympyStepper.hpp"
#include "ActsTests/CommonHelpers/PredefinedMaterials.hpp"

#include <utility>

#include "PropagationDatasets.hpp"
#include "PropagationTests.hpp"

namespace {

namespace ds = ActsTests::PropagationDatasets;

using namespace Acts;
using namespace UnitLiterals;

using MagneticField = ConstantBField;
using Stepper = SympyStepper;
using TestPropagator = Propagator<Stepper, Navigator>;
using RiddersPropagator = RiddersPropagator<TestPropagator>;

// absolute parameter tolerances for position, direction, and absolute momentum
constexpr auto epsPos = 10_um;
constexpr auto epsDir = 1_mrad;
constexpr auto epsMom = 5_MeV;

const auto geoCtx = GeometryContext::dangerouslyDefaultConstruct();
const MagneticFieldContext magCtx;

inline TestPropagator makePropagator(
    double bz, std::shared_ptr<const TrackingGeometry> geo) {
  auto magField = std::make_shared<MagneticField>(Vector3(0.0, 0.0, bz));
  Stepper stepper(std::move(magField));
  return TestPropagator(std::move(stepper), Navigator({std::move(geo)}));
}

}  // namespace

BOOST_AUTO_TEST_SUITE(PropagationSympyDenseConstant)

// check that the propagation is reversible and self-consistent

BOOST_DATA_TEST_CASE(ForwardBackward,
                     ds::phi* ds::thetaWithoutBeam* ds::absMomentum*
                         ds::chargeNonZero* ds::pathLength* ds::magneticField,
                     phi, theta, p, q, s, bz) {
  runForwardBackwardTest(makePropagator(bz, createDenseBlock(geoCtx)), geoCtx,
                         magCtx, makeParametersCurvilinear(phi, theta, p, q), s,
                         epsPos, epsDir, epsMom);
}

// check effects on the parameter covariance

BOOST_DATA_TEST_CASE(DenseTelescopeCovariance,
                     ds::absMomentum* ds::chargeNonZero, p, q) {
  const double bz = 0_T;
  const Material material = ActsTests::makeLiquidArgon();
  const double thickness = 1_m;

  auto [geo, surfaces] = createDenseTelescope(geoCtx, material, thickness);

  auto propagator = makePropagator(bz, std::move(geo));

  runDenseForwardTest(propagator, geoCtx, magCtx, p, q, *surfaces.back(),
                      material, thickness);
}

BOOST_AUTO_TEST_SUITE_END()
