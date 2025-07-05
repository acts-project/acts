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
#include "Acts/Tests/CommonHelpers/PredefinedMaterials.hpp"

#include <utility>

#include "PropagationDatasets.hpp"
#include "PropagationTests.hpp"

namespace {

namespace ds = ActsTests::PropagationDatasets;
using namespace Acts::UnitLiterals;

using MagneticField = Acts::ConstantBField;
using Stepper = Acts::SympyStepper;
using Propagator = Acts::Propagator<Stepper, Acts::Navigator>;
using RiddersPropagator = Acts::RiddersPropagator<Propagator>;

// absolute parameter tolerances for position, direction, and absolute momentum
constexpr auto epsPos = 10_um;
constexpr auto epsDir = 1_mrad;
constexpr auto epsMom = 5_MeV;

const Acts::GeometryContext geoCtx;
const Acts::MagneticFieldContext magCtx;

inline Propagator makePropagator(
    double bz, std::shared_ptr<const Acts::TrackingGeometry> geo) {
  auto magField = std::make_shared<MagneticField>(Acts::Vector3(0.0, 0.0, bz));
  Stepper stepper(std::move(magField));
  return Propagator(std::move(stepper), Acts::Navigator({std::move(geo)}));
}

}  // namespace

BOOST_AUTO_TEST_SUITE(PropagationSympyDenseConstant)

// check that the propagation is reversible and self-consistent

BOOST_DATA_TEST_CASE(ForwardBackward,
                     ds::phi* ds::thetaWithoutBeam* ds::absMomentum*
                         ds::chargeNonZero* ds::pathLength* ds::magneticField,
                     phi, theta, p, q, s, bz) {
  runForwardBackwardTest<Propagator>(
      makePropagator(bz, createDenseBlock(geoCtx)), geoCtx, magCtx,
      makeParametersCurvilinear(phi, theta, p, q), s, epsPos, epsDir, epsMom);
}

// check effects on the parameter covariance

BOOST_DATA_TEST_CASE(DenseTelescopeCovariance,
                     ds::absMomentum* ds::chargeNonZero, p, q) {
  const double bz = 0_T;
  const Acts::Material material = Acts::Test::makeLiquidArgon();
  const double thickness = 1_m;

  auto [geo, surfaces] = createDenseTelescope(geoCtx, material, thickness);

  auto propagator = makePropagator(bz, std::move(geo));

  runDenseForwardTest(propagator, geoCtx, magCtx, p, q, *surfaces.back(),
                      material, thickness);
}

BOOST_AUTO_TEST_SUITE_END()
