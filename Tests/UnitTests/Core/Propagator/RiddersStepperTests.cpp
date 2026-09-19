// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Direction.hpp"
#include "Acts/Definitions/Tolerance.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/BoundTrackParameters.hpp"
#include "Acts/EventData/ParticleHypothesis.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Propagator/ConstrainedStep.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/RiddersStepper.hpp"
#include "Acts/Surfaces/BoundaryTolerance.hpp"
#include "Acts/Surfaces/CurvilinearSurface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Utilities/Intersection.hpp"
#include "ActsTests/CommonHelpers/FloatComparisons.hpp"

#include <memory>

using namespace Acts;
using namespace Acts::UnitLiterals;

namespace ActsTests {

using Stepper = Experimental::RiddersStepper<EigenStepper<>>;

BOOST_AUTO_TEST_SUITE(PropagatorSuite)

// A transport anchors the state, so a second transport at the same point does
// not change the covariance again
BOOST_AUTO_TEST_CASE(ridders_stepper_repeated_transport) {
  const GeometryContext geoCtx = GeometryContext::dangerouslyDefaultConstruct();
  const MagneticFieldContext magCtx;

  Stepper stepper(std::make_shared<ConstantBField>(Vector3(0, 0, 2_T)));

  BoundMatrix cov = BoundMatrix::Identity();
  cov(eBoundQOverP, eBoundQOverP) = 1e-4;
  const auto start = BoundTrackParameters::createCurvilinear(
      Vector4::Zero(), 0.3, 1.2, 1. / 1_GeV, cov, ParticleHypothesis::pion());

  Stepper::Options options(geoCtx, magCtx);
  options.maxStepSize = 10_mm;
  auto state = stepper.makeState(options);
  stepper.initialize(state, start);

  for (int i = 0; i < 20; ++i) {
    BOOST_REQUIRE(stepper.step(state, Direction::Forward(), nullptr).ok());
  }

  const auto surface = CurvilinearSurface(stepper.position(state) +
                                              5_mm * stepper.direction(state),
                                          stepper.direction(state))
                           .planeSurface();

  // Step all states onto the surface, like the navigator does
  bool onSurface = false;
  for (int i = 0; i < 10 && !onSurface; ++i) {
    onSurface =
        stepper.updateSurfaceStatus(
            state, *surface, 0, Direction::Forward(),
            BoundaryTolerance::Infinite(), s_onSurfaceTolerance,
            ConstrainedStep::Type::Navigator) == IntersectionStatus::onSurface;
    if (!onSurface) {
      BOOST_REQUIRE(stepper.step(state, Direction::Forward(), nullptr).ok());
    }
  }
  BOOST_REQUIRE(onSurface);

  const auto firstJacobian = stepper.transportToBound(state, *surface);
  BOOST_REQUIRE(firstJacobian.ok());
  const BoundMatrix firstCovariance = stepper.covariance(state);
  BOOST_CHECK(!firstJacobian->isIdentity(1e-3));

  const auto secondJacobian = stepper.transportToBound(state, *surface);
  BOOST_REQUIRE(secondJacobian.ok());
  CHECK_CLOSE_ABS(*secondJacobian, BoundMatrix(BoundMatrix::Identity()), 1e-6);
  CHECK_CLOSE_COVARIANCE(stepper.covariance(state), firstCovariance, 1e-6);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
