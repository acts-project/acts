// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/BoundTrackParameters.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"
#include "ActsTests/CommonHelpers/CylindricalTrackingGeometry.hpp"
#include "ActsTests/CommonHelpers/FloatComparisons.hpp"

#include <cstdint>
#include <optional>

using namespace Acts;
using namespace Acts::UnitLiterals;

namespace ActsTests {

namespace {

const auto gctx = GeometryContext::dangerouslyDefaultConstruct();
const MagneticFieldContext mctx;

using TestPropagator = Propagator<StraightLineStepper, Navigator>;

// The builder owns the detector elements, so it must outlive the geometry
CylindricalTrackingGeometry geometryBuilder(gctx);
const auto geometry = geometryBuilder();

TestPropagator makePropagator() {
  return TestPropagator(StraightLineStepper(), Navigator({geometry}));
}

// A straight track from the origin out of the barrel
const BoundTrackParameters start = BoundTrackParameters::createCurvilinear(
    Vector4::Zero(), Vector3::UnitX(), 1 / 1_GeV, std::nullopt,
    ParticleHypothesis::pion());

// The path limit is far outside the world
constexpr double farPathLimit = 50_m;

double propagatePathLength(const TestPropagator::Options<>& options) {
  auto result = makePropagator().propagate(start, options);
  BOOST_REQUIRE(result.ok());
  return result->pathLength;
}

}  // namespace

BOOST_AUTO_TEST_SUITE(PropagatorSuite)

BOOST_AUTO_TEST_CASE(DefaultActorPathLimit) {
  TestPropagator::Options<> options(gctx, mctx);
  options.pathLimit = 10_cm;

  CHECK_CLOSE_ABS(propagatePathLength(options), 10_cm, 1_um);
}

BOOST_AUTO_TEST_CASE(DefaultActorEndOfWorld) {
  TestPropagator::Options<> options(gctx, mctx);
  options.pathLimit = farPathLimit;

  // By default the propagation goes on after the end of the world
  CHECK_CLOSE_ABS(propagatePathLength(options), farPathLimit, 1_um);

  options.stopAtEndOfWorld = true;
  BOOST_CHECK_LT(propagatePathLength(options), farPathLimit);
}

BOOST_AUTO_TEST_CASE(DefaultActorVolumeConstraint) {
  TestPropagator::Options<> options(gctx, mctx);
  options.pathLimit = farPathLimit;
  options.stopAtEndOfWorld = true;
  const double endOfWorldPathLength = propagatePathLength(options);

  // Stop when the track leaves the volume it starts in
  const auto startVolume =
      geometry->resolveLowestTrackingVolume(gctx, start.position(gctx));
  BOOST_REQUIRE(startVolume.ok());
  options.constrainToVolumeIds = {
      static_cast<std::uint32_t>((*startVolume)->geometryId().volume())};

  BOOST_CHECK_LT(propagatePathLength(options), endOfWorldPathLength);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
