// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/SpacePointFormation/SpacePointFormationError.hpp"
#include "Acts/SpacePointFormation/StripSpacePointBuilder.hpp"
#include "ActsTests/CommonHelpers/FloatComparisons.hpp"

#include <cmath>

using namespace Acts;
using namespace Acts::StripSpacePointBuilder;

namespace {

/// Stereo angle between the two strip layers, as in a barrel strip module.
constexpr double stereo = 0.04;

/// Inner strip: at x = 100, along z, centred on the x axis.
StripEnds innerStrip(double halfLength) {
  return StripEnds{Vector3(100, 0, halfLength), Vector3(100, 0, -halfLength)};
}

/// Outer strip: at x = 110, tilted by the stereo angle, centred on
/// (110, y, z). With a vertex at the origin this gives
///   m = (10 / (11 * halfLengthInner)) * (z - y / tan(stereo))
///   n = -y / (halfLengthOuter * sin(stereo))
StripEnds outerStrip(double y, double z, double halfLength) {
  const Vector3 centre(110, y, z);
  const Vector3 dir(0, std::sin(stereo), std::cos(stereo));
  return StripEnds{centre + halfLength * dir, centre - halfLength * dir};
}

/// Place the outer strip so that the space point parameters take the requested
/// values. Inverts the two relations above.
StripEnds outerStripFor(double m, double n, double halfLengthInner,
                        double halfLengthOuter) {
  const double y = -n * halfLengthOuter * std::sin(stereo);
  const double z = m * (11. * halfLengthInner) / 10. + y / std::tan(stereo);
  return outerStrip(y, z, halfLengthOuter);
}

}  // namespace

BOOST_AUTO_TEST_SUITE(StripSpacePointBuilderTests)

// A track from the origin through the centre of both strips puts the space
// point at the centre of the inner strip.
BOOST_AUTO_TEST_CASE(ConstrainedCentralHit) {
  const StripEnds first = innerStrip(30);
  const StripEnds second = outerStripFor(0., 0., 30, 30);

  ConstrainedOptions options;
  const Result<Vector3> sp =
      computeConstrainedSpacePoint(first, second, options);

  BOOST_REQUIRE(sp.ok());
  CHECK_CLOSE_ABS(sp->x(), 100., 1e-6);
  CHECK_CLOSE_ABS(sp->z(), 0., 1e-6);
}

// Well inside both strips: accepted, and the space point sits at m along the
// inner strip.
BOOST_AUTO_TEST_CASE(ConstrainedInsideBothStrips) {
  const StripEnds first = innerStrip(30);
  const StripEnds second = outerStripFor(0.5, -0.3, 30, 30);

  ConstrainedOptions options;
  const Result<Vector3> sp =
      computeConstrainedSpacePoint(first, second, options);

  BOOST_REQUIRE(sp.ok());
  CHECK_CLOSE_ABS(sp->z(), 0.5 * 30., 1e-6);
}

// Far outside stays rejected however generous the gap tolerance is.
BOOST_AUTO_TEST_CASE(ConstrainedFarOutsideIsRejected) {
  const StripEnds first = innerStrip(30);
  const StripEnds second = outerStripFor(10., 0.2, 30, 30);

  ConstrainedOptions options;
  options.stripLengthGapTolerance = 5.;

  const Result<Vector3> sp =
      computeConstrainedSpacePoint(first, second, options);

  BOOST_CHECK(!sp.ok());
  BOOST_CHECK_EQUAL(sp.error(), SpacePointFormationError::OutsideRelaxedLimits);
}

// Just beyond the end of the inner strip, with the outer one well inside.
BOOST_AUTO_TEST_CASE(ConstrainedRecoversOneSidedOvershoot) {
  const StripEnds first = innerStrip(30);
  const StripEnds second = outerStripFor(1.05, -0.5, 30, 30);

  ConstrainedOptions options;
  options.stripLengthTolerance = 0.01;

  // No gap tolerance: correctly outside the limits.
  options.stripLengthGapTolerance = 0.;
  const Result<Vector3> strict =
      computeConstrainedSpacePoint(first, second, options);
  BOOST_CHECK(!strict.ok());
  BOOST_CHECK_EQUAL(strict.error(),
                    SpacePointFormationError::OutsideRelaxedLimits);

  // With a gap tolerance the pair has to be recovered onto the strip end.
  options.stripLengthGapTolerance = 5.;
  const Result<Vector3> relaxed =
      computeConstrainedSpacePoint(first, second, options);
  BOOST_REQUIRE(relaxed.ok());
  CHECK_CLOSE_ABS(relaxed->z(), 30., 1e-6);
}

// The gap tolerance is a length, so a short strip gets a proportionally larger
// tolerance on its parameter than a long one.
BOOST_AUTO_TEST_CASE(GapToleranceUsesEachStripLength) {
  const StripEnds first = innerStrip(50);
  const StripEnds second = outerStripFor(0., -1.05, 50, 5);

  ConstrainedOptions options;
  options.stripLengthTolerance = 0.01;
  options.stripLengthGapTolerance = 1.;

  // 1 mm is 1% of the 100 mm inner strip but 10% of the 10 mm outer strip
  const Result<Vector3> sp =
      computeConstrainedSpacePoint(first, second, options);
  BOOST_REQUIRE(sp.ok());

  // Too small for either strip: rejected.
  options.stripLengthGapTolerance = 0.1;
  const Result<Vector3> tight =
      computeConstrainedSpacePoint(first, second, options);
  BOOST_CHECK(!tight.ok());
  BOOST_CHECK_EQUAL(tight.error(),
                    SpacePointFormationError::OutsideRelaxedLimits);
}

// After the shift both parameters still have to be on their strips.
BOOST_AUTO_TEST_CASE(ConstrainedRecoveryStillChecksLimits) {
  const StripEnds first = innerStrip(30);
  const StripEnds second = outerStripFor(1.05, -1.05, 30, 30);

  ConstrainedOptions options;
  options.stripLengthTolerance = 0.01;
  options.stripLengthGapTolerance = 5.;

  const Result<Vector3> sp =
      computeConstrainedSpacePoint(first, second, options);
  BOOST_CHECK(!sp.ok());
  BOOST_CHECK_EQUAL(sp.error(), SpacePointFormationError::OutsideLimits);
}

BOOST_AUTO_TEST_SUITE_END()
