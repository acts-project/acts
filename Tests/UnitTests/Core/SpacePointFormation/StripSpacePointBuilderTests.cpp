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
#include "Acts/SpacePointFormation/SpacePointFormationError.hpp"
#include "Acts/SpacePointFormation/StripSpacePointBuilder.hpp"
#include "ActsTests/CommonHelpers/FloatComparisons.hpp"

#include <cmath>

using namespace Acts;
using namespace Acts::UnitLiterals;

namespace ActsTests {

namespace {

StripSpacePointBuilder::StripEnds makeStrip(const Vector3& center,
                                            const Vector3& direction,
                                            double length) {
  const Vector3 half = 0.5 * length * direction.normalized();
  return {center + half, center - half};
}

}  // namespace

BOOST_AUTO_TEST_SUITE(StripSpacePointBuilderSuite)

/// Two orthogonal unit length strips with an analytically known answer.
BOOST_AUTO_TEST_CASE(CosmicOrthogonalUnitStrips) {
  const StripSpacePointBuilder::StripEnds strip1{Vector3(0, 0, 0),
                                                 Vector3(-1, 0, 0)};
  const StripSpacePointBuilder::StripEnds strip2{Vector3(-0.5, 0.5, 1),
                                                 Vector3(-0.5, -0.5, 1)};

  const StripSpacePointBuilder::CosmicOptions options;
  const auto result =
      StripSpacePointBuilder::computeCosmicSpacePoint(strip1, strip2, options);

  BOOST_REQUIRE(result.ok());
  // The second strip runs along y at x = -0.5, so the closest point on the
  // first strip is at x = -0.5.
  CHECK_CLOSE_ABS(*result, Vector3(-0.5, 0, 0), 1e-9);
}

/// The closest approach is defined by the connecting vector being perpendicular
/// to both strips. Checked on a realistic module pair: long strips with a small
/// stereo angle, which is where a missing normalisation shows up most strongly.
BOOST_AUTO_TEST_CASE(CosmicPerpendicularityRealisticStereoPair) {
  const double stereo = 26_mrad;
  const StripSpacePointBuilder::StripEnds strip1 =
      makeStrip(Vector3(0, 0, 0), Vector3(0, 1, 0), 50_mm);
  const StripSpacePointBuilder::StripEnds strip2 =
      makeStrip(Vector3(0.3_mm, 0, 1_mm),
                Vector3(std::sin(stereo), std::cos(stereo), 0), 50_mm);

  const StripSpacePointBuilder::CosmicOptions options;
  const auto onFirst =
      StripSpacePointBuilder::computeCosmicSpacePoint(strip1, strip2, options);
  const auto onSecond =
      StripSpacePointBuilder::computeCosmicSpacePoint(strip2, strip1, options);

  BOOST_REQUIRE(onFirst.ok());
  BOOST_REQUIRE(onSecond.ok());

  const Vector3 connection = *onFirst - *onSecond;
  const Vector3 dir1 = (strip1.top - strip1.bottom).normalized();
  const Vector3 dir2 = (strip2.top - strip2.bottom).normalized();

  CHECK_SMALL(connection.dot(dir1), 1e-9);
  CHECK_SMALL(connection.dot(dir2), 1e-9);
  // Both points have to stay on their strip
  BOOST_CHECK_LE(std::abs((*onFirst).y()), 25_mm);
  BOOST_CHECK_LE(std::abs((*onSecond).y()), 25_mm);
}

/// The result must scale with the strips, i.e. not depend on the length unit.
BOOST_AUTO_TEST_CASE(CosmicScaleInvariance) {
  const StripSpacePointBuilder::StripEnds strip1{Vector3(0, 0, 0),
                                                 Vector3(-1, 0, 0)};
  const StripSpacePointBuilder::StripEnds strip2{Vector3(-0.5, 0.5, 1),
                                                 Vector3(-0.5, -0.5, 1)};

  constexpr double scale = 50;
  const StripSpacePointBuilder::StripEnds scaled1{scale * strip1.top,
                                                  scale * strip1.bottom};
  const StripSpacePointBuilder::StripEnds scaled2{scale * strip2.top,
                                                  scale * strip2.bottom};

  const StripSpacePointBuilder::CosmicOptions options;
  const auto result =
      StripSpacePointBuilder::computeCosmicSpacePoint(strip1, strip2, options);
  const auto scaledResult = StripSpacePointBuilder::computeCosmicSpacePoint(
      scaled1, scaled2, options);

  BOOST_REQUIRE(result.ok());
  BOOST_REQUIRE(scaledResult.ok());
  CHECK_CLOSE_ABS(*scaledResult, (scale * (*result)).eval(), 1e-6);
}

/// Two strips crossing at a known point, separated only along z.
BOOST_AUTO_TEST_CASE(CosmicCrossingStrips) {
  const StripSpacePointBuilder::StripEnds strip1 =
      makeStrip(Vector3(0, 0, 0), Vector3(0, 1, 0), 50_mm);
  const StripSpacePointBuilder::StripEnds strip2 =
      makeStrip(Vector3(0, 0, 2_mm), Vector3(1, 0, 0), 50_mm);

  const StripSpacePointBuilder::CosmicOptions options;
  const auto onFirst =
      StripSpacePointBuilder::computeCosmicSpacePoint(strip1, strip2, options);
  const auto onSecond =
      StripSpacePointBuilder::computeCosmicSpacePoint(strip2, strip1, options);

  BOOST_REQUIRE(onFirst.ok());
  BOOST_REQUIRE(onSecond.ok());
  CHECK_CLOSE_ABS(*onFirst, Vector3(0, 0, 0), 1e-9);
  CHECK_CLOSE_ABS(*onSecond, Vector3(0, 0, 2_mm), 1e-9);
}

/// A crossing beyond the end of the first strip is not a valid space point.
BOOST_AUTO_TEST_CASE(CosmicCrossingOffFirstStripRejected) {
  const StripSpacePointBuilder::StripEnds strip1{Vector3(0, 0, 0),
                                                 Vector3(-1, 0, 0)};
  // The strips cross at x = 0.5, half a strip length past the top end
  const StripSpacePointBuilder::StripEnds strip2{Vector3(0.5, 0.5, 1),
                                                 Vector3(0.5, -0.5, 1)};

  const StripSpacePointBuilder::CosmicOptions options;
  const auto result =
      StripSpacePointBuilder::computeCosmicSpacePoint(strip1, strip2, options);

  BOOST_REQUIRE(!result.ok());
  BOOST_CHECK(result.error() == SpacePointFormationError::OutsideLimits);
}

/// The crossing has to lie on the second strip as well, even though the space
/// point is reported on the first one.
BOOST_AUTO_TEST_CASE(CosmicCrossingOffSecondStripRejected) {
  const StripSpacePointBuilder::StripEnds strip1 =
      makeStrip(Vector3(0, 0, 0), Vector3(1, 0, 0), 50_mm);
  // Crosses the first strip at its centre, but far off its own ends
  const StripSpacePointBuilder::StripEnds strip2 =
      makeStrip(Vector3(0, 100_mm, 1_mm), Vector3(0, 1, 0), 50_mm);

  const StripSpacePointBuilder::CosmicOptions options;
  const auto result =
      StripSpacePointBuilder::computeCosmicSpacePoint(strip1, strip2, options);

  BOOST_REQUIRE(!result.ok());
  BOOST_CHECK(result.error() == SpacePointFormationError::OutsideLimits);
}

/// The tolerance is a fraction of the strip length, so a crossing just past the
/// end is still accepted while a larger overshoot is not.
BOOST_AUTO_TEST_CASE(CosmicStripLengthTolerance) {
  const StripSpacePointBuilder::CosmicOptions options;
  const double halfLength = 25_mm;

  const StripSpacePointBuilder::StripEnds strip1 =
      makeStrip(Vector3(0, 0, 0), Vector3(1, 0, 0), 2 * halfLength);

  // Overshoots of half and twice the tolerance
  for (const double overshoot :
       {0.5 * options.stripLengthTolerance, 2 * options.stripLengthTolerance}) {
    const double x = (1 + overshoot) * halfLength;
    const StripSpacePointBuilder::StripEnds strip2 =
        makeStrip(Vector3(x, 0, 1_mm), Vector3(0, 1, 0), 50_mm);

    const auto result = StripSpacePointBuilder::computeCosmicSpacePoint(
        strip1, strip2, options);

    if (overshoot < options.stripLengthTolerance) {
      BOOST_REQUIRE(result.ok());
      CHECK_CLOSE_ABS(*result, Vector3(x, 0, 0), 1e-9);
    } else {
      BOOST_REQUIRE(!result.ok());
      BOOST_CHECK(result.error() == SpacePointFormationError::OutsideLimits);
    }
  }
}

/// A degenerate strip of zero length has no direction to cross with.
BOOST_AUTO_TEST_CASE(CosmicZeroLengthStripRejected) {
  const StripSpacePointBuilder::StripEnds strip1{Vector3(0, 0, 0),
                                                 Vector3(0, 0, 0)};
  const StripSpacePointBuilder::StripEnds strip2 =
      makeStrip(Vector3(0, 0, 1_mm), Vector3(0, 1, 0), 50_mm);

  const StripSpacePointBuilder::CosmicOptions options;
  const auto result =
      StripSpacePointBuilder::computeCosmicSpacePoint(strip1, strip2, options);

  BOOST_REQUIRE(!result.ok());
  BOOST_CHECK(result.error() ==
              SpacePointFormationError::CosmicToleranceNotMet);
}

/// Parallel strips have no well defined closest approach and are rejected,
/// independently of how long the strips are.
BOOST_AUTO_TEST_CASE(CosmicParallelStripsRejected) {
  const StripSpacePointBuilder::CosmicOptions options;

  for (const double length : {1_mm, 50_mm, 500_mm}) {
    const StripSpacePointBuilder::StripEnds strip1 =
        makeStrip(Vector3(0, 0, 0), Vector3(0, 1, 0), length);
    const StripSpacePointBuilder::StripEnds strip2 =
        makeStrip(Vector3(0, 0, 1_mm), Vector3(0, 1, 0), length);

    const auto result = StripSpacePointBuilder::computeCosmicSpacePoint(
        strip1, strip2, options);

    BOOST_REQUIRE(!result.ok());
    BOOST_CHECK(result.error() ==
                SpacePointFormationError::CosmicToleranceNotMet);
  }
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
