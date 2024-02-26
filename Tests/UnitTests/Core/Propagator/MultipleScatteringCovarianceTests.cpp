// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/ParticleHypothesis.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"
#include "Acts/Propagator/VoidNavigator.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/VectorHelpers.hpp"

#include <optional>

namespace bdata = boost::unit_test::data;

namespace Acts {
namespace Test {

using namespace Acts::UnitLiterals;

Acts::GeometryContext gctx;
Acts::MagneticFieldContext mctx;

auto makeEigenPropagator(const Vector3& bField) {
  return Propagator{EigenStepper<>{std::make_shared<ConstantBField>(bField)},
                    VoidNavigator{}};
}

auto makeStraightLinePropagator() {
  return Propagator{StraightLineStepper{}, VoidNavigator{}};
}

auto makeDist(double a, double b, int seed) {
  return bdata::random(
      (bdata::engine = std::mt19937{}, bdata::seed = seed,
       bdata::distribution = std::uniform_real_distribution<double>(a, b)));
}

const auto bFieldDist = makeDist(0, 3_T, 11);
const auto directionDist = makeDist(0, 1, 12);
const auto distanceDist = makeDist(100_mm, 1000_mm, 13);
const auto momentumDist = makeDist(0.1_GeV, 10_GeV, 14);
const auto chargeDist = makeDist(-1, 1, 15);
const auto theta0Dist = makeDist(10_mrad, 100_mrad, 16);

// Test the invariance of the angle covariance with and without a magnetic
// field.
//
// The test is performed by propagating a track with a given direction and
// momentum through a magnetic field and comparing the covariance of the
// direction at the end of the propagation with the covariance of the direction
// at the end of a straight line propagation.
BOOST_DATA_TEST_CASE(angle_cov_field_invariance,
                     (bFieldDist ^ bFieldDist ^ bFieldDist ^ directionDist ^
                      directionDist ^ directionDist ^ distanceDist ^
                      momentumDist ^ chargeDist ^ theta0Dist) ^
                         bdata::xrange(1),
                     Bx, By, Bz, Dx, Dy, Dz, s, p, q_, theta0, index) {
  (void)index;
  const Vector3 bField{Bx, By, Bz};
  const Vector3 direction = Vector3{Dx, Dy, Dz}.normalized();
  const ParticleHypothesis particle = ParticleHypothesis::pion();
  const double q = q_ > 0 ? 1 : -1;
  const double qOverP = particle.qOverP(p, q);

  PropagatorOptions<> options(gctx, mctx);
  options.direction = Direction::Forward;
  options.pathLimit = s;

  // Using optionals as there is no default constructor for BoundTrackParameters
  std::optional<BoundTrackParameters> eigenParameters;
  std::optional<BoundTrackParameters> straightLineBackwardParameters;
  std::optional<BoundTrackParameters> straightLineParameters;

  // Reference propagation
  {
    // sigmaPhi = theta0 / sin(theta)
    const double sigmaPhi =
        theta0 * (direction.norm() / VectorHelpers::perp(direction));
    const double sigmaTheta = theta0;

    // Initial track parameters
    auto startSurface =
        Surface::makeShared<PlaneSurface>(Vector3{0, 0, 0}, direction);
    auto startParameters = BoundTrackParameters(
        startSurface,
        (BoundVector() << 0_mm, 0_mm, VectorHelpers::phi(direction),
         VectorHelpers::theta(direction), qOverP, 0_ns)
            .finished(),
        (BoundVector() << 1_um, 1_um, sigmaPhi, sigmaTheta, 0.01 / 1_GeV, 1_ps)
            .finished()
            .array()
            .square()
            .matrix()
            .asDiagonal(),
        particle);

    auto eigenPropagator = makeEigenPropagator(bField);
    auto eigenRes = eigenPropagator.propagate(startParameters, options);
    BOOST_CHECK(eigenRes.ok());
    BOOST_CHECK(eigenRes->endParameters.has_value());
    BOOST_CHECK(eigenRes->endParameters.value().covariance().has_value());
    eigenParameters = eigenRes->endParameters.value();

    std::cout << startParameters << std::endl;
    std::cout << eigenParameters.value() << std::endl;
  }

  // Now we propagate backwards with a straight line propagator so we can later
  // propagate forwards again to compare the covariances.
  {
    PropagatorOptions<> straightLineBackwardOptions(gctx, mctx);
    straightLineBackwardOptions.direction = Direction::Backward;
    straightLineBackwardOptions.pathLimit = -s;

    BoundTrackParameters startParameters = eigenParameters.value();

    auto straightLinePropagator = makeStraightLinePropagator();
    auto straightLineBackwardRes = straightLinePropagator.propagate(
        startParameters, straightLineBackwardOptions);
    BOOST_CHECK(straightLineBackwardRes.ok());
    BOOST_CHECK(straightLineBackwardRes->endParameters.has_value());
    BOOST_CHECK(straightLineBackwardRes->endParameters.value()
                    .covariance()
                    .has_value());
    straightLineBackwardParameters =
        straightLineBackwardRes->endParameters.value();
  }

  // Straight line propagation
  {
    // sigmaPhi = theta0 / sin(theta)
    const double sigmaPhi =
        theta0 * (eigenParameters->direction().norm() /
                  VectorHelpers::perp(eigenParameters->direction()));
    const double sigmaTheta = theta0;

    BoundTrackParameters startParameters =
        straightLineBackwardParameters.value();
    startParameters.covariance() =
        (BoundVector() << 1_um, 1_um, sigmaPhi, sigmaTheta, 1 / 1_MeV, 1_ps)
            .finished()
            .array()
            .square()
            .matrix()
            .asDiagonal();

    auto straightLinePropagator = makeStraightLinePropagator();
    auto straightLineForwardRes =
        straightLinePropagator.propagate(startParameters, options);
    BOOST_CHECK(straightLineForwardRes.ok());
    BOOST_CHECK(straightLineForwardRes->endParameters.has_value());
    BOOST_CHECK(
        straightLineForwardRes->endParameters.value().covariance().has_value());
    straightLineParameters = straightLineForwardRes->endParameters.value();
  }

  // Optional was just a hack and should always be set
  BOOST_CHECK(eigenParameters.has_value());
  BOOST_CHECK(straightLineParameters.has_value());

  // Check if we get the same parameters
  {
    BoundVector exp = eigenParameters->parameters();
    BoundVector obs = straightLineParameters->parameters();
    CHECK_CLOSE_ABS(obs, exp, 1e-7);
  }

  // Check if we get the same direction covariance
  {
    SquareMatrix2 exp =
        eigenParameters->covariance().value().block<2, 2>(eBoundPhi, eBoundPhi);
    SquareMatrix2 obs =
        straightLineParameters->covariance().value().block<2, 2>(eBoundPhi,
                                                                 eBoundPhi);
    CHECK_CLOSE_ABS(obs, exp, 1e-7);
  }

  // Check if we get the same position covariance
  {
    SquareMatrix2 exp = eigenParameters->covariance().value().block<2, 2>(
        eBoundLoc0, eBoundLoc0);
    SquareMatrix2 obs =
        straightLineParameters->covariance().value().block<2, 2>(eBoundLoc0,
                                                                 eBoundLoc0);
    CHECK_CLOSE_ABS(obs, exp, 1e-7);
  }
}

}  // namespace Test
}  // namespace Acts
