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

SquareMatrix2 makeAngleCov(double sigmaPhi, double sigmaTheta) {
  return (Vector2() << sigmaPhi, sigmaTheta)
      .finished()
      .array()
      .square()
      .matrix()
      .asDiagonal();
}

SquareMatrix2 makeMscAngleCov(double theta0, const Vector3& direction) {
  // sigmaPhi = theta0 / sin(theta)
  const double sigmaPhi =
      theta0 * (direction.norm() / VectorHelpers::perp(direction));
  const double sigmaTheta = theta0;
  return makeAngleCov(sigmaPhi, sigmaTheta);
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
                         bdata::xrange(1000),
                     Bx, By, Bz, Dx, Dy, Dz, s, p, q_, theta0, index) {
  (void)index;
  const Vector3 bField{Bx, By, Bz};
  const Vector3 direction = Vector3{Dx, Dy, Dz}.normalized();
  const ParticleHypothesis particle = ParticleHypothesis::pion();
  const double q = q_ > 0 ? 1 : -1;
  const double qOverP = particle.qOverP(p, q);

  std::cout << "Test " << index << " with B = " << bField.transpose() / 1_T
            << " direction = " << direction.transpose() << " s = " << s / 1_mm
            << " p = " << p / 1_GeV << " q = " << q
            << " theta0 = " << theta0 / 1_mrad << std::endl;

  std::cout << bField.normalized().dot(direction) << std::endl;

  PropagatorOptions<> options(gctx, mctx);
  options.direction = Direction::Forward;
  options.pathLimit = s;

  // Using optionals as there is no default constructor for BoundTrackParameters
  std::optional<BoundTrackParameters> referenceParameters;
  std::optional<BoundTrackParameters> otherParameters;

  // Initial track parameters
  auto startSurface =
      Surface::makeShared<PlaneSurface>(Vector3{0, 0, 0}, direction);
  auto startParameters = BoundTrackParameters(
      startSurface,
      (BoundVector() << 0_mm, 0_mm, VectorHelpers::phi(direction),
       VectorHelpers::theta(direction), qOverP, 0_ns)
          .finished(),
      (BoundVector() << 1_um, 1_um, 10_mrad, 10_mrad, 0.01 / 1_GeV, 1_ps)
          .finished()
          .array()
          .square()
          .matrix()
          .asDiagonal(),
      particle);

  // Reference propagation
  {
    auto referenceStartParameters = startParameters;
    referenceStartParameters.covariance().value().block<2, 2>(
        eBoundPhi, eBoundPhi) += makeMscAngleCov(theta0, direction);

    auto propagator = makeEigenPropagator(bField);
    auto res = propagator.propagate(referenceStartParameters, options);
    BOOST_CHECK(res.ok());
    BOOST_CHECK(res->endParameters.has_value());
    BOOST_CHECK(res->endParameters.value().covariance().has_value());
    referenceParameters = res->endParameters.value();
  }

  // Propagation without MSC and correction at the end
  {
    auto propagator = makeEigenPropagator(bField);
    auto res = propagator.propagate(startParameters, options);
    BOOST_CHECK(res.ok());
    BOOST_CHECK(res->endParameters.has_value());
    BOOST_CHECK(res->endParameters.value().covariance().has_value());
    otherParameters = res->endParameters.value();

    // Correct the covariance

    SquareMatrix2 locCov =
        (Vector2() << std::sin(theta0) * s, std::sin(theta0) * s)
            .finished()
            .array()
            .square()
            .matrix()
            .asDiagonal();

    otherParameters->covariance().value().block<2, 2>(eBoundPhi, eBoundPhi) +=
        makeMscAngleCov(theta0, otherParameters->direction());
    otherParameters->covariance().value().block<2, 2>(eBoundLoc0, eBoundLoc0) +=
        locCov;
  }

  // Check if we get the same parameters
  {
    BoundVector exp = referenceParameters->parameters();
    BoundVector obs = otherParameters->parameters();
    CHECK_CLOSE_ABS(obs, exp, 1e-12);
  }

  // Check if we get the same direction covariance
  {
    SquareMatrix2 ref = referenceParameters->covariance().value().block<2, 2>(
        eBoundPhi, eBoundPhi);
    SquareMatrix2 obs =
        otherParameters->covariance().value().block<2, 2>(eBoundPhi, eBoundPhi);
    // diagonal should be close
    CHECK_CLOSE_REL(obs.diagonal(), ref.diagonal(), 1e-1);
    CHECK_CLOSE_COVARIANCE(obs, ref, 1e-1);
  }

  // Check if we get the same position covariance
  {
    SquareMatrix2 ref = referenceParameters->covariance().value().block<2, 2>(
        eBoundLoc0, eBoundLoc0);
    SquareMatrix2 obs = otherParameters->covariance().value().block<2, 2>(
        eBoundLoc0, eBoundLoc0);
    // diagonal should be close
    CHECK_CLOSE_REL(obs.diagonal(), ref.diagonal(), 1e-1);
    CHECK_CLOSE_COVARIANCE(obs, ref, 1e-1);
  }
}

}  // namespace Test
}  // namespace Acts
