// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/tools/old/interface.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/detail/TestSourceLink.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Geometry/LayerCreator.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Seeding/EstimateTrackParamsFromSeed.hpp"
#include "Acts/Seeding/TrackParamsEstimationError.hpp"
#include "Acts/Seeding/detail/CircleFit.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsTests/CommonHelpers/CylindricalTrackingGeometry.hpp"
#include "ActsTests/CommonHelpers/FloatComparisons.hpp"
#include "ActsTests/CommonHelpers/MeasurementsCreator.hpp"
#include "ActsTests/CommonHelpers/TestLogger.hpp"

#include <array>
#include <cmath>
#include <memory>
#include <optional>
#include <random>
#include <set>
#include <span>
#include <utility>
#include <vector>

using namespace Acts;
using namespace Acts::UnitLiterals;

namespace ActsTests {

using ConstantFieldStepper = EigenStepper<>;
using ConstantFieldPropagator = Propagator<ConstantFieldStepper, Navigator>;

const auto geoCtx = GeometryContext::dangerouslyDefaultConstruct();
const MagneticFieldContext magCtx;

// detector geometry
CylindricalTrackingGeometry geometryStore(geoCtx);
const auto geometry = geometryStore();

// Two dimensional measurement with zero resolution
const MeasurementResolutionMap resolutions = {
    {GeometryIdentifier(),
     MeasurementResolution{MeasurementType::eLoc01, {0, 0}}}};

// Construct initial track parameters.
BoundTrackParameters makeParameters(double phi, double theta, double p,
                                    double q) {
  // create covariance matrix from reasonable standard deviations
  BoundVector stddev;
  stddev[eBoundLoc0] = 100_um;
  stddev[eBoundLoc1] = 100_um;
  stddev[eBoundTime] = 25_ns;
  stddev[eBoundPhi] = 2_degree;
  stddev[eBoundTheta] = 2_degree;
  stddev[eBoundQOverP] = 1 / 100_GeV;
  BoundMatrix cov = stddev.cwiseProduct(stddev).asDiagonal();
  // Let the particle starts from the origin
  Vector4 mPos4(0., 0., 0., 0.);
  return BoundTrackParameters::createCurvilinear(
      mPos4, phi, theta, q / p, cov, ParticleHypothesis::pionLike(std::abs(q)));
}

std::default_random_engine rng(42);

BOOST_AUTO_TEST_SUITE(SeedingSuite)

BOOST_AUTO_TEST_CASE(trackparameters_estimation_test) {
  // Construct a propagator with the cylinderal geometry and a constant magnetic
  // field along z
  Navigator navigator({
      geometry,
      true,  // sensitive
      true,  // material
      false  // passive
  });
  const Vector3 bField(0, 0, 2._T);
  auto field = std::make_shared<ConstantBField>(bField);
  ConstantFieldStepper stepper(std::move(field));

  ConstantFieldPropagator propagator(std::move(stepper), std::move(navigator));

  std::array<double, 2> pArray = {0.5_GeV, 1.0_GeV};
  std::array<double, 3> phiArray = {20._degree, 0._degree - 20._degree};
  std::array<double, 3> thetaArray = {80._degree, 90.0_degree, 100._degree};
  std::array<double, 2> qArray = {1, -1};

  auto logger = getTestLogger("estimateTrackParamsFromSeed", Logging::INFO);

  for (const auto& p : pArray) {
    for (const auto& phi : phiArray) {
      for (const auto& theta : thetaArray) {
        for (const auto& q : qArray) {
          BOOST_TEST_INFO("Test track with p = " << p << ", phi = " << phi
                                                 << ", theta = " << theta
                                                 << ", q = " << q);
          auto start = makeParameters(phi, theta, p, q);
          auto measurements = createMeasurements(propagator, geoCtx, magCtx,
                                                 start, resolutions, rng);

          // Create space points from different detector layers
          std::vector<Vector3> spacePoints;
          std::set<GeometryIdentifier::Value> usedLayers;
          const Surface* bottomSurface = nullptr;
          for (const auto& sl : measurements.sourceLinks) {
            const auto geoId = sl.m_geometryId;
            const auto& layer = geoId.layer();
            // Avoid to use space point from the same layers
            if (const auto it = usedLayers.find(layer);
                it != usedLayers.end()) {
              continue;
            }
            const auto surface = geometry->findSurface(geoId);
            const auto& localPos = sl.parameters;
            const Vector3 globalFakeMom(1, 1, 1);
            const Vector3 globalPos =
                surface->localToGlobal(geoCtx, localPos, globalFakeMom);
            spacePoints.push_back(globalPos);
            usedLayers.insert(layer);
            if (spacePoints.size() == 1) {
              bottomSurface = surface;
            }
          }

          // Check if there are at least 3 space points
          if (spacePoints.size() < 3) {
            BOOST_TEST_WARN("Number of space points less than 3.");
            continue;
          }

          // The truth track parameters at the bottom space point
          const auto& expBoundParams = measurements.truthParameters[0];
          BOOST_TEST_INFO(
              "The truth track parameters at the bottom space point: \n"
              << expBoundParams.transpose());

          // Test the free track parameters estimator
          FreeVector estFreeParams = estimateTrackParamsFromSeed(
              spacePoints[0], 0, spacePoints[1], spacePoints[2], bField);
          BOOST_CHECK(!estFreeParams.hasNaN());

          // Test the bound track parameters estimator
          auto estBoundParamsResult = estimateTrackParamsFromSeed(
              geoCtx, *bottomSurface, spacePoints[0], 0, spacePoints[1],
              spacePoints[2], bField);
          BOOST_CHECK(estBoundParamsResult.ok());
          const auto& estBoundParams = estBoundParamsResult.value();
          BOOST_TEST_INFO(
              "The estimated full track parameters at the bottom space point: "
              "\n"
              << estBoundParams.transpose());

          CHECK_CLOSE_ABS(estBoundParams[eBoundLoc0],
                          expBoundParams[eBoundLoc0], 1e-5);
          CHECK_CLOSE_ABS(estBoundParams[eBoundLoc1],
                          expBoundParams[eBoundLoc1], 1e-5);
          // @todo Understand why the estimated phi has a limited precision
          CHECK_CLOSE_ABS(estBoundParams[eBoundPhi], expBoundParams[eBoundPhi],
                          1e-1);
          CHECK_CLOSE_ABS(estBoundParams[eBoundTheta],
                          expBoundParams[eBoundTheta], 1e-2);
          CHECK_CLOSE_ABS(estBoundParams[eBoundQOverP],
                          expBoundParams[eBoundQOverP], 1e-2);
          // time is not estimated so we check if it is default zero
          CHECK_CLOSE_ABS(estBoundParams[eBoundTime], 0, 1e-6);
        }
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(trackparm_estimate_aligined) {
  Vector3 sp0{-72.775, -0.325, -615.6};
  Vector3 sp1{-84.325, -0.325, -715.6};
  Vector3 sp2{-98.175, -0.325, -835.6};
  Vector3 bField{0, 0, 0.000899377};

  FreeVector params = estimateTrackParamsFromSeed(sp0, 0, sp1, sp2, bField);
  BOOST_CHECK_EQUAL(params[eFreeQOverP], 0);
}

BOOST_AUTO_TEST_CASE(trackparm_estimate_degenerate) {
  const Vector3 bField{0, 0, 2_T};

  // bottom and middle share a transverse position, so the estimation frame has
  // no x axis
  {
    const Vector3 sp0{10, 0, 0};
    const Vector3 sp1{10, 0, 100};
    const Vector3 sp2{10, 0, 200};

    BOOST_CHECK(
        !estimateTrackParamsFromSeed(sp0, 0, sp1, sp2, bField).allFinite());
  }

  // middle and top share a transverse position, so the conformal mapping maps
  // them onto each other
  {
    const Vector3 sp0{10, 0, 0};
    const Vector3 sp1{20, 5, 0};
    const Vector3 sp2{20, 5, 50};

    BOOST_CHECK(
        !estimateTrackParamsFromSeed(sp0, 0, sp1, sp2, bField).allFinite());
  }
}

// Multi space point estimator against truth, and against the closed-form
// three-point estimator at N=3.
BOOST_AUTO_TEST_CASE(spacepoint_estimator_vs_truth) {
  Navigator navigator({
      geometry,
      true,   // sensitive
      true,   // material
      false,  // passive
  });
  const Vector3 bField(0, 0, 2._T);
  auto field = std::make_shared<ConstantBField>(bField);
  ConstantFieldStepper stepper(std::move(field));
  ConstantFieldPropagator propagator(std::move(stepper), std::move(navigator));

  std::array<double, 2> pArray = {0.5_GeV, 1.0_GeV};
  std::array<double, 3> phiArray = {20._degree, 0._degree, -20._degree};
  std::array<double, 3> thetaArray = {80._degree, 90._degree, 100._degree};
  std::array<double, 2> qArray = {1, -1};

  for (const auto& p : pArray) {
    for (const auto& phi : phiArray) {
      for (const auto& theta : thetaArray) {
        for (const auto& q : qArray) {
          BOOST_TEST_INFO("Test track with p = " << p << ", phi = " << phi
                                                 << ", theta = " << theta
                                                 << ", q = " << q);
          auto start = makeParameters(phi, theta, p, q);
          auto measurements = createMeasurements(propagator, geoCtx, magCtx,
                                                 start, resolutions, rng);

          // Collect one space point per detector layer, in track order.
          std::vector<Vector3> spacePoints;
          std::set<GeometryIdentifier::Value> usedLayers;
          const Surface* bottomSurface = nullptr;
          for (const auto& sl : measurements.sourceLinks) {
            const auto geoId = sl.m_geometryId;
            const auto& layer = geoId.layer();
            if (const auto it = usedLayers.find(layer);
                it != usedLayers.end()) {
              continue;
            }
            const auto surface = geometry->findSurface(geoId);
            const Vector3 globalFakeMom(1, 1, 1);
            const Vector3 globalPos =
                surface->localToGlobal(geoCtx, sl.parameters, globalFakeMom);
            spacePoints.push_back(globalPos);
            usedLayers.insert(layer);
            if (spacePoints.size() == 1) {
              bottomSurface = surface;
            }
          }

          if (spacePoints.size() < 3) {
            BOOST_TEST_WARN("Number of space points less than 3.");
            continue;
          }

          const auto& expBoundParams = measurements.truthParameters[0];

          auto freeResult =
              estimateTrackParamsFromSpacePoints(spacePoints, bField);
          BOOST_REQUIRE(freeResult.ok());
          auto boundResult = transformFreeToBoundParameters(
              *freeResult, *bottomSurface, geoCtx);
          BOOST_REQUIRE(boundResult.ok());
          const auto& estBoundParams = boundResult.value();

          CHECK_CLOSE_ABS(estBoundParams[eBoundLoc0],
                          expBoundParams[eBoundLoc0], 1e-4);
          CHECK_CLOSE_ABS(estBoundParams[eBoundLoc1],
                          expBoundParams[eBoundLoc1], 1e-4);
          CHECK_CLOSE_ABS(estBoundParams[eBoundPhi], expBoundParams[eBoundPhi],
                          1e-1);
          CHECK_CLOSE_ABS(estBoundParams[eBoundTheta],
                          expBoundParams[eBoundTheta], 1e-2);
          CHECK_CLOSE_ABS(estBoundParams[eBoundQOverP],
                          expBoundParams[eBoundQOverP], 1e-2);

          // At N=3 both estimators must agree, which pins the sign and
          // momentum conventions.
          auto free3 = estimateTrackParamsFromSpacePoints(
              std::span<const Vector3>(spacePoints).first(3), bField, 0);
          BOOST_CHECK(free3.ok());
          const FreeVector ref = estimateTrackParamsFromSeed(
              spacePoints[0], 0, spacePoints[1], spacePoints[2], bField);
          for (int i = 0; i < 3; ++i) {
            CHECK_CLOSE_ABS((*free3)[eFreeDir0 + i], ref[eFreeDir0 + i], 1e-6);
          }
          CHECK_CLOSE_ABS((*free3)[eFreeQOverP], ref[eFreeQOverP],
                          1e-6 * std::abs(ref[eFreeQOverP]) + 1e-9);
        }
      }
    }
  }
}

// More than half a turn, to exercise the phi unwrapping and the R-Z fit.
BOOST_AUTO_TEST_CASE(spacepoint_estimator_large_arc) {
  const double bMag = 2._T;
  const Vector3 bField(0, 0, bMag);
  const double radius = 300;  // mm
  const double theta = 70._degree;
  const double lambda = 1 / std::tan(theta);  // dz/ds
  const double alpha0 = 0.3;
  const double dAlpha = 40._degree;  // 6 points -> 200 degrees
  const std::size_t nsp = 6;

  auto runWith = [&](std::size_t refineIter) {
    std::vector<Vector3> spacePoints;
    for (std::size_t i = 0; i < nsp; ++i) {
      const double a = alpha0 + i * dAlpha;
      const double s = radius * (a - alpha0);
      spacePoints.emplace_back(radius * std::cos(a), radius * std::sin(a),
                               100 + lambda * s);
    }
    const Vector3 first = spacePoints.front();

    auto result =
        estimateTrackParamsFromSpacePoints(spacePoints, bField, 0, refineIter);
    BOOST_CHECK(result.ok());
    const FreeVector& fp = *result;
    BOOST_CHECK(!fp.hasNaN());

    for (int i = 0; i < 3; ++i) {
      CHECK_CLOSE_ABS(fp[eFreePos0 + i], first[i], 1e-6);
    }

    // Tangent at the first point: transverse (-sin a0, cos a0), z slope lambda.
    Vector3 dir(-std::sin(alpha0), std::cos(alpha0), lambda);
    dir.normalize();
    for (int i = 0; i < 3; ++i) {
      CHECK_CLOSE_ABS(fp[eFreeDir0 + i], dir[i], 1e-6);
    }

    // |q/p| = sin(theta) / (R * |B|).
    const double expAbsQOverP = std::sin(theta) / (radius * bMag);
    CHECK_CLOSE_ABS(std::abs(fp[eFreeQOverP]), expAbsQOverP,
                    1e-6 * expAbsQOverP);
  };

  runWith(3);
  runWith(0);
}

// The reference point can be any space point: the fitted helix is the same and
// only the point it is evaluated at moves.
BOOST_AUTO_TEST_CASE(spacepoint_estimator_reference_index) {
  const double bMag = 2._T;
  const Vector3 bField(0, 0, bMag);
  const double radius = 300;  // mm
  const double theta = 70._degree;
  const double lambda = 1 / std::tan(theta);  // dz/ds
  const double alpha0 = 0.3;
  const double dAlpha = 40._degree;  // 6 points -> 200 degrees
  const std::size_t nsp = 6;

  std::vector<Vector3> spacePoints;
  for (std::size_t i = 0; i < nsp; ++i) {
    const double a = alpha0 + i * dAlpha;
    const double s = radius * (a - alpha0);
    spacePoints.emplace_back(radius * std::cos(a), radius * std::sin(a),
                             100 + lambda * s);
  }

  const auto atFirst = estimateTrackParamsFromSpacePoints(spacePoints, bField);
  BOOST_REQUIRE(atFirst.ok());

  for (std::size_t ref = 0; ref < nsp; ++ref) {
    BOOST_TEST_CONTEXT("referenceIndex = " << ref) {
      auto result = estimateTrackParamsFromSpacePoints(spacePoints, bField, 0,
                                                       0, {}, ref);
      BOOST_REQUIRE(result.ok());
      const FreeVector& fp = *result;

      // The parameters sit on the chosen space point.
      for (int i = 0; i < 3; ++i) {
        CHECK_CLOSE_ABS(fp[eFreePos0 + i], spacePoints[ref][i], 1e-6);
      }

      // Tangent of the truth helix there: transverse (-sin a, cos a), z slope
      // lambda.
      const double a = alpha0 + ref * dAlpha;
      Vector3 dir(-std::sin(a), std::cos(a), lambda);
      dir.normalize();
      for (int i = 0; i < 3; ++i) {
        CHECK_CLOSE_ABS(fp[eFreeDir0 + i], dir[i], 1e-6);
      }

      // The helix itself does not depend on the reference, so charge and
      // momentum are identical to the estimate at the first point.
      BOOST_CHECK_EQUAL(fp[eFreeQOverP], (*atFirst)[eFreeQOverP]);
    }
  }

  // The default is the first space point.
  const auto atZero =
      estimateTrackParamsFromSpacePoints(spacePoints, bField, 0, 0, {}, 0);
  BOOST_REQUIRE(atZero.ok());
  BOOST_CHECK_EQUAL(*atZero, *atFirst);
}

// Zero field: a straight-line fit with vanishing q/p.
BOOST_AUTO_TEST_CASE(spacepoint_estimator_zero_field) {
  const Vector3 bField = Vector3::Zero();
  const Vector3 dir = Vector3(1, 2, 3).normalized();
  const Vector3 p0(10, -5, 7);
  std::vector<Vector3> spacePoints;
  for (int i = 0; i < 5; ++i) {
    spacePoints.push_back(p0 + i * 20 * dir);
  }

  auto result = estimateTrackParamsFromSpacePoints(spacePoints, bField, 0);
  BOOST_CHECK(result.ok());
  const FreeVector& fp = *result;
  BOOST_CHECK_EQUAL(fp[eFreeQOverP], 0);
  for (int i = 0; i < 3; ++i) {
    CHECK_CLOSE_ABS(fp[eFreeDir0 + i], dir[i], 1e-9);
    CHECK_CLOSE_ABS(fp[eFreePos0 + i], p0[i], 1e-9);
  }
}

// Fewer than three space points is an error.
BOOST_AUTO_TEST_CASE(spacepoint_estimator_not_enough_points) {
  const std::array<Vector3, 2> spacePoints{Vector3(1, 0, 0), Vector3(2, 0, 0)};
  auto result =
      estimateTrackParamsFromSpacePoints(spacePoints, Vector3(0, 0, 2._T), 0);
  BOOST_CHECK(!result.ok());
  BOOST_CHECK(result.error() ==
              TrackParamsEstimationError::NotEnoughSpacePoints);
}

// Coincident space points must report a degenerate fit.
BOOST_AUTO_TEST_CASE(spacepoint_estimator_degenerate) {
  const std::array<Vector3, 4> spacePoints{Vector3(3, 4, 5), Vector3(3, 4, 5),
                                           Vector3(3, 4, 5), Vector3(3, 4, 5)};
  auto result =
      estimateTrackParamsFromSpacePoints(spacePoints, Vector3(0, 0, 2._T), 0);
  BOOST_CHECK(!result.ok());
  BOOST_CHECK(result.error() == TrackParamsEstimationError::DegenerateFit);
}

// The Taubin fit recovers a known circle and rejects collinear input.
BOOST_AUTO_TEST_CASE(taubin_circle_fit) {
  const Vector2 center(3.5, -2);
  const double radius = 12;
  std::vector<Vector3> points;
  for (int i = 0; i < 7; ++i) {
    const double a = 0.2 + i * 0.35;
    const Vector2 xy = center + radius * Vector2(std::cos(a), std::sin(a));
    points.emplace_back(xy.x(), xy.y(), 0);
  }
  const auto circle = Acts::detail::fitCircleTaubin(points);
  BOOST_REQUIRE(circle.has_value());
  CHECK_CLOSE_ABS(circle->center.x(), center.x(), 1e-6);
  CHECK_CLOSE_ABS(circle->center.y(), center.y(), 1e-6);
  CHECK_CLOSE_ABS(circle->radius, radius, 1e-6);

  // Collinear points do not define a finite circle.
  std::vector<Vector3> line;
  for (int i = 0; i < 5; ++i) {
    line.emplace_back(2 * i, 1 + 6 * i, 0);
  }
  BOOST_CHECK(!Acts::detail::fitCircleTaubin(line).has_value());
}

// A near-zero weight on a displaced point recovers the true circle.
BOOST_AUTO_TEST_CASE(taubin_circle_fit_weighted) {
  const Vector2 center(3.5, -2);
  const double radius = 12;
  std::vector<Vector3> points;
  std::vector<double> weights;
  for (int i = 0; i < 7; ++i) {
    const double a = 0.2 + i * 0.35;
    const Vector2 xy = center + radius * Vector2(std::cos(a), std::sin(a));
    points.emplace_back(xy.x(), xy.y(), 0);
    weights.push_back(1);
  }
  // Displace the last point radially outward and give it a tiny weight.
  const double aOut = 0.2 + 6 * 0.35;
  const Vector2 xyOut =
      center + (radius + 5) * Vector2(std::cos(aOut), std::sin(aOut));
  points.back() = Vector3(xyOut.x(), xyOut.y(), 0);
  weights.back() = 1e-6;

  const auto unweighted = Acts::detail::fitCircleTaubin(points);
  BOOST_REQUIRE(unweighted.has_value());
  BOOST_CHECK(std::abs(unweighted->radius - radius) > 0.1);

  const auto weighted = Acts::detail::fitCircleTaubin(points, weights);
  BOOST_REQUIRE(weighted.has_value());
  CHECK_CLOSE_ABS(weighted->center.x(), center.x(), 1e-3);
  CHECK_CLOSE_ABS(weighted->center.y(), center.y(), 1e-3);
  CHECK_CLOSE_ABS(weighted->radius, radius, 1e-3);
}

// The same through the estimator: a helix with one point displaced in z, where
// uniform weights bias theta and a small weight does not.
BOOST_AUTO_TEST_CASE(spacepoint_estimator_weighted) {
  const double bMag = 2._T;
  const Vector3 bField(0, 0, bMag);
  const double radius = 300;  // mm
  const double theta = 70._degree;
  const double lambda = 1 / std::tan(theta);  // dz/ds
  const double alpha0 = 0.3;
  const double dAlpha = 40._degree;
  const std::size_t nsp = 6;

  std::vector<Vector3> spacePoints;
  for (std::size_t i = 0; i < nsp; ++i) {
    const double a = alpha0 + i * dAlpha;
    const double s = radius * (a - alpha0);
    spacePoints.emplace_back(radius * std::cos(a), radius * std::sin(a),
                             100 + lambda * s);
  }
  // Displace one point in z only (transverse position untouched, so the circle
  // fit stays on the true circle) and down-weight it.
  spacePoints[3].z() += 50;
  std::vector<double> weights(nsp, 1);
  weights[3] = 1e-6;

  const std::span<const Vector3> sp(spacePoints);

  auto biased = estimateTrackParamsFromSpacePoints(sp, bField, 0, 0);
  BOOST_CHECK(biased.ok());
  auto corrected =
      estimateTrackParamsFromSpacePoints(sp, bField, 0, 0, weights);
  BOOST_CHECK(corrected.ok());

  Vector3 dir(-std::sin(alpha0), std::cos(alpha0), lambda);
  dir.normalize();

  // The weighted fit recovers the true direction; the unweighted one is pulled
  // off by the displaced point.
  double biasedErr = 0, correctedErr = 0;
  for (int i = 0; i < 3; ++i) {
    biasedErr += std::abs((*biased)[eFreeDir0 + i] - dir[i]);
    correctedErr += std::abs((*corrected)[eFreeDir0 + i] - dir[i]);
  }
  BOOST_CHECK(correctedErr < biasedErr);
  CHECK_CLOSE_ABS((*corrected)[eFreeDir0 + 2], dir[2], 1e-3);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
