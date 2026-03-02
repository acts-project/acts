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
#include "Acts/EventData/TrackParameters.hpp"
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
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsTests/CommonHelpers/CylindricalTrackingGeometry.hpp"
#include "ActsTests/CommonHelpers/FloatComparisons.hpp"
#include "ActsTests/CommonHelpers/MeasurementsCreator.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <iterator>
#include <map>
#include <memory>
#include <optional>
#include <random>
#include <utility>
#include <vector>

#include "SpacePoint.hpp"

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

  auto logger = getDefaultLogger("estimateTrackParamsFromSeed", Logging::INFO);

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
          std::map<GeometryIdentifier::Value, SpacePoint> spacePoints;
          const Surface* bottomSurface = nullptr;
          for (const auto& sl : measurements.sourceLinks) {
            const auto geoId = sl.m_geometryId;
            const auto& layer = geoId.layer();
            auto it = spacePoints.find(layer);
            // Avoid to use space point from the same layers
            if (it != spacePoints.end()) {
              continue;
            }
            const auto surface = geometry->findSurface(geoId);
            const auto& localPos = sl.parameters;
            Vector3 globalFakeMom(1, 1, 1);
            Vector3 globalPos =
                surface->localToGlobal(geoCtx, localPos, globalFakeMom);
            // Create a space point (varianceR and varianceZ are lazily set to
            // zero since they are not important for the test)
            float r = std::hypot(globalPos.x(), globalPos.y());
            spacePoints.emplace(
                layer, SpacePoint{static_cast<float>(globalPos.x()),
                                  static_cast<float>(globalPos.y()),
                                  static_cast<float>(globalPos.z()), r,
                                  static_cast<int>(geoId.layer()), 0., 0.,
                                  std::nullopt, std::nullopt});
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
          const auto& expParams = measurements.truthParameters[0];
          BOOST_TEST_INFO(
              "The truth track parameters at the bottom space point: \n"
              << expParams.transpose());

          // The space point pointers
          std::array<const SpacePoint*, 3> spacePointPtrs{};
          std::transform(spacePoints.begin(), std::next(spacePoints.begin(), 3),
                         spacePointPtrs.begin(),
                         [](const auto& sp) { return &sp.second; });

          // Test the free track parameters estimator
          FreeVector estFreeParams =
              estimateTrackParamsFromSeed(spacePointPtrs, bField);
          BOOST_CHECK(!estFreeParams.hasNaN());

          // Test the bound track parameters estimator
          auto estFullParamsResult = estimateTrackParamsFromSeed(
              geoCtx, spacePointPtrs, *bottomSurface, bField);
          BOOST_CHECK(estFullParamsResult.ok());
          const auto& estFullParams = estFullParamsResult.value();
          BOOST_TEST_INFO(
              "The estimated full track parameters at the bottom space point: "
              "\n"
              << estFullParams.transpose());

          CHECK_CLOSE_ABS(estFullParams[eBoundLoc0], expParams[eBoundLoc0],
                          1e-5);
          CHECK_CLOSE_ABS(estFullParams[eBoundLoc1], expParams[eBoundLoc1],
                          1e-5);
          // @todo Understand why the estimated phi has a limited precision
          CHECK_CLOSE_ABS(estFullParams[eBoundPhi], expParams[eBoundPhi], 1e-1);
          CHECK_CLOSE_ABS(estFullParams[eBoundTheta], expParams[eBoundTheta],
                          1e-2);
          CHECK_CLOSE_ABS(estFullParams[eBoundQOverP], expParams[eBoundQOverP],
                          1e-2);
          // time is not estimated so we check if it is default zero
          CHECK_CLOSE_ABS(estFullParams[eBoundTime], 0, 1e-6);
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

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
