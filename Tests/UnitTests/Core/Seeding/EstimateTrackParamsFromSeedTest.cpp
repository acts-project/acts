// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/EventData/Measurement.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Seeding/EstimateTrackParamsFromSeed.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Tests/CommonHelpers/CylindricalTrackingGeometry.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Tests/CommonHelpers/MeasurementsCreator.hpp"
#include "Acts/Utilities/CalibrationContext.hpp"

#include <algorithm>
#include <array>
#include <vector>

#include "SpacePoint.hpp"

namespace {

using namespace Acts;
using namespace Acts::Test;
using namespace Acts::UnitLiterals;

using ConstantFieldStepper = Acts::EigenStepper<>;
using ConstantFieldPropagator =
    Acts::Propagator<ConstantFieldStepper, Acts::Navigator>;

const GeometryContext geoCtx;
const MagneticFieldContext magCtx;
const CalibrationContext calCtx;

// detector geometry
CylindricalTrackingGeometry geometryStore(geoCtx);
const auto geometry = geometryStore();

// Two dimensional measurement with zero resolution
const MeasurementResolutionMap resolutions = {
    {GeometryIdentifier(),
     MeasurementResolution{MeasurementType::eLoc01, {0, 0}}}};

// Construct initial track parameters.
CurvilinearTrackParameters makeParameters(double phi, double theta, double p,
                                          double q) {
  // create covariance matrix from reasonable standard deviations
  Acts::BoundVector stddev;
  stddev[Acts::eBoundLoc0] = 100_um;
  stddev[Acts::eBoundLoc1] = 100_um;
  stddev[Acts::eBoundTime] = 25_ns;
  stddev[Acts::eBoundPhi] = 2_degree;
  stddev[Acts::eBoundTheta] = 2_degree;
  stddev[Acts::eBoundQOverP] = 1 / 100_GeV;
  BoundSymMatrix cov = stddev.cwiseProduct(stddev).asDiagonal();
  // Let the particle starts from the origin
  Vector4 mPos4(0., 0., 0., 0.);
  return CurvilinearTrackParameters(mPos4, phi, theta, p, q, cov);
}

std::default_random_engine rng(42);

}  // namespace

BOOST_AUTO_TEST_CASE(trackparameters_estimation_test) {
  // Construct a propagator with the cylinderal geometry and a constant magnetic
  // field along z
  Acts::Navigator navigator({
      geometry,
      true,  // sensitive
      true,  // material
      false  // passive
  });
  auto field =
      std::make_shared<Acts::ConstantBField>(Acts::Vector3(0.0, 0.0, 2._T));
  ConstantFieldStepper stepper(std::move(field));

  ConstantFieldPropagator propagator(std::move(stepper), std::move(navigator));

  std::array<double, 2> pArray = {0.5_GeV, 1.0_GeV};
  std::array<double, 3> phiArray = {20._degree, 0._degree - 20._degree};
  std::array<double, 3> thetaArray = {80._degree, 90.0_degree, 100._degree};
  std::array<double, 2> qArray = {1, -1};

  auto logger = Acts::getDefaultLogger("estimateTrackParamsFromSeed",
                                       Acts::Logging::INFO);

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
            const auto geoId = sl.geometryId();
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
                                  static_cast<int>(geoId.layer()), 0., 0.});
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
          // The curvature of track projection on the transverse plane in unit
          // of 1/mm
          double rho = expParams[eBoundQOverP] * 0.3 * 2. / UnitConstants::m;

          // The space point pointers
          std::array<const SpacePoint*, 3> spacePointPtrs{};
          std::transform(spacePoints.begin(), std::next(spacePoints.begin(), 3),
                         spacePointPtrs.begin(),
                         [](const auto& sp) { return &sp.second; });

          // Test the partial track parameters estimator
          auto partialParamsOpt = estimateTrackParamsFromSeed(
              spacePointPtrs.begin(), spacePointPtrs.end(), *logger);
          BOOST_REQUIRE(partialParamsOpt.has_value());
          const auto& estPartialParams = partialParamsOpt.value();
          BOOST_TEST_INFO(
              "The estimated track parameters at the transverse plane: \n"
              << estPartialParams.transpose());

          // The particle starting position is (0, 0, 0). Hence, d0 is zero; the
          // phi at the point of cloest approach is exactly the phi of the truth
          // particle
          CHECK_CLOSE_ABS(estPartialParams[eBoundLoc0], 0., 1e-5);
          CHECK_CLOSE_ABS(estPartialParams[eBoundPhi], phi, 1e-5);
          CHECK_CLOSE_ABS(estPartialParams[eBoundQOverP], rho, 1e-4);
          // The loc1, theta and time are set to zero in the estimator
          CHECK_CLOSE_ABS(estPartialParams[eBoundLoc1], 0., 1e-10);
          CHECK_CLOSE_ABS(estPartialParams[eBoundTheta], 0., 1e-10);
          CHECK_CLOSE_ABS(estPartialParams[eBoundTime], 0., 1e-10);

          // Test the full track parameters estimator
          auto fullParamsOpt = estimateTrackParamsFromSeed(
              geoCtx, spacePointPtrs.begin(), spacePointPtrs.end(),
              *bottomSurface, Vector3(0, 0, 2._T), 0.1_T, *logger);
          BOOST_REQUIRE(fullParamsOpt.has_value());
          const auto& estFullParams = fullParamsOpt.value();
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
          CHECK_CLOSE_ABS(estFullParams[eBoundTime], expParams[eBoundTime], 1.);
        }
      }
    }
  }
}
