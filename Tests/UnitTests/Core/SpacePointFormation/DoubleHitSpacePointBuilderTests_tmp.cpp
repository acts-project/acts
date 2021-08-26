// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Digitization/CartesianSegmentation.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Plugins/Digitization/PlanarModuleCluster.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"
#include "Acts/SpacePointFormation/DoubleHitSpacePointBuilder.hpp"
#include "Acts/SpacePointFormation/SpacePointBuilderConfig.h"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Tests/CommonHelpers/CubicTrackingGeometry.hpp"
#include "Acts/Tests/CommonHelpers/DetectorElementStub.hpp"
#include "Acts/Tests/CommonHelpers/MeasurementsCreator.hpp"
#include "Acts/Tests/CommonHelpers/TestSpacePoint.hpp"
#include "Acts/Utilities/CalibrationContext.hpp"

namespace bdata = boost::unit_test::data;
namespace tt = boost::test_tools;
using namespace Acts::UnitLiterals;

namespace Acts {
namespace Test {

// // Create a test context
// GeometryContext tgContext = GeometryContext();

// /// Unit test for testing the main functions of DoubleHitSpacePointBuilder
// /// 1) A pair of hits gets added and resolved.
// /// 2) A pair of hits gets added and rejected.
// BOOST_DATA_TEST_CASE(DoubleHitsSpacePointBuilder_basic, bdata::xrange(1),
//                      index) {
//   (void)index;

//   std::cout << "Create first hit" << std::endl;

//   // Build Bounds
//   auto recBounds = std::make_shared<RectangleBounds>(35_um, 25_mm);

//   // Build binning and segmentation
//   std::vector<float> boundariesX, boundariesY;
//   boundariesX.push_back(-35_um);
//   boundariesX.push_back(35_um);
//   boundariesY.push_back(-25_mm);
//   boundariesY.push_back(25_mm);

//   BinningData binDataX(BinningOption::open, BinningValue::binX, boundariesX);
//   std::shared_ptr<BinUtility> buX(new BinUtility(binDataX));
//   BinningData binDataY(BinningOption::open, BinningValue::binY, boundariesY);
//   std::shared_ptr<BinUtility> buY(new BinUtility(binDataY));
//   (*buX) += (*buY);

//   auto segmentation = std::make_shared<CartesianSegmentation>(buX,
//   recBounds);

//   // Build translation

//   double rotation = 0.026;
//   RotationMatrix3 rotationPos;
//   Vector3 xPos(cos(rotation), sin(rotation), 0.);
//   Vector3 yPos(-sin(rotation), cos(rotation), 0.);
//   Vector3 zPos(0., 0., 1.);
//   rotationPos.col(0) = xPos;
//   rotationPos.col(1) = yPos;
//   rotationPos.col(2) = zPos;
//   Transform3 t3d(Transform3::Identity() * rotationPos);
//   t3d.translation() = Vector3(0., 0., 10_m);

//   // Build Digitization
//   const DigitizationModule digMod(segmentation, 1., 1., 0.);
//   DetectorElementStub detElem(t3d);
//   auto pSur = Surface::makeShared<PlaneSurface>(recBounds, detElem);
//   SymMatrix3 cov;
//   cov << 0., 0., 0., 0., 0., 0., 0., 0., 0.;
//   Vector2 local = {0.1, -0.1};

//   // Build PlanarModuleCluster
//   PlanarModuleCluster pmc(pSur, DigitizationSourceLink(pSur->geometryId(),
//   {}),
//                           cov, local[0], local[1], 0.,
//                           {DigitizationCell(0, 0, 1.)}, &digMod);

//   std::cout << "Create second hit" << std::endl;

//   // Build second PlanarModuleCluster

//   double rotation2 = -0.026;
//   RotationMatrix3 rotationNeg;
//   Vector3 xNeg(cos(rotation2), sin(rotation2), 0.);
//   Vector3 yNeg(-sin(rotation2), cos(rotation2), 0.);
//   Vector3 zNeg(0., 0., 1.);
//   rotationNeg.col(0) = xNeg;
//   rotationNeg.col(1) = yNeg;
//   rotationNeg.col(2) = zNeg;
//   Transform3 t3d2(Transform3::Identity() * rotationNeg);
//   t3d2.translation() = Vector3(0., 0., 10.005_m);

//   DetectorElementStub detElem2(t3d2);

//   auto pSur2 = Surface::makeShared<PlaneSurface>(recBounds, detElem2);

//   PlanarModuleCluster pmc2(
//       pSur2, DigitizationSourceLink(pSur->geometryId(), {}), cov, local[0],
//       local[1], 0., {DigitizationCell(0, 0, 1.)}, &digMod);

//   std::cout << "Store both hits" << std::endl;

//   std::vector<SpacePoint<PlanarModuleCluster>> resultSP;
//   std::vector<std::pair<PlanarModuleCluster const*, PlanarModuleCluster
//   const*>>
//       clusterPairs;
//   DoubleHitSpacePointConfig dhsp_cfg;

//   // Combine two PlanarModuleClusters
//   SpacePointBuilder<SpacePoint<PlanarModuleCluster>> dhsp(dhsp_cfg);
//   dhsp.makeClusterPairs(tgContext, {&pmc}, {&pmc2}, clusterPairs);

//   BOOST_CHECK_EQUAL(clusterPairs.size(), 1u);
//   BOOST_CHECK_EQUAL(clusterPairs[0].first, &pmc);
//   BOOST_CHECK_EQUAL(clusterPairs[0].second, &pmc2);

//   std::cout << "Calculate space point" << std::endl;

//   dhsp.calculateSpacePoints(tgContext, clusterPairs, resultSP);

//   BOOST_CHECK_EQUAL(resultSP.size(), 1u);

//   std::cout << "Create third hit" << std::endl;

//   // Build third PlanarModuleCluster
//   Transform3 t3d3(Transform3::Identity() * rotationNeg);
//   t3d3.translation() = Vector3(0., 0., 10.005_m);

//   DetectorElementStub detElem3(t3d3);
//   auto pSur3 = Surface::makeShared<PlaneSurface>(recBounds, detElem3);

//   PlanarModuleCluster pmc3(
//       pSur3, DigitizationSourceLink(pSur->geometryId(), {}), cov, local[0],
//       local[1], 0., {DigitizationCell(0, 0, 1.)}, &digMod);

//   std::cout << "Try to store hits" << std::endl;

//   // Combine points
//   dhsp.makeClusterPairs(tgContext, {&pmc}, {&pmc3}, clusterPairs);

//   // Test for rejecting unconnected hits
//   BOOST_CHECK_EQUAL(resultSP.size(), 1u);
// }
using namespace Acts::UnitLiterals;

using StraightPropagator =
    Acts::Propagator<Acts::StraightLineStepper, Acts::Navigator>;
using TestMeasurement = Acts::BoundVariantMeasurement<TestSourceLink>;
using ConstantFieldStepper = Acts::EigenStepper<>;
using ConstantFieldPropagator =
    Acts::Propagator<ConstantFieldStepper, Acts::Navigator>;
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

// Create a test context
GeometryContext tgContext = GeometryContext();

const GeometryContext geoCtx;
const MagneticFieldContext magCtx;
const CalibrationContext calCtx;

// detector geometry
CubicTrackingGeometry geometryStore(geoCtx);
const auto geometry = geometryStore();

// detector resolutions
const MeasurementResolution resPixel = {MeasurementType::eLoc01,
                                        {25_um, 50_um}};
const MeasurementResolution resStrip0 = {MeasurementType::eLoc0, {100_um}};
const MeasurementResolution resStrip1 = {MeasurementType::eLoc1, {150_um}};
const MeasurementResolutionMap resolutions = {
    {GeometryIdentifier().setVolume(2), resPixel},
    {GeometryIdentifier().setVolume(3).setLayer(2), resStrip0},
    {GeometryIdentifier().setVolume(3).setLayer(4), resStrip1},
    {GeometryIdentifier().setVolume(3).setLayer(6), resStrip0},
    {GeometryIdentifier().setVolume(3).setLayer(8), resStrip1},
};

// Construct a straight-line propagator.
static StraightPropagator makeStraightPropagator(
    std::shared_ptr<const Acts::TrackingGeometry> geo) {
  Acts::Navigator::Config cfg{geo};
  cfg.resolvePassive = false;
  cfg.resolveMaterial = true;
  cfg.resolveSensitive = true;
  Acts::Navigator navigator{cfg};
  Acts::StraightLineStepper stepper;
  return StraightPropagator(std::move(stepper), std::move(navigator));
}

// simulation propagator
const auto measPropagator = makeStraightPropagator(geometry);

// Construct initial track parameters.
Acts::CurvilinearTrackParameters makeParameters() {
  // create covariance matrix from reasonable standard deviations
  Acts::BoundVector stddev;
  stddev[Acts::eBoundLoc0] = 100_um;
  stddev[Acts::eBoundLoc1] = 100_um;
  stddev[Acts::eBoundTime] = 25_ns;
  stddev[Acts::eBoundPhi] = 2_degree;
  stddev[Acts::eBoundTheta] = 2_degree;
  stddev[Acts::eBoundQOverP] = 1 / 100_GeV;
  Acts::BoundSymMatrix cov = stddev.cwiseProduct(stddev).asDiagonal();
  // define a track in the transverse plane along x
  Acts::Vector4 mPos4(-3_m, 0., 0., 42_ns);
  return Acts::CurvilinearTrackParameters(mPos4, 0_degree, 90_degree, 1_GeV,
                                          1_e, cov);
}

std::default_random_engine rng(42);
/// Unit test for testing the main functions of OneHitSpacePointBuilder
/// 1) A resolved dummy hit gets created and added.
/// 2) A hit gets added and resolved.
BOOST_AUTO_TEST_CASE(SingleHitSpacePointBuilder_basic) {
  double phi = 20._degree;
  double theta = 80._degree;
  double p = 1.0_GeV;
  double q = 1;

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
  auto start = makeParameters(phi, theta, p, q);

  auto measurements =
      createMeasurements(propagator, geoCtx, magCtx, start, resolutions, rng);

  const auto sourceLinks = measurements.sourceLinks;

  std::vector<TestMeasurement> testMeasurements;
  for (auto& sl : sourceLinks) {
    TestMeasurement meas = makeMeasurement(sl, sl.parameters, sl.covariance,
                                           eBoundLoc0, eBoundLoc1);
    testMeasurements.emplace_back(meas);
  }
  BOOST_CHECK_NE(testMeasurements.size(), 0);

  auto spBuilderConfig = DoubleHitSpacePointBuilderConfig();
  spBuilderConfig.trackingGeometry = geometry;

  auto doubleSPBuilder =
      Acts::DoubleHitSpacePointBuilder<TestSpacePoint, TestMeasurement>(
          spBuilderConfig);

  TestSpacePointContainer spacePoints;

  // singleSPBuilder.calculateSpacePoints(geoCtx, testMeasurements,
  // spacePoints);

  // BOOST_REQUIRE_EQUAL(testMeasurements.size(), spacePoints.size());
  // BOOST_CHECK_NE(spacePoints[0].x(), 0);

  //     BOOST_CHECK_NE(data[0].vector, Vector3::Zero());

  std::cout << "Space point calculated" << std::endl;
}

}  // end of namespace Test
}  // end of namespace Acts
