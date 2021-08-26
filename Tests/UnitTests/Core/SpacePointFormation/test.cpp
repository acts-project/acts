// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/Measurement.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Tests/CommonHelpers/GenerateParameters.hpp"
#include "Acts/Tests/CommonHelpers/TestSourceLink.hpp"
//#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Digitization/CartesianSegmentation.hpp"
#include "Acts/EventData/Measurement.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/SpacePointFormation/SingleHitSpacePointBuilder.hpp"
#include "Acts/SpacePointFormation/SpacePointBuilderConfig.h"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Tests/CommonHelpers/CubicTrackingGeometry.hpp"
#include "Acts/Tests/CommonHelpers/DetectorElementStub.hpp"
#include "Acts/Tests/CommonHelpers/TestSpacePoint.hpp"

#include <limits>
#include <random>
#include <tuple>
#include <vector>

using namespace Acts;
using namespace Acts::Test;
using namespace Acts::UnitLiterals;
using SourceLink = Acts::Test::TestSourceLink;
using TestMeasurement = Acts::BoundVariantMeasurement<TestSourceLink>;
namespace bd = boost::unit_test::data;

namespace {
constexpr BoundIndices boundIndices[] = {
    eBoundLoc0, eBoundLoc1, eBoundTime, eBoundPhi, eBoundTheta, eBoundQOverP,
};
constexpr FreeIndices freeIndices[] = {
    eFreePos0, eFreePos1, eFreePos2, eFreeTime,
    eFreeDir0, eFreeDir1, eFreeDir2, eFreeQOverP,
};
const TestSourceLink source;
// fix seed for reproducible tests
std::default_random_engine rng(123);
}  // namespace

// the underlying subspace impl
// separate unit test. here we only test concrete extreme cases and
// measurement-specific functionality.

// Create a test context
GeometryContext tgContext = GeometryContext();

const GeometryContext geoCtx;
// const MagneticFieldContext magCtx;
// const CalibrationContext calCtx;

// detector geometry
CubicTrackingGeometry cGeometry(geoCtx);
const auto geometry = cGeometry();
// geometry.boundarySurfaces();
// const auto geometry = geometryStore(geoCtx);

BOOST_AUTO_TEST_SUITE(EventDataMeasurement)

// BOOST_DATA_TEST_CASE(FixedBoundOne, bd::make(boundIndices), index) {
//   auto [params, cov] = generateParametersCovariance<ActsScalar, 1u>(rng);
//   auto meas = makeMeasurement(source, params, cov, index);

//   BOOST_CHECK_EQUAL(meas.size(), 1);
//   for (auto i : boundIndices) {
//     if (i == index) {
//       BOOST_CHECK(meas.contains(i));
//     } else {
//       BOOST_CHECK(not meas.contains(i));
//     }
//   }
//   BOOST_CHECK_EQUAL(meas.parameters(), params);
//   BOOST_CHECK_EQUAL(meas.covariance(), cov);
//   BOOST_CHECK_EQUAL(meas.sourceLink(), source);
// }

// BOOST_AUTO_TEST_CASE(FixedBoundAll) {
//   auto [params, cov] = generateBoundParametersCovariance(rng);
//   auto meas = makeMeasurement(source, params, cov, eBoundLoc0, eBoundLoc1,
//                               eBoundPhi, eBoundTheta, eBoundQOverP,
//                               eBoundTime);

//   BOOST_CHECK_EQUAL(meas.size(), eBoundSize);
//   for (auto i : boundIndices) {
//     BOOST_CHECK(meas.contains(i));
//   }
//   BOOST_CHECK_EQUAL(meas.parameters(), params);
//   BOOST_CHECK_EQUAL(meas.covariance(), cov);
//   BOOST_CHECK_EQUAL(meas.sourceLink(), source);
// }

namespace {
// example data for phi residual tests. each entry contains
//
//     measured, reference, expected residual
//
const std::vector<std::tuple<double, double, double>> kPhiDataset = {
    // measurement and reference in bounds and close
    {0.5, 0.75, -0.25},
    // measurement and reference in bounds but at different edges
    {0.25, 2 * M_PI - 0.25, 0.5},
    {2 * M_PI - 0.125, 0.125, -0.25},
    // measurement in bounds, reference ouf-of-bounds, both near lower edge
    {0.25, -0.25, 0.5},
    // measurement in bounds, reference ouf-of-bounds, both near upper edge
    {2 * M_PI - 0.25, 2 * M_PI + 0.25, -0.5},
    // measurement out-of-bounds, reference in bounds, both near lower edge
    {-0.25, 0.25, -0.5},
    // measurement out-of-bounds, reference in bounds, both near upper edge
    {2 * M_PI + 0.25, 2 * M_PI - 0.25, 0.5},
};
}  // namespace

BOOST_DATA_TEST_CASE(BoundResidualsPhi, bd::make(kPhiDataset), phiMea, phiRef,
                     phiRes) {
  using MeasurementVector = Acts::ActsVector<1>;
  using MeasurementCovariance = Acts::ActsSymMatrix<1>;

  // prepare measurement
  MeasurementVector params = MeasurementVector::Zero();
  MeasurementCovariance cov = MeasurementCovariance::Zero();
  params[0] = phiMea;
  auto measurement = makeMeasurement(source, params, cov, eBoundPhi);
  // prepare reference parameters
  Acts::BoundVector reference = Acts::BoundVector::Zero();
  reference[eBoundPhi] = phiRef;

  // compute and check residual
  auto res = measurement.residuals(reference);
  CHECK_CLOSE_ABS(res[0], phiRes, std::numeric_limits<ActsScalar>::epsilon());
}

BOOST_DATA_TEST_CASE(SingleHitSpacePointBuilder_test, bd::xrange(1), index) {
  // (void)index;
  std::cout << "test" << std::endl;
  // auto [params_gen, cov_gen] = generateParametersCovariance<ActsScalar,
  // 2u>(rng);
  /// auto meas = makeMeasurement(source, params, cov, 0);
  // auto meas = makeMeasurement(source, params, cov, Acts::eBoundLoc0);

  // auto [params1, cov1] = generateParametersCovariance<ActsScalar, 2u>(rng);
  // std::cout << " " << std::endl;
  // std::cout << params1 << std::endl;
  // std::cout << " " << std::endl;

  // Build bounds
  std::shared_ptr<const RectangleBounds> recBounds(
      new RectangleBounds(35_um, 25_mm));

  // Build binning and segmentation
  std::vector<float> boundariesX, boundariesY;
  boundariesX.push_back(-35_um);
  boundariesX.push_back(35_um);
  boundariesY.push_back(-25_mm);
  boundariesY.push_back(25_mm);

  BinningData binDataX(BinningOption::open, BinningValue::binX, boundariesX);
  std::shared_ptr<BinUtility> buX(new BinUtility(binDataX));
  BinningData binDataY(BinningOption::open, BinningValue::binY, boundariesY);
  std::shared_ptr<BinUtility> buY(new BinUtility(binDataY));
  (*buX) += (*buY);

  std::shared_ptr<const Segmentation> segmentation(
      new CartesianSegmentation(buX, recBounds));

  // Build translation

  double rotation = 0.026_rad;
  RotationMatrix3 rotationPos;
  Vector3 xPos(cos(rotation), sin(rotation), 0.);
  Vector3 yPos(-sin(rotation), cos(rotation), 0.);
  Vector3 zPos(0., 0., 1.);
  rotationPos.col(0) = xPos;
  rotationPos.col(1) = yPos;
  rotationPos.col(2) = zPos;
  Transform3 t3d(Transform3::Identity() * rotationPos);
  t3d.translation() = Vector3(0., 0., 10_m);

  // Build Digitization
  // const DigitizationModule digMod(segmentation, 1., 1., 0.);
  DetectorElementStub detElem(t3d);
  auto pSur = Surface::makeShared<PlaneSurface>(recBounds, detElem);
  pSur->assignGeometryId(GeometryIdentifier(13u));
  //  SymMatrix3 cov;
  //  cov << 0., 0., 0., 0., 0., 0., 0., 0., 0.;
  Vector2 local = {0.1, -0.1};

  BoundVector vec = BoundVector::Zero();
  vec[eBoundLoc0] = local[0];
  vec[eBoundLoc1] = local[1];

  const auto param = BoundTrackParameters(pSur, vec);
  // std::vector<TestMeasurement> testMeasurements;
  // const TestSourceLink slink = TestSourceLink();
  auto param_test = Eigen::Matrix<ActsScalar, 2, 1>(local[0], local[1]);
  auto cov_test = Eigen::Matrix<ActsScalar, 2, 2>();
  cov_test << 0., 0., 0., 0.;
  // param_test << 0., 1.;
  // param_test[eBoundLoc0] = local[0];
  // param_test[eBoundLoc1] = local[1];
  const auto geoId = pSur->geometryId();

  TestSourceLink sl = TestSourceLink(Acts::eBoundLoc0, Acts::eBoundLoc1,
                                     param_test, cov_test, geoId);

  auto meas = makeMeasurement(sl, param_test, cov_test, Acts::eBoundLoc0,
                              Acts::eBoundLoc1);
  // std::cout << meas.size() << std::endl;
  std::vector<TestMeasurement> testMeasurements = {meas};

  auto param0 = sl.parameters;
  auto gid = sl.geoId;
  auto gid2 = pSur->geometryId();
  std::cout << "geoId2 :" << gid2 << std::endl;
  std::cout << param << std::endl;

  auto spBuilderConfig = SingleHitSpacePointBuilderConfig();
  spBuilderConfig.trackingGeometry = geometry;

  auto singleSPBuilder =
      Acts::SingleHitSpacePointBuilder<TestSpacePoint, TestSourceLink>(
          spBuilderConfig);
  TestSpacePointContainer spacePoints;

  singleSPBuilder.calculateSpacePoints(geoCtx, testMeasurements, spacePoints);

  BOOST_REQUIRE_EQUAL(testMeasurements.size(), spacePoints.size());
  BOOST_CHECK_NE(spacePoints[0].x(), 0);
}

// BOOST_DATA_TEST_CASE(FixedFreeOne, bd::make(freeIndices), index) {
//   //std::cout << index << std::endl;
//   auto [params, cov] = generateParametersCovariance<ActsScalar, 1u>(rng);
//   //std::cout << params << std::endl;
//   //std::cout << std::endl;
//   auto meas = makeMeasurement(source, params, cov, index);

//   BOOST_CHECK_EQUAL(meas.size(), 1);
//   for (auto i : freeIndices) {
//     if (i == index) {
//       BOOST_CHECK(meas.contains(i));
//     } else {
//       BOOST_CHECK(not meas.contains(i));
//     }
//   }
//   BOOST_CHECK_EQUAL(meas.parameters(), params);
//   BOOST_CHECK_EQUAL(meas.covariance(), cov);
//   BOOST_CHECK_EQUAL(meas.sourceLink(), source);

//   // all free parameters are unrestricted and we know the expected residual.
//   constexpr auto tol = std::numeric_limits<ActsScalar>::epsilon();
//   auto [ref, refCov] = generateFreeParametersCovariance(rng);
//   auto res = meas.residuals(ref);
//   CHECK_CLOSE_ABS(res[0], params[0] - ref[index], tol);
// }

// BOOST_AUTO_TEST_CASE(FixedFreeAll) {
//   auto [params, cov] = generateFreeParametersCovariance(rng);
//   auto meas =
//       makeMeasurement(source, params, cov, eFreePos0, eFreePos1, eFreePos2,
//                       eFreeTime, eFreeDir0, eFreeDir1, eFreeDir2,
//                       eFreeQOverP);

//   BOOST_CHECK_EQUAL(meas.size(), eFreeSize);
//   for (auto i : freeIndices) {
//     BOOST_CHECK(meas.contains(i));
//   }
//   BOOST_CHECK_EQUAL(meas.parameters(), params);
//   BOOST_CHECK_EQUAL(meas.covariance(), cov);
//   BOOST_CHECK_EQUAL(meas.sourceLink(), source);

//   // all free parameters are unrestricted and we know the expected residual.
//   constexpr auto tol = std::numeric_limits<ActsScalar>::epsilon();
//   auto [ref, refCov] = generateFreeParametersCovariance(rng);
//   CHECK_CLOSE_ABS(meas.residuals(ref), params - ref, tol);
// }

// BOOST_AUTO_TEST_CASE(VariantBound) {
//   // generate w/ a single parameter
//   auto [par1, cov1] = generateParametersCovariance<ActsScalar, 1u>(rng);
//   BoundVariantMeasurement<SourceLink> meas =
//       makeMeasurement(source, par1, cov1, eBoundTheta);
//   std::visit(
//       [](const auto& m) {
//         BOOST_CHECK_EQUAL(m.size(), 1);
//         BOOST_CHECK(not m.contains(eBoundLoc0));
//         BOOST_CHECK(not m.contains(eBoundLoc1));
//         BOOST_CHECK(not m.contains(eBoundTime));
//         BOOST_CHECK(not m.contains(eBoundPhi));
//         BOOST_CHECK(m.contains(eBoundTheta));
//         BOOST_CHECK(not m.contains(eBoundQOverP));
//       },
//       meas);

//   // reassign w/ all parameters
//   auto [parN, covN] = generateBoundParametersCovariance(rng);
//   meas = makeMeasurement(source, parN, covN, eBoundLoc0, eBoundLoc1,
//   eBoundPhi,
//                          eBoundTheta, eBoundQOverP, eBoundTime);
//   std::visit(
//       [](const auto& m) {
//         BOOST_CHECK_EQUAL(m.size(), eBoundSize);
//         BOOST_CHECK(m.contains(eBoundLoc0));
//         BOOST_CHECK(m.contains(eBoundLoc1));
//         BOOST_CHECK(m.contains(eBoundTime));
//         BOOST_CHECK(m.contains(eBoundPhi));
//         BOOST_CHECK(m.contains(eBoundTheta));
//         BOOST_CHECK(m.contains(eBoundQOverP));
//       },
//       meas);
// }

// BOOST_AUTO_TEST_CASE(VariantFree) {
//   // generate w/ two parameters
//   auto [par2, cov2] = generateParametersCovariance<ActsScalar, 2u>(rng);
//   FreeVariantMeasurement<SourceLink> meas =
//       makeMeasurement(source, par2, cov2, eFreePos2, eFreeTime);
//   std::visit(
//       [](const auto& m) {
//         BOOST_CHECK_EQUAL(m.size(), 2);
//         BOOST_CHECK(not m.contains(eFreePos0));
//         BOOST_CHECK(not m.contains(eFreePos1));
//         BOOST_CHECK(m.contains(eFreePos2));
//         BOOST_CHECK(m.contains(eFreeTime));
//         BOOST_CHECK(not m.contains(eFreeDir0));
//         BOOST_CHECK(not m.contains(eFreeDir1));
//         BOOST_CHECK(not m.contains(eFreeDir2));
//         BOOST_CHECK(not m.contains(eFreeQOverP));
//       },
//       meas);

//   // reassign w/ all parameters
//   auto [parN, covN] = generateFreeParametersCovariance(rng);
//   meas =
//       makeMeasurement(source, parN, covN, eFreePos0, eFreePos1, eFreePos2,
//                       eFreeTime, eFreeDir0, eFreeDir1, eFreeDir2,
//                       eFreeQOverP);
//   std::visit(
//       [](const auto& m) {
//         BOOST_CHECK_EQUAL(m.size(), eFreeSize);
//         BOOST_CHECK(m.contains(eFreePos0));
//         BOOST_CHECK(m.contains(eFreePos1));
//         BOOST_CHECK(m.contains(eFreePos2));
//         BOOST_CHECK(m.contains(eFreeTime));
//         BOOST_CHECK(m.contains(eFreeDir0));
//         BOOST_CHECK(m.contains(eFreeDir1));
//         BOOST_CHECK(m.contains(eFreeDir2));
//         BOOST_CHECK(m.contains(eFreeQOverP));
//       },
//       meas);
// }

BOOST_AUTO_TEST_SUITE_END()
