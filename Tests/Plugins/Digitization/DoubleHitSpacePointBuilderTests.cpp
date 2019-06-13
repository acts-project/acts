// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#define BOOST_TEST_MODULE Strip Space Point Builder Tests
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <boost/test/data/test_case.hpp>
#include "Acts/Plugins/Digitization/DoubleHitSpacePointBuilder.hpp"
#include "Acts/Utilities/Definitions.hpp"

#include "Acts/Plugins/Digitization/PlanarModuleCluster.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Tests/CommonHelpers/DetectorElementStub.hpp"

namespace bdata = boost::unit_test::data;
namespace tt = boost::test_tools;
using namespace Acts::UnitLiterals;

namespace Acts {
namespace Test {

// Create a test context
GeometryContext tgContext = GeometryContext();

/// Unit test for testing the main functions of DoubleHitSpacePointBuilder
/// 1) A pair of hits gets added and resolved.
/// 2) A pair of hits gets added and rejected.
BOOST_DATA_TEST_CASE(DoubleHitsSpacePointBuilder_basic, bdata::xrange(1),
                     index) {
  (void)index;

  std::cout << "Create first hit" << std::endl;

  // Build Bounds
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

  double rotation = 0.026;
  RotationMatrix3D rotationPos;
  Vector3D xPos(cos(rotation), sin(rotation), 0.);
  Vector3D yPos(-sin(rotation), cos(rotation), 0.);
  Vector3D zPos(0., 0., 1.);
  rotationPos.col(0) = xPos;
  rotationPos.col(1) = yPos;
  rotationPos.col(2) = zPos;
  Transform3D t3d(Transform3D::Identity() * rotationPos);
  t3d.translation() = Vector3D(0., 0., 10_m);

  // Build Digitization
  const DigitizationModule digMod(segmentation, 1., 1., 0.);
  DetectorElementStub detElem(std::make_shared<const Transform3D>(t3d));
  auto pSur = Surface::makeShared<PlaneSurface>(recBounds, detElem);
  ActsSymMatrixD<2> cov;
  cov << 0., 0., 0., 0.;
  Vector2D local = {0.1, -0.1};

  // Build PlanarModuleCluster
  PlanarModuleCluster* pmc = new PlanarModuleCluster(
      pSur, {}, cov, local[0], local[1], {DigitizationCell(0, 0, 1.)}, &digMod);

  std::cout << "Create second hit" << std::endl;

  // Build second PlanarModuleCluster

  double rotation2 = -0.026;
  RotationMatrix3D rotationNeg;
  Vector3D xNeg(cos(rotation2), sin(rotation2), 0.);
  Vector3D yNeg(-sin(rotation2), cos(rotation2), 0.);
  Vector3D zNeg(0., 0., 1.);
  rotationNeg.col(0) = xNeg;
  rotationNeg.col(1) = yNeg;
  rotationNeg.col(2) = zNeg;
  Transform3D t3d2(Transform3D::Identity() * rotationNeg);
  t3d2.translation() = Vector3D(0., 0., 10.005_m);

  DetectorElementStub detElem2(std::make_shared<const Transform3D>(t3d2));

  auto pSur2 = Surface::makeShared<PlaneSurface>(recBounds, detElem2);

  PlanarModuleCluster* pmc2 =
      new PlanarModuleCluster(pSur2, {}, cov, local[0], local[1],
                              {DigitizationCell(0, 0, 1.)}, &digMod);

  std::cout << "Store both hits" << std::endl;

  std::vector<SpacePoint<PlanarModuleCluster>> resultSP;
  std::vector<std::pair<PlanarModuleCluster const*, PlanarModuleCluster const*>>
      clusterPairs;
  DoubleHitSpacePointConfig dhsp_cfg;

  // Combine two PlanarModuleClusters
  SpacePointBuilder<SpacePoint<PlanarModuleCluster>> dhsp(dhsp_cfg);
  dhsp.makeClusterPairs(tgContext, {pmc}, {pmc2}, clusterPairs);

  BOOST_CHECK_EQUAL(clusterPairs.size(), 1);
  BOOST_CHECK_EQUAL(*(clusterPairs[0].first), *pmc);
  BOOST_CHECK_EQUAL(*(clusterPairs[0].second), *pmc2);

  std::cout << "Calculate space point" << std::endl;

  dhsp.calculateSpacePoints(tgContext, clusterPairs, resultSP);

  BOOST_CHECK_EQUAL(resultSP.size(), 1);

  std::cout << "Create third hit" << std::endl;

  // Build third PlanarModuleCluster
  Transform3D t3d3(Transform3D::Identity() * rotationNeg);
  t3d3.translation() = Vector3D(0., 0., 10.005_m);

  DetectorElementStub detElem3(std::make_shared<const Transform3D>(t3d3));
  auto pSur3 = Surface::makeShared<PlaneSurface>(recBounds, detElem3);

  PlanarModuleCluster* pmc3 =
      new PlanarModuleCluster(pSur3, {}, cov, local[0], local[1],
                              {DigitizationCell(0, 0, 1.)}, &digMod);

  std::cout << "Try to store hits" << std::endl;

  // Combine points
  dhsp.makeClusterPairs(tgContext, {pmc}, {pmc3}, clusterPairs);

  // Test for rejecting unconnected hits
  BOOST_CHECK_EQUAL(resultSP.size(), 1);
}

}  // end of namespace Test
}  // end of namespace Acts
