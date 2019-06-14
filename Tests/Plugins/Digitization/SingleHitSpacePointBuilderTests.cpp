// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#define BOOST_TEST_MODULE Pixel Space Point Builder Tests
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <boost/test/data/test_case.hpp>

#include "Acts/Plugins/Digitization/PlanarModuleCluster.hpp"
#include "Acts/Plugins/Digitization/SingleHitSpacePointBuilder.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Tests/CommonHelpers/DetectorElementStub.hpp"
#include "Acts/Utilities/Units.hpp"

namespace bdata = boost::unit_test::data;
namespace tt = boost::test_tools;

namespace Acts {
namespace Test {

using namespace Acts::UnitLiterals;

// Create a test context
GeometryContext tgContext = GeometryContext();

/// Unit test for testing the main functions of OneHitSpacePointBuilder
/// 1) A resolved dummy hit gets created and added.
/// 2) A hit gets added and resolved.
BOOST_DATA_TEST_CASE(SingleHitSpacePointBuilder_basic, bdata::xrange(1),
                     index) {
  (void)index;

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
      pSur, {}, cov, local[0], local[1], {DigitizationCell(0, 0, 1.)});

  std::cout << "Hit created" << std::endl;

  std::vector<SpacePoint<PlanarModuleCluster>> data;
  SpacePointBuilder<SpacePoint<PlanarModuleCluster>> shsp;

  std::cout << "Hit added to storage" << std::endl;

  shsp.calculateSpacePoints(tgContext, {pmc}, data);
  BOOST_CHECK_NE(data[0].spacePoint, SpacePointVector::Zero());

  std::cout << "Space point calculated" << std::endl;
}

}  // end of namespace Test
}  // end of namespace Acts
