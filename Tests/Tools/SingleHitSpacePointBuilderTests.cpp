// This file is part of the ACTS project.
//
// Copyright (C) 2018 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#define BOOST_TEST_MODULE Pixel Space Point Builder Tests
#include <boost/test/included/unit_test.hpp>

#include <boost/test/data/test_case.hpp>
#include "ACTS/Tools/SingleHitSpacePointBuilder.hpp"
#include "ACTS/Utilities/Definitions.hpp"

#include "ACTS/Digitization/PlanarModuleCluster.hpp"
#include "ACTS/Surfaces/PlaneSurface.hpp"
#include "ACTS/Surfaces/RectangleBounds.hpp"
#include "DetectorElementStub.hpp"

namespace bdata = boost::unit_test::data;
namespace tt    = boost::test_tools;

namespace Acts {

namespace Test {

  /// Unit test for testing the main functions of OneHitSpacePointBuilder
  /// 1) A resolved dummy hit gets created and added.
  /// 2) A hit gets added and resolved.
  BOOST_DATA_TEST_CASE(SingleHitSpacePointBuilder_basic, bdata::xrange(1), index)
  {
    (void)index;

    // Build bounds
    std::shared_ptr<const RectangleBounds> recBounds(
        new RectangleBounds(35. * units::_um, 25. * units::_mm));

    // Build binning and segmentation
    std::vector<float> boundariesX, boundariesY;
    boundariesX.push_back(-35. * units::_um);
    boundariesX.push_back(35. * units::_um);
    boundariesY.push_back(-25. * units::_mm);
    boundariesY.push_back(25. * units::_mm);

    BinningData binDataX(BinningOption::open, BinningValue::binX, boundariesX);
    std::shared_ptr<BinUtility> buX(new BinUtility(binDataX));
    BinningData binDataY(BinningOption::open, BinningValue::binY, boundariesY);
    std::shared_ptr<BinUtility> buY(new BinUtility(binDataY));
    (*buX) += (*buY);

    std::shared_ptr<const Segmentation> segmentation(
        new CartesianSegmentation(buX, recBounds));

    // Build translation
    const Identifier id(0);

    double           rotation = 0.026;
    RotationMatrix3D rotationPos;
    Vector3D         xPos(cos(rotation), sin(rotation), 0.);
    Vector3D         yPos(-sin(rotation), cos(rotation), 0.);
    Vector3D         zPos(0., 0., 1.);
    rotationPos.col(0) = xPos;
    rotationPos.col(1) = yPos;
    rotationPos.col(2) = zPos;
    Transform3D t3d    = getTransformFromRotTransl(
        rotationPos, Vector3D(0., 0., 10. * units::_m));

    // Build Digitization
    const DigitizationModule digMod(segmentation, 1., 1., 0.);
    DetectorElementStub      detElem(
        id,
        std::make_shared<const Transform3D>(t3d),
        std::make_shared<const DigitizationModule>(digMod));
    PlaneSurface      pSur(recBounds, detElem, detElem.identify());
    ActsSymMatrixD<2> cov;
    cov << 0., 0., 0., 0.;
    Vector2D local = {0.1, -0.1};

    // Build PlanarModuleCluster
    PlanarModuleCluster* pmc
        = new PlanarModuleCluster(pSur,
                                  Identifier(0),
                                  cov,
                                  local[0],
                                  local[1],
                                  {DigitizationCell(0, 0, 1.)});

    std::cout << "Hit created" << std::endl;

    std::vector<SingleHitSpacePoint>                              data;
    SPB::addHits<SingleHitSpacePoint>(data, {pmc});

    // Test for adding a SingleHitSpacePoint
    BOOST_TEST(data.size() == 1, "Failed to add element");
    BOOST_TEST(*(data[0].hitModule) == *pmc, "Wrong element added");

    std::cout << "Hit added to storage" << std::endl;

    SPB::calculateSpacePoints<SingleHitSpacePoint>(data);
    BOOST_TEST(data[0].spacePoint != Vector3D::Zero(3),
               "Failed to calculate space point");

    std::cout << "Space point calculated" << std::endl;
  }
}  // end of namespace Test

}  // end of namespace Acts
