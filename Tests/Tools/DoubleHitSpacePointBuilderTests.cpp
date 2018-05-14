// This file is part of the ACTS project.
//
// Copyright (C) 2018 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#define BOOST_TEST_MODULE Strip Space Point Builder Tests
#include <boost/test/included/unit_test.hpp>

#include <boost/test/data/test_case.hpp>
#include "ACTS/Tools/DoubleHitSpacePointBuilder.hpp"
#include "ACTS/Utilities/Definitions.hpp"

#include "ACTS/Digitization/PlanarModuleCluster.hpp"
#include "ACTS/Surfaces/PlaneSurface.hpp"
#include "ACTS/Surfaces/RectangleBounds.hpp"
#include "DetectorElementStub.hpp"

namespace bdata = boost::unit_test::data;
namespace tt    = boost::test_tools;

namespace Acts {

namespace Test {

  /// Unit test for testing the main functions of DoubleHitSpacePointBuilder
  /// 1) A pair of hits gets added and resolved.
  /// 2) A pair of hits gets added and rejected.
  BOOST_DATA_TEST_CASE(DoubleHitsSpacePointBuilder_basic,
                       bdata::xrange(1),
                       index)
  {
    (void)index;

    std::cout << "Create first hit" << std::endl;

    // Build Bounds
    std::shared_ptr<const RectangleBounds> recBounds(
        new RectangleBounds(35. * Acts::units::_um, 25. * units::_mm));

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

    std::cout << "Create second hit" << std::endl;

    // Build second PlanarModuleCluster
    const Identifier id2(1);

    double           rotation2 = -0.026;
    RotationMatrix3D rotationNeg;
    Vector3D         xNeg(cos(rotation2), sin(rotation2), 0.);
    Vector3D         yNeg(-sin(rotation2), cos(rotation2), 0.);
    Vector3D         zNeg(0., 0., 1.);
    rotationNeg.col(0) = xNeg;
    rotationNeg.col(1) = yNeg;
    rotationNeg.col(2) = zNeg;
    Transform3D t3d2   = getTransformFromRotTransl(
        rotationNeg, Vector3D(0., 0., 10.005 * units::_m));

    DetectorElementStub detElem2(
        id2,
        std::make_shared<const Transform3D>(t3d2),
        std::make_shared<const DigitizationModule>(digMod));
    PlaneSurface pSur2(recBounds, detElem2, detElem2.identify());

    PlanarModuleCluster* pmc2
        = new PlanarModuleCluster(pSur2,
                                  Identifier(1),
                                  cov,
                                  local[0],
                                  local[1],
                                  {DigitizationCell(0, 0, 1.)});

    std::cout << "Store both hits" << std::endl;

    std::vector<DoubleHitSpacePoint>           data;
    std::shared_ptr<DoubleHitSpacePointConfig> dhsp_cfg
        = std::make_shared<DoubleHitSpacePointConfig>(
            DoubleHitSpacePointConfig());
    ;

    // Combine two PlanarModuleClusters
    SPB::addHits<DoubleHitSpacePoint, DoubleHitSpacePointConfig>(
        data, {pmc}, {pmc2}, dhsp_cfg);

    BOOST_TEST(data.size() == 1, "Failed to add element");
    BOOST_TEST(*(data[0].hitModule1) == *pmc, "Failed to set hit");
    BOOST_TEST(*(data[0].hitModule2) == *pmc2, "Failed to set hit");

    std::cout << "Calculate space point" << std::endl;

    SPB::calculateSpacePoints<DoubleHitSpacePoint, DoubleHitSpacePointConfig>(
        data, dhsp_cfg);

    BOOST_TEST(data[0].spacePoint != Vector3D::Zero(3),
               "Failed to calculate space point");

    std::cout << "Create third hit" << std::endl;

    // Build third PlanarModuleCluster
    const Identifier id3(2);
    Transform3D      t3d3 = getTransformFromRotTransl(
        rotationNeg, Vector3D(0., 10. * units::_m, 10.005 * units::_m));

    DetectorElementStub detElem3(
        id3,
        std::make_shared<const Transform3D>(t3d3),
        std::make_shared<const DigitizationModule>(digMod));
    PlaneSurface pSur3(recBounds, detElem3, detElem3.identify());

    PlanarModuleCluster* pmc3
        = new PlanarModuleCluster(pSur3,
                                  Identifier(2),
                                  cov,
                                  local[0],
                                  local[1],
                                  {DigitizationCell(0, 0, 1.)});

    std::cout << "Try to store hits" << std::endl;

    // Combine points
    SPB::addHits<DoubleHitSpacePoint, DoubleHitSpacePointConfig>(
        data, {pmc}, {pmc3}, dhsp_cfg);

    // Test for rejecting unconnected hits
    BOOST_TEST(data.size() == 1, "Failed to reject potential combination");
  }
}  // end of namespace Test

}  // end of namespace Acts
