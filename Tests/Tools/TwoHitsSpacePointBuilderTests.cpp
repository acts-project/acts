// This file is part of the ACTS project.
//
// Copyright (C) 2017-2018 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#define BOOST_TEST_MODULE Strip Space Point Builder Tests
#include <boost/test/included/unit_test.hpp>

#include <boost/test/data/test_case.hpp>
#include "ACTS/Tools/TwoHitsSpacePointBuilder.hpp"
#include "ACTS/Utilities/Definitions.hpp"

#include "ACTS/Digitization/PlanarModuleCluster.hpp"
#include "ACTS/Surfaces/PlaneSurface.hpp"
#include "ACTS/Surfaces/RectangleBounds.hpp"
#include "DetectorElementStub.hpp"

namespace bdata = boost::unit_test::data;
namespace tt    = boost::test_tools;

namespace Acts {

namespace Test {

  /// Unit test for testing the main functions of TwoHitsSpacePointBuilder
  BOOST_DATA_TEST_CASE(TwoHitsSpacePointBuilder_basic, bdata::xrange(1), index)
  {
    (void)index;
    TwoHitsSpacePointBuilder::CombinedHits combHits;

    std::shared_ptr<const RectangleBounds> recBounds(
        new RectangleBounds(35. * Acts::units::_um, 25. * units::_mm));

    // Build Segmentation
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

    const DigitizationModule digMod(segmentation, 1., 1., 0.);
    DetectorElementStub      detElem(
        id,
        std::make_shared<const Transform3D>(t3d),
        std::make_shared<const DigitizationModule>(digMod));
    PlaneSurface      pSur(recBounds, detElem, detElem.identify());
    ActsSymMatrixD<2> cov;
    cov << 0., 0., 0., 0.;
    Vector2D             local = {0.1, -0.1};
    PlanarModuleCluster* pmc
        = new PlanarModuleCluster(pSur,
                                  Identifier(0),
                                  cov,
                                  local[0],
                                  local[1],
                                  {DigitizationCell(0, 0, 1.)});

    // Test for setting a TwoHitsSpacePointBuilder::CombinedHits
    combHits.hitModule1 = pmc;
    BOOST_TEST(combHits.hitModule1 == pmc,
               "Failed to set element in combHits.hitModule1");

    combHits.hitModule2 = pmc;
    BOOST_TEST(combHits.hitModule2 == pmc,
               "Failed to set element in combHits.hitModule2");

    double diff   = 0.;
    combHits.diff = diff;
    BOOST_TEST(combHits.diff == diff, "Failed to set element in combHits.diff");

    Vector3D spacePoint = {1., 1., 1.};
    combHits.spacePoint = spacePoint;
    BOOST_TEST(combHits.spacePoint == spacePoint,
               "Failed to set element in combHits.spacePoint");

    TwoHitsSpacePointBuilder::Config sspb_cfg;
    TwoHitsSpacePointBuilder         sspb(sspb_cfg);
    sspb.addCombinedHit(combHits);

    // Test for adding a TwoHitsSpacePointBuilder::CombinedHits
    const std::vector<TwoHitsSpacePointBuilder::CombinedHits> vecCombHits
        = sspb.combinedHits();
    BOOST_TEST(vecCombHits.size() == 1,
               "Failed to add element to SpacePointBuilder");
    BOOST_TEST(vecCombHits[0].hitModule1 == combHits.hitModule1,
               "Wrong element added");
    BOOST_TEST(vecCombHits[0].hitModule2 == combHits.hitModule2,
               "Wrong element added");
    BOOST_TEST(vecCombHits[0].diff == combHits.diff, "Wrong element added");
    BOOST_TEST(vecCombHits[0].spacePoint == combHits.spacePoint,
               "Wrong element added");

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

    sspb.combineHits({*pmc}, {*pmc2});
    sspb.calculateSpacePoints();

    const std::vector<TwoHitsSpacePointBuilder::CombinedHits> vecCombHits2
        = sspb.combinedHits();
    // Test for creating a new TwoHitsSpacePointBuilder::CombinedHits element
    // with
    // PlanarModuleClusters
    BOOST_TEST(vecCombHits2.size() == 2,
               "Failed to add element to SpacePointBuilder");

    // Test for calculating space points
    BOOST_TEST(vecCombHits2.back().spacePoint != Vector3D::Zero(3),
               "Failed to calculate space point");

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

    sspb.combineHits({*pmc}, {*pmc3});

    // Test for rejecting unconnected hits
    BOOST_TEST(sspb.combinedHits().size() == 2,
               "Failed to reject potential combination");
  }
}  // end of namespace Test

}  // end of namespace Acts
