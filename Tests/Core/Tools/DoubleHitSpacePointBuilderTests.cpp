// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#define BOOST_TEST_MODULE Strip Space Point Builder Tests
#include <boost/test/included/unit_test.hpp>

#include <boost/test/data/test_case.hpp>
#include "Acts/EventData/Measurement.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Plugins/Digitization/DoubleHitSpacePointBuilder.hpp"
#include "Acts/Utilities/Definitions.hpp"

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

    PlaneSurface      pSur(std::make_shared<const Transform3D>(t3d));
    ActsSymMatrixD<2> cov;
    cov << 0., 0., 0., 0.;

    // Build Clusters
    auto* cc = new Measurement<int, eLOC_0, eLOC_1>(pSur, 0, cov, 0.1, -0.1);

    std::cout << "Create second hit" << std::endl;

    // Build second Cluster
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

    PlaneSurface pSur2(std::make_shared<const Transform3D>(t3d2));

    auto* cc2 = new Measurement<int, eLOC_0, eLOC_1>(pSur2, 1, cov, 0.1, -0.1);

    std::cout << "Store both hits" << std::endl;

    std::vector<DoubleHitSpacePoint<Measurement<int, eLOC_0, eLOC_1>>> resultSP;
    std::vector<std::pair<const Measurement<int, eLOC_0, eLOC_1>*,
                          const Measurement<int, eLOC_0, eLOC_1>*>>
        clusterPairs;
    SpacePointBuilder<DoubleHitSpacePoint<Measurement<int, eLOC_0, eLOC_1>>>::
        DoubleHitSpacePointConfig dhsp_cfg;

    // Combine two PlanarModuleClusters
    SpacePointBuilder<DoubleHitSpacePoint<Measurement<int, eLOC_0, eLOC_1>>>
        dhsp(dhsp_cfg);
    dhsp.makeClusterPairs({cc}, {cc2}, clusterPairs);

    BOOST_TEST(clusterPairs.size() == 1, "Failed to add element");
    BOOST_TEST(*(clusterPairs[0].first) == *cc, "Failed to set hit");
    BOOST_TEST(*(clusterPairs[0].second) == *cc2, "Failed to set hit");

    Vector3D stripTop1(0., 25. * units::_mm, 1. * units::_m);
    Vector3D stripBottom1(0., -25. * units::_mm, 1. * units::_m);
    auto     strip1 = std::make_pair(stripTop1, stripBottom1);

    Vector3D stripTop2(0.1 * units::_mm, 25. * units::_mm, 1. * units::_m);
    Vector3D stripBottom2(-0.1 * units::_mm, -25. * units::_mm, 1. * units::_m);
    auto     strip2 = std::make_pair(stripTop2, stripBottom2);

    std::cout << "Calculate space point" << std::endl;

    dhsp.calculateSpacePoints(clusterPairs, {strip1}, {strip2}, resultSP);

    BOOST_TEST(resultSP.size() == 1, "Failed to calculate space point");

    std::cout << "Create third hit" << std::endl;

    // Build third Cluster
    Transform3D t3d3 = getTransformFromRotTransl(
        rotationNeg, Vector3D(0., 10. * units::_m, 10.005 * units::_m));
    PlaneSurface pSur3(std::make_shared<const Transform3D>(t3d3));

    auto* cc3 = new Measurement<int, eLOC_0, eLOC_1>(pSur3, 2, cov, 0.1, -0.1);

    std::cout << "Try to store hits" << std::endl;

    // Combine points
    dhsp.makeClusterPairs({cc}, {cc3}, clusterPairs);

    // Test for rejecting unconnected hits
    BOOST_TEST(resultSP.size() == 1, "Failed to reject potential combination");
  }
}  // end of namespace Test

}  // end of namespace Acts
