// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#define BOOST_TEST_MODULE Pixel Space Point Builder Tests
#include <boost/test/included/unit_test.hpp>

#include <boost/test/data/test_case.hpp>
#include "Acts/EventData/Measurement.hpp"
#include "Acts/Plugins/Digitzation/SingleHitSpacePointBuilder.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Utilities/Definitions.hpp"

namespace bdata = boost::unit_test::data;
namespace tt    = boost::test_tools;

namespace Acts {

namespace Test {

  /// Unit test for testing the main functions of OneHitSpacePointBuilder
  /// 1) A resolved dummy hit gets created and added.
  /// 2) A hit gets added and resolved.
  BOOST_DATA_TEST_CASE(SingleHitSpacePointBuilder_basic,
                       bdata::xrange(1),
                       index)
  {
    (void)index;

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
    PlaneSurface      ps(std::make_shared<const Transform3D>(t3d));
    ActsSymMatrixD<2> cov;
    cov << 0., 0., 0., 0.;

    // Build Cluster
    auto cc = new Measurement<int, eLOC_0, eLOC_1>(ps, 0, cov, 0.1, -0.1);

    std::cout << "Hit created" << std::endl;

    std::vector<SingleHitSpacePoint<Measurement<int, eLOC_0, eLOC_1>>> data;
    SpacePointBuilder<SingleHitSpacePoint<Measurement<int, eLOC_0, eLOC_1>>>
        shsp;

    std::cout << "Hit added to storage" << std::endl;

    shsp.calculateSpacePoints({cc}, data);
    BOOST_TEST(data[0].spacePoint != Vector3D::Zero(3),
               "Failed to calculate space point");

    std::cout << "Space point calculated" << std::endl;
  }
}  // end of namespace Test

}  // end of namespace Acts
