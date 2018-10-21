// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// clang-format off
#define BOOST_TEST_MODULE Cylinder Volume Bounds Tests
#include <boost/test/included/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
// clang-format on

#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Volumes/CylinderVolumeBounds.hpp"

namespace bdata = boost::unit_test::data;
namespace tt    = boost::test_tools;

namespace Acts {

namespace Test {

  /// Unit test for testing the decomposeToSurfaces() function
  BOOST_DATA_TEST_CASE(CylinderVolumeBounds_decomposeToSurfaces,
                       bdata::random(-M_PI, M_PI) ^ bdata::random(-M_PI, M_PI)
                           ^ bdata::random(-M_PI, M_PI)
                           ^ bdata::random(-10., 10.)
                           ^ bdata::random(-10., 10.)
                           ^ bdata::random(-10., 10.)
                           ^ bdata::xrange(100),
                       alpha,
                       beta,
                       gamma,
                       posX,
                       posY,
                       posZ,
                       index)
  {
    (void)index;
    // position of volume
    const Vector3D pos(posX, posY, posZ);
    // rotation around x axis
    AngleAxis3D rotX(alpha, Vector3D(1., 0., 0.));
    // rotation around y axis
    AngleAxis3D rotY(beta, Vector3D(0., 1., 0.));
    // rotation around z axis
    AngleAxis3D rotZ(gamma, Vector3D(0., 0., 1.));

    // create the cylinder bounds
    CylinderVolumeBounds cylBounds(1., 2., 3.);
    // create the transformation matrix
    auto mutableTransformPtr
        = std::make_shared<Transform3D>(Translation3D(pos));
    (*mutableTransformPtr) *= rotZ;
    (*mutableTransformPtr) *= rotY;
    (*mutableTransformPtr) *= rotX;
    auto transformPtr
        = std::const_pointer_cast<const Transform3D>(mutableTransformPtr);
    // get the boundary surfaces
    std::vector<std::shared_ptr<const Acts::Surface>> boundarySurfaces
        = cylBounds.decomposeToSurfaces(transformPtr);
    // Test

    // check if difference is halfZ - sign and direction independent
    CHECK_CLOSE_REL((pos - boundarySurfaces.at(0)->center()).norm(),
                    cylBounds.halflengthZ(),
                    1e-12);
    CHECK_CLOSE_REL((pos - boundarySurfaces.at(1)->center()).norm(),
                    cylBounds.halflengthZ(),
                    1e-12);
    // transform to local
    double posDiscPosZ
        = (transformPtr->inverse() * boundarySurfaces.at(1)->center()).z();
    double centerPosZ = (transformPtr->inverse() * pos).z();
    double negDiscPosZ
        = (transformPtr->inverse() * boundarySurfaces.at(0)->center()).z();
    // check if center of disc boundaries lies in the middle in z
    BOOST_CHECK_LT(centerPosZ, posDiscPosZ);
    BOOST_CHECK_GT(centerPosZ, negDiscPosZ);
    // check positions of disc boundarysurfaces
    // checks for zero value. double precision value is not exact.
    CHECK_CLOSE_ABS(negDiscPosZ + cylBounds.halflengthZ(), centerPosZ, 1e-12);
    CHECK_CLOSE_ABS(posDiscPosZ - cylBounds.halflengthZ(), centerPosZ, 1e-12);
    // orientation of disc surfaces
    // positive disc durface should point in positive direction in the frame of
    // the volume
    CHECK_CLOSE_REL(transformPtr->rotation().col(2).dot(
                        boundarySurfaces.at(1)->normal(Acts::Vector2D(0., 0.))),
                    1.,
                    1e-12);
    // negative disc durface should point in negative direction in the frame of
    // the volume
    CHECK_CLOSE_REL(transformPtr->rotation().col(2).dot(
                        boundarySurfaces.at(0)->normal(Acts::Vector2D(0., 0.))),
                    -1.,
                    1e-12);
    // test in r
    CHECK_CLOSE_REL(boundarySurfaces.at(3)->center(), pos, 1e-12);
    CHECK_CLOSE_REL(boundarySurfaces.at(2)->center(), pos, 1e-12);
  }

}  // namespace Test

}  // namespace Acts
