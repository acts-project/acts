// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/Polyhedron.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Tests/CommonHelpers/ObjTestWriter.hpp"
#include "Acts/Utilities/BoundingBox.hpp"
#include "Acts/Utilities/Definitions.hpp"

namespace bdata = boost::unit_test::data;
namespace tt = boost::test_tools;

namespace Acts {

namespace Test {

BOOST_AUTO_TEST_SUITE(Volumes)

/// Unit test for testing the decomposeToSurfaces() function
BOOST_DATA_TEST_CASE(CylinderVolumeBounds_decomposeToSurfaces,
                     bdata::random(-M_PI, M_PI) ^ bdata::random(-M_PI, M_PI) ^
                         bdata::random(-M_PI, M_PI) ^ bdata::random(-10., 10.) ^
                         bdata::random(-10., 10.) ^ bdata::random(-10., 10.) ^
                         bdata::xrange(100),
                     alpha, beta, gamma, posX, posY, posZ, index) {
  (void)index;

  // Create a test context
  GeometryContext tgContext = GeometryContext();

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
  auto mutableTransformPtr = std::make_shared<Transform3D>(Translation3D(pos));
  (*mutableTransformPtr) *= rotZ;
  (*mutableTransformPtr) *= rotY;
  (*mutableTransformPtr) *= rotX;
  auto transformPtr =
      std::const_pointer_cast<const Transform3D>(mutableTransformPtr);
  // get the boundary surfaces
  std::vector<std::shared_ptr<const Acts::Surface>> boundarySurfaces =
      cylBounds.decomposeToSurfaces(transformPtr.get());
  // Test

  // check if difference is halfZ - sign and direction independent
  CHECK_CLOSE_REL((pos - boundarySurfaces.at(0)->center(tgContext)).norm(),
                  cylBounds.halflengthZ(), 1e-12);
  CHECK_CLOSE_REL((pos - boundarySurfaces.at(1)->center(tgContext)).norm(),
                  cylBounds.halflengthZ(), 1e-12);
  // transform to local
  double posDiscPosZ =
      (transformPtr->inverse() * boundarySurfaces.at(1)->center(tgContext)).z();
  double centerPosZ = (transformPtr->inverse() * pos).z();
  double negDiscPosZ =
      (transformPtr->inverse() * boundarySurfaces.at(0)->center(tgContext)).z();
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
  CHECK_CLOSE_REL(
      transformPtr->rotation().col(2).dot(
          boundarySurfaces.at(1)->normal(tgContext, Acts::Vector2D(0., 0.))),
      1., 1e-12);
  // negative disc durface should point in negative direction in the frame of
  // the volume
  CHECK_CLOSE_REL(
      transformPtr->rotation().col(2).dot(
          boundarySurfaces.at(0)->normal(tgContext, Acts::Vector2D(0., 0.))),
      -1., 1e-12);
  // test in r
  CHECK_CLOSE_REL(boundarySurfaces.at(3)->center(tgContext), pos, 1e-12);
  CHECK_CLOSE_REL(boundarySurfaces.at(2)->center(tgContext), pos, 1e-12);
}

BOOST_AUTO_TEST_CASE(bounding_box_creation) {
  GeometryContext tgContext = GeometryContext();
  std::vector<IdentifiedPolyderon> tPolyhedrons;

  auto combineAndDecompose = [&](const SurfacePtrVector& surfaces,
                                 const std::string& name) -> void {
    std::string writeBase = std::string("CylinderVolumeBounds_") + name;

    Polyhedron phCombined;
    size_t is = 0;
    for (const auto& sf : surfaces) {
      Polyhedron phComponent = sf->polyhedronRepresentation(tgContext, 72);
      phCombined.merge(phComponent);
      tPolyhedrons.push_back(
          {writeBase + std::string("_comp_") + std::to_string(is++), false,
           phComponent});
    }
    tPolyhedrons.push_back({writeBase, false, phCombined});
  };

  float tol = 1e-4;

  CylinderVolumeBounds cvb(5, 10);
  auto cvbSurfaces = cvb.decomposeToSurfaces();
  combineAndDecompose(cvbSurfaces, "Solid");
  auto bb = cvb.boundingBox();
  ObjTestWriter::writeObj("CylinderVolumeBounds_Solid_BB", bb);

  Transform3D rot;
  rot = AngleAxis3D(M_PI / 2., Vector3D::UnitX());

  BOOST_CHECK_EQUAL(bb.entity(), nullptr);
  BOOST_CHECK_EQUAL(bb.max(), Vector3D(5, 5, 10));
  BOOST_CHECK_EQUAL(bb.min(), Vector3D(-5, -5, -10));

  bb = cvb.boundingBox(&rot);
  BOOST_CHECK_EQUAL(bb.entity(), nullptr);
  CHECK_CLOSE_ABS(bb.max(), Vector3D(5, 10, 5), tol);
  CHECK_CLOSE_ABS(bb.min(), Vector3D(-5, -10, -5), tol);

  cvb = CylinderVolumeBounds(5, 8, 12);
  bb = cvb.boundingBox();
  BOOST_CHECK_EQUAL(bb.entity(), nullptr);
  BOOST_CHECK_EQUAL(bb.max(), Vector3D(8, 8, 12));
  BOOST_CHECK_EQUAL(bb.min(), Vector3D(-8, -8, -12));
  cvbSurfaces = cvb.decomposeToSurfaces();
  combineAndDecompose(cvbSurfaces, "Tube");
  ObjTestWriter::writeObj("CylinderVolumeBounds_Tube_BB", bb);

  double angle = M_PI / 8.;
  cvb = CylinderVolumeBounds(5, 8, angle, 13);
  bb = cvb.boundingBox();
  BOOST_CHECK_EQUAL(bb.entity(), nullptr);
  CHECK_CLOSE_ABS(bb.max(), Vector3D(8, 8 * std::sin(angle), 13), tol);
  CHECK_CLOSE_ABS(
      bb.min(), Vector3D(5 * std::cos(angle), -8 * std::sin(angle), -13), tol);
  cvbSurfaces = cvb.decomposeToSurfaces();
  combineAndDecompose(cvbSurfaces, "TubeSector");
  ObjTestWriter::writeObj("CylinderVolumeBounds_TubeSector_BB", bb);

  rot = AngleAxis3D(M_PI / 2., Vector3D::UnitZ());
  bb = cvb.boundingBox(&rot);
  BOOST_CHECK_EQUAL(bb.entity(), nullptr);
  CHECK_CLOSE_ABS(bb.max(), Vector3D(8 * std::sin(angle), 8, 13), tol);
  CHECK_CLOSE_ABS(
      bb.min(), Vector3D(-8 * std::sin(angle), 5 * std::cos(angle), -13), tol);

  rot = AngleAxis3D(M_PI / 2., Vector3D(-2, 4, 5).normalized());
  bb = cvb.boundingBox(&rot);
  BOOST_CHECK_EQUAL(bb.entity(), nullptr);
  CHECK_CLOSE_ABS(bb.max(), Vector3D(8.40007, 15.2828, 3.88911), tol);
  CHECK_CLOSE_ABS(bb.min(), Vector3D(-7.27834, -8.12028, -14.2182), tol);
  ObjTestWriter::writeObj(tPolyhedrons);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Test

}  // namespace Acts
