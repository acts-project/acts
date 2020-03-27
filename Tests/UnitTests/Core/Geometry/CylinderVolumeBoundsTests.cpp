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
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Tests/CommonHelpers/ObjTestWriter.hpp"
#include "Acts/Utilities/BoundingBox.hpp"
#include "Acts/Utilities/Definitions.hpp"

namespace bdata = boost::unit_test::data;
namespace tt = boost::test_tools;

namespace Acts {

namespace Test {

BOOST_AUTO_TEST_SUITE(Geometry)

BOOST_AUTO_TEST_CASE(CylinderVolumeBoundsConstruction) {
  double rmin{10.}, rmax{20.}, halfz{30.}, halfphi{M_PI / 4}, avgphi{0.};

  // Test different construciton modes: solid
  CylinderVolumeBounds solidCylinder(0., rmax, halfz);
  BOOST_CHECK_EQUAL(solidCylinder.decomposeToSurfaces().size(), 3);

  // Test different construciton modes: sectoral solid
  CylinderVolumeBounds solidCylinderSector(0., rmax, halfz, halfphi);
  BOOST_CHECK_EQUAL(solidCylinderSector.decomposeToSurfaces().size(), 5);

  // Test different construciton modes: tube
  CylinderVolumeBounds tubeCylinder(rmin, rmax, halfz);
  BOOST_CHECK_EQUAL(tubeCylinder.decomposeToSurfaces().size(), 4);

  // Test different construciton modes: sectoral tube
  CylinderVolumeBounds tubeCylinderSector(rmin, rmax, halfz, halfphi);
  BOOST_CHECK_EQUAL(tubeCylinderSector.decomposeToSurfaces().size(), 6);

  CylinderVolumeBounds original(rmin, rmax, halfz, halfphi, avgphi);

  // Test construction from CylinderBounds and thickness
  double rmed = 0.5 * (rmin + rmax);
  double rthickness = (rmax - rmin);
  CylinderBounds cBounds(rmed, halfz, halfphi, avgphi);
  CylinderVolumeBounds fromCylinder(cBounds, rthickness);
  BOOST_CHECK_EQUAL(original, fromCylinder);

  // Test construction from RadialBounds and thickness
  RadialBounds rBounds(rmin, rmax, halfphi, avgphi);
  CylinderVolumeBounds fromDisc(rBounds, 2 * halfz);
  BOOST_CHECK_EQUAL(original, fromDisc);

  // Test the copy construction
  CylinderVolumeBounds copied(original);
  BOOST_CHECK_EQUAL(original, copied);

  // Test the assignment
  CylinderVolumeBounds assigned = original;
  BOOST_CHECK_EQUAL(original, assigned);
}

BOOST_AUTO_TEST_CASE(CylinderVolumeBoundsRecreation) {
  double rmin{10.}, rmax{20.}, halfz{30.}, halfphi{M_PI / 4}, avgphi{0.};

  CylinderVolumeBounds original(rmin, rmax, halfz, halfphi, avgphi);
  std::array<double, CylinderVolumeBounds::eSize> values;
  std::vector<double> valvector = original.values();
  std::copy_n(valvector.begin(), CylinderVolumeBounds::eSize, values.begin());
  CylinderVolumeBounds recreated(values);
  BOOST_CHECK_EQUAL(original, recreated);
}

BOOST_AUTO_TEST_CASE(CylinderVolumeBoundsExceptions) {
  double rmin{10.}, rmax{20.}, halfz{30.}, halfphi{M_PI / 4}, avgphi{0.};

  // Negative inner radius
  BOOST_CHECK_THROW(CylinderVolumeBounds(-rmin, rmax, halfz, halfphi, avgphi),
                    std::logic_error);

  // Negative outer radius
  BOOST_CHECK_THROW(CylinderVolumeBounds(rmin, -rmax, halfz, halfphi, avgphi),
                    std::logic_error);

  // Swapped radii
  BOOST_CHECK_THROW(CylinderVolumeBounds(rmax, rmin, halfz, halfphi, avgphi),
                    std::logic_error);

  // Zero half length
  BOOST_CHECK_THROW(CylinderVolumeBounds(rmax, rmin, 0., halfphi, avgphi),
                    std::logic_error);

  // Negative half length
  BOOST_CHECK_THROW(CylinderVolumeBounds(rmax, rmin, -halfz, halfphi, avgphi),
                    std::logic_error);

  // Out of bounds half phi
  BOOST_CHECK_THROW(CylinderVolumeBounds(rmax, rmin, halfz, -4., avgphi),
                    std::logic_error);

  // Wrong positioning phi
  BOOST_CHECK_THROW(CylinderVolumeBounds(rmax, rmin, halfz, halfphi, 4.),
                    std::logic_error);

  // Test construction from CylinderBounds and thickness
  double rmed = 0.5 * (rmin + rmax);
  CylinderBounds cBounds(rmed, halfz, halfphi, avgphi);
  RadialBounds rBounds(rmin, rmax, halfphi, avgphi);

  // Negative thickness
  BOOST_CHECK_THROW(CylinderVolumeBounds(cBounds, -1.), std::logic_error);

  // Wrong thickness
  BOOST_CHECK_THROW(CylinderVolumeBounds(cBounds, 1000.), std::logic_error);

  // Test construction from RadialBounds and thickness
  BOOST_CHECK_THROW(CylinderVolumeBounds(rBounds, -1), std::logic_error);
}

BOOST_AUTO_TEST_CASE(CylinderVolumeBoundsAccess) {
  double rmin{10.}, rmax{20.}, halfz{30.}, halfphi{M_PI / 4}, avgphi{0.};
  CylinderVolumeBounds cvBounds(rmin, rmax, halfz, halfphi, avgphi);

  // Test the accessors
  BOOST_CHECK_EQUAL(cvBounds.get(CylinderVolumeBounds::eMinR), rmin);
  BOOST_CHECK_EQUAL(cvBounds.get(CylinderVolumeBounds::eMaxR), rmax);
  BOOST_CHECK_EQUAL(cvBounds.get(CylinderVolumeBounds::eHalfLengthZ), halfz);
  BOOST_CHECK_EQUAL(cvBounds.get(CylinderVolumeBounds::eHalfPhiSector),
                    halfphi);
  BOOST_CHECK_EQUAL(cvBounds.get(CylinderVolumeBounds::eAveragePhi), avgphi);
}

/// Unit test for testing the decomposeToSurfaces() function
BOOST_DATA_TEST_CASE(CylinderVolumeBoundsDecomposeToSurfaces,
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
  double rmin = 1.;
  double rmax = 2.;
  double halfz = 3.;
  CylinderVolumeBounds cylBounds(rmin, rmax, halfz);
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
                  cylBounds.get(CylinderVolumeBounds::eHalfLengthZ), 1e-12);
  CHECK_CLOSE_REL((pos - boundarySurfaces.at(1)->center(tgContext)).norm(),
                  cylBounds.get(CylinderVolumeBounds::eHalfLengthZ), 1e-12);
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
  CHECK_CLOSE_ABS(
      negDiscPosZ + cylBounds.get(CylinderVolumeBounds::eHalfLengthZ),
      centerPosZ, 1e-12);
  CHECK_CLOSE_ABS(
      posDiscPosZ - cylBounds.get(CylinderVolumeBounds::eHalfLengthZ),
      centerPosZ, 1e-12);
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

BOOST_AUTO_TEST_CASE(CylinderVolumeBoundsBoundingBox) {
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

  CylinderVolumeBounds cvb(0., 5, 10);
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
  cvb = CylinderVolumeBounds(5, 8, 13, angle);
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
