// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/DiamondVolumeBounds.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Visualization/ObjVisualization3D.hpp"

using namespace Acts;
using namespace Acts::UnitLiterals;

namespace ActsTests {

GeometryContext gctx = GeometryContext();

BOOST_AUTO_TEST_SUITE(GeometrySuite)
BOOST_AUTO_TEST_CASE(DiamondVolumeBoundsCreation) {
  // Create a DiamondVolumeBounds
  double halfX1 = 10.0;
  double halfX2 = 12.0;
  double halfX3 = 12.0;
  double halfY1 = 5.0;
  double halfY2 = 10.0;
  double halfZ = 2.0;

  DiamondVolumeBounds polygonBounds(halfX1, halfX2, halfX3, halfY1, halfY2,
                                    halfZ);

  // Test correct parameter return
  BOOST_CHECK_EQUAL(polygonBounds.get(DiamondVolumeBounds::eHalfLengthX1),
                    halfX1);
  BOOST_CHECK_EQUAL(polygonBounds.get(DiamondVolumeBounds::eHalfLengthX2),
                    halfX2);
  BOOST_CHECK_EQUAL(polygonBounds.get(DiamondVolumeBounds::eHalfLengthX3),
                    halfX3);
  BOOST_CHECK_EQUAL(polygonBounds.get(DiamondVolumeBounds::eLengthY1), halfY1);
  BOOST_CHECK_EQUAL(polygonBounds.get(DiamondVolumeBounds::eLengthY2), halfY2);
  BOOST_CHECK_EQUAL(polygonBounds.get(DiamondVolumeBounds::eHalfLengthZ),
                    halfZ);

  BOOST_CHECK_EQUAL(polygonBounds.type(), VolumeBounds::eDiamond);

  // try to create a volume with this bounds and visualize it

  Transform3 transform = Transform3(Translation3(100., 25., -250.));
  // add a rotation also
  transform.rotate(AngleAxis3(30._degree, Vector3::UnitZ()));

  auto trackingVolume = std::make_unique<TrackingVolume>(
      transform, std::make_shared<DiamondVolumeBounds>(polygonBounds),
      "TestDiamondVolume");

  Acts::ObjVisualization3D helper;
  trackingVolume->visualize(helper, gctx, {.visible = true}, {.visible = true},
                            {.visible = true});
  helper.write(trackingVolume->volumeName() + ".obj");
}
BOOST_AUTO_TEST_CASE(DiamondVolumeBoundsInside) {
  // Create a DiamondVolumeBounds
  double halfX1 = 10.0;
  double halfX2 = 12.0;
  double halfX3 = 5.0;
  double halfY1 = 5.0;
  double halfY2 = 10.0;
  double halfZ = 2.0;

  DiamondVolumeBounds polygonBounds(halfX1, halfX2, halfX3, halfY1, halfY2,
                                    halfZ);

  // Points inside
  BOOST_CHECK(polygonBounds.inside(Vector3::Zero(), 0.1));
  BOOST_CHECK(polygonBounds.inside(Vector3(5.0, 2.0, 1.0), 0.1));
  BOOST_CHECK(polygonBounds.inside(Vector3(-5.0, -2.0, -1.0), 0.1));
  BOOST_CHECK(polygonBounds.inside(Vector3(8.5, 4.0, 0.0),
                                   0.0));  // check side inclided faces

  // Points outside
  BOOST_CHECK(!polygonBounds.inside(Vector3(15.0, 0.0, 0.0), 0.1));
  BOOST_CHECK(!polygonBounds.inside(Vector3(-13., 15.0, 0.0), 0.1));
  BOOST_CHECK(!polygonBounds.inside(Vector3(0.0, 0.0, 5.0), 0.1));
  BOOST_CHECK(!polygonBounds.inside(Vector3(-12., 5.0, 0.0), 0.1));
  BOOST_CHECK(!polygonBounds.inside(Vector3(8.5, 7.0, 0.0),
                                    0.0));  // check side inclided faces
}

BOOST_AUTO_TEST_CASE(DiamondBoundarySurfaces) {
  // Create a DiamondVolumeBounds
  double halfX1 = 10.0;
  double halfX2 = 12.0;
  double halfX3 = 8.0;
  double halfY1 = 5.0;
  double halfY2 = 10.0;
  double halfZ = 2.0;

  DiamondVolumeBounds polygonBounds(halfX1, halfX2, halfX3, halfY1, halfY2,
                                    halfZ);

  Transform3 transform = Transform3::Identity();

  auto surfaces = polygonBounds.orientedSurfaces(transform);

  // There should be 8 boundary surfaces
  BOOST_CHECK_EQUAL(surfaces.size(), 8);

  // Check if the normal vector is correct
  for (auto& os : surfaces) {
    auto osCenter = os.surface->center(gctx);
    // cast into PlaneSurface
    const auto* pSurface = dynamic_cast<const PlaneSurface*>(os.surface.get());
    BOOST_CHECK(pSurface != nullptr);
    auto osNormal = pSurface->normal(gctx);
    Vector3 insideTvb = osCenter + os.direction * osNormal;
    Vector3 outsideTvb = osCenter - os.direction * osNormal;
    BOOST_CHECK(polygonBounds.inside(insideTvb));
    BOOST_CHECK(!polygonBounds.inside(outsideTvb));
  }

  // Test orientation of the boundary surfaces
  const Vector3 xaxis = Vector3::UnitX();
  const Vector3 yaxis = Vector3::UnitY();
  const Vector3 zaxis = Vector3::UnitZ();

  using enum DiamondVolumeBounds::Face;
  using enum DiamondVolumeBounds::BoundValues;

  auto pFaceXY = surfaces[toUnderlying(PositiveZFaceXY)]
                     .surface->transform(gctx)
                     .rotation();
  BOOST_CHECK(pFaceXY.col(0).isApprox(xaxis));
  BOOST_CHECK(pFaceXY.col(1).isApprox(yaxis));
  BOOST_CHECK(pFaceXY.col(2).isApprox(zaxis));

  auto nFaceXY = surfaces[toUnderlying(NegativeZFaceXY)]
                     .surface->transform(gctx)
                     .rotation();
  BOOST_CHECK(nFaceXY.col(0).isApprox(xaxis));
  BOOST_CHECK(nFaceXY.col(1).isApprox(yaxis));
  BOOST_CHECK(nFaceXY.col(2).isApprox(zaxis));

  auto pFaceXZ = surfaces[toUnderlying(PositiveYFaceZX)]
                     .surface->transform(gctx)
                     .rotation();
  BOOST_CHECK(pFaceXZ.col(0).isApprox(zaxis));
  BOOST_CHECK(pFaceXZ.col(1).isApprox(xaxis));
  BOOST_CHECK(pFaceXZ.col(2).isApprox(yaxis));

  auto nFaceXZ = surfaces[toUnderlying(NegativeYFaceZX)]
                     .surface->transform(gctx)
                     .rotation();
  BOOST_CHECK(nFaceXZ.col(0).isApprox(zaxis));
  BOOST_CHECK(nFaceXZ.col(1).isApprox(xaxis));
  BOOST_CHECK(nFaceXZ.col(2).isApprox(yaxis));

  // these are the surfaces in positive y attached to x2 and x3 (rotation with
  // beta angle)

  auto pFaceYZ23 = surfaces[toUnderlying(PositiveXFaceYZ23)]
                       .surface->transform(gctx)
                       .rotation();

  // use the vertices to check the expected rotation
  Vector3 vertA(polygonBounds.get(eHalfLengthX3), polygonBounds.get(eLengthY2),
                0.);
  Vector3 vertB(polygonBounds.get(eHalfLengthX2), 0., 0.);

  Vector3 vecAB = (vertA - vertB).normalized();
  BOOST_CHECK(pFaceYZ23.col(0).isApprox(vecAB));
  BOOST_CHECK(pFaceYZ23.col(1).isApprox(zaxis));
  BOOST_CHECK(pFaceYZ23.col(2).isApprox(vecAB.cross(zaxis)));

  // this is expected to be rotated by the beta angle (angle between x3 and x2
  // edges) along z axis
  auto nFaceYZ23 = surfaces[toUnderlying(NegativeXFaceYZ23)]
                       .surface->transform(gctx)
                       .rotation();
  vertA = Vector3(-polygonBounds.get(eHalfLengthX3),
                  polygonBounds.get(eLengthY2), 0.);
  vertB = Vector3(-polygonBounds.get(eHalfLengthX2), 0., 0.);
  vecAB = (vertA - vertB).normalized();
  BOOST_CHECK(nFaceYZ23.col(0).isApprox(vecAB));
  BOOST_CHECK(nFaceYZ23.col(1).isApprox(zaxis));
  BOOST_CHECK(nFaceYZ23.col(2).isApprox(vecAB.cross(zaxis)));

  // these are the face surfaces in negative y attached to x1 and x2 (rotation
  // with alpha angle)
  auto pFaceYZ12 = surfaces[toUnderlying(PositiveXFaceYZ12)]
                       .surface->transform(gctx)
                       .rotation();

  vertA = Vector3(polygonBounds.get(eHalfLengthX1),
                  -polygonBounds.get(eLengthY1), 0.);
  vertB = Vector3(polygonBounds.get(eHalfLengthX2), 0., 0.);
  vecAB = (vertA - vertB).normalized();
  BOOST_CHECK(pFaceYZ12.col(0).isApprox(-vecAB));
  BOOST_CHECK(pFaceYZ12.col(1).isApprox(zaxis));
  BOOST_CHECK(pFaceYZ12.col(2).isApprox(vecAB.cross(-zaxis)));

  auto nFaceYZ12 = surfaces[toUnderlying(NegativeXFaceYZ12)]
                       .surface->transform(gctx)
                       .rotation();
  vertA = Vector3(-polygonBounds.get(eHalfLengthX1),
                  -polygonBounds.get(eLengthY1), 0.);
  vertB = Vector3(-polygonBounds.get(eHalfLengthX2), 0., 0.);
  vecAB = (vertA - vertB).normalized();
  BOOST_CHECK(nFaceYZ12.col(0).isApprox(-vecAB));
  BOOST_CHECK(nFaceYZ12.col(1).isApprox(zaxis));
  BOOST_CHECK(nFaceYZ12.col(2).isApprox(vecAB.cross(-zaxis)));
}
BOOST_AUTO_TEST_SUITE_END()
}  // namespace ActsTests
