// This file is part of the ACTS project.
//
// Copyright (C) 2025 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Geometry/ConvexPolygonVolumeBounds.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Visualization/ObjVisualization3D.hpp"

using namespace Acts;

namespace ActsTests {

GeometryContext gctx = GeometryContext();

BOOST_AUTO_TEST_SUITE(GeometrySuite)
BOOST_AUTO_TEST_CASE(ConvexPolygonVolumeBoundsCreation) {
  // Create a ConvexPolygonVolumeBounds
  double halfX1 = 10.0;
  double halfX2 = 12.0;
  double halfX3 = 12.0;
  double halfY1 = 5.0;
  double halfY2 = 10.0;
  double halfZ = 2.0;

  ConvexPolygonVolumeBounds polygonBounds(halfX1, halfX2, halfX3, halfY1,
                                          halfY2, halfZ);

  // Test correct parameter return
  BOOST_CHECK_EQUAL(polygonBounds.get(ConvexPolygonVolumeBounds::eHalfLengthX1),
                    halfX1);
  BOOST_CHECK_EQUAL(polygonBounds.get(ConvexPolygonVolumeBounds::eHalfLengthX2),
                    halfX2);
  BOOST_CHECK_EQUAL(polygonBounds.get(ConvexPolygonVolumeBounds::eHalfLengthX3),
                    halfX3);
  BOOST_CHECK_EQUAL(polygonBounds.get(ConvexPolygonVolumeBounds::eHalfLengthY1),
                    halfY1);
  BOOST_CHECK_EQUAL(polygonBounds.get(ConvexPolygonVolumeBounds::eHalfLengthY2),
                    halfY2);
  BOOST_CHECK_EQUAL(polygonBounds.get(ConvexPolygonVolumeBounds::eHalfLengthZ),
                    halfZ);

  BOOST_CHECK_EQUAL(polygonBounds.type(), VolumeBounds::eConvexPolygon);

  // try to create a volume with this bounds and visualize it

  Transform3 transform = Transform3(Translation3(100., 25., -250.));
  // add a rotation also
  transform.rotate(AngleAxis3(M_PI / 6, Vector3::UnitZ()));

  std::unique_ptr<TrackingVolume> trackingVolume =
      std::make_unique<TrackingVolume>(
          transform, std::make_shared<ConvexPolygonVolumeBounds>(polygonBounds),
          "TestConvexPolygonVolume");

  Acts::ObjVisualization3D helper;
  trackingVolume->visualize(helper, gctx, {.visible = true}, {.visible = true},
                            {.visible = true});
  helper.write(trackingVolume->volumeName() + ".obj");
}
BOOST_AUTO_TEST_CASE(ConvexPolygonVolumeBoundsInside) {
  // Create a ConvexPolygonVolumeBounds
  double halfX1 = 10.0;
  double halfX2 = 12.0;
  double halfX3 = 5.0;
  double halfY1 = 5.0;
  double halfY2 = 10.0;
  double halfZ = 2.0;

  ConvexPolygonVolumeBounds polygonBounds(halfX1, halfX2, halfX3, halfY1,
                                          halfY2, halfZ);

  // Points inside
  BOOST_CHECK(polygonBounds.inside(Vector3(0.0, 0.0, 0.0), 0.1));
  BOOST_CHECK(polygonBounds.inside(Vector3(5.0, 2.0, 1.0), 0.1));
  BOOST_CHECK(polygonBounds.inside(Vector3(-5.0, -2.0, -1.0), 0.1));

  // Points outside
  BOOST_CHECK(!polygonBounds.inside(Vector3(15.0, 0.0, 0.0), 0.1));
  BOOST_CHECK(!polygonBounds.inside(Vector3(-13., 15.0, 0.0), 0.1));
  BOOST_CHECK(!polygonBounds.inside(Vector3(0.0, 0.0, 5.0), 0.1));
}

BOOST_AUTO_TEST_CASE(ConvexPolygonBoundarySurfaces) {
  // Create a ConvexPolygonVolumeBounds
  double halfX1 = 10.0;
  double halfX2 = 12.0;
  double halfX3 = 12.0;
  double halfY1 = 5.0;
  double halfY2 = 10.0;
  double halfZ = 2.0;

  ConvexPolygonVolumeBounds polygonBounds(halfX1, halfX2, halfX3, halfY1,
                                          halfY2, halfZ);

  Transform3 transform = Transform3(Translation3(0., 0., 0.));

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
  Vector3 xaxis(1., 0., 0.);
  Vector3 yaxis(0., 1., 0.);
  Vector3 zaxis(0., 0., 1.);

  using enum ConvexPolygonVolumeBounds::Face;

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
}
BOOST_AUTO_TEST_SUITE_END()
}  // namespace ActsTests