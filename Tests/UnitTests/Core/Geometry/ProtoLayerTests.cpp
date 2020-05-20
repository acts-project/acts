// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/ProtoLayer.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"

#include <cmath>

namespace Acts {

namespace Test {
namespace Layers {
BOOST_AUTO_TEST_SUITE(Geometry)

BOOST_AUTO_TEST_CASE(ProtoLayerTests) {
  GeometryContext tgContext = GeometryContext();

  // Create a proto layer with 4 surfaces on the x/y grid
  auto recBounds = std::make_shared<RectangleBounds>(3., 6.);

  // Planar definitions to help construct the boundary surfaces
  static const Transform3D planeYZ =
      AngleAxis3D(0.5 * M_PI, Vector3D::UnitY()) *
      AngleAxis3D(0.5 * M_PI, Vector3D::UnitZ()) * Transform3D::Identity();
  static const Transform3D planeZX =
      AngleAxis3D(-0.5 * M_PI, Vector3D::UnitX()) *
      AngleAxis3D(-0.5 * M_PI, Vector3D::UnitZ()) * Transform3D::Identity();

  auto createProtoLayer = [&](const Transform3D& trf,
                              bool shared = false) -> ProtoLayer {
    auto atNegX = Surface::makeShared<PlaneSurface>(
        std::make_shared<Transform3D>(
            trf * Translation3D(Vector3D(-3., 0., 0.)) * planeYZ),
        recBounds);

    auto atPosX = Surface::makeShared<PlaneSurface>(
        std::make_shared<Transform3D>(
            trf * Translation3D(Vector3D(3., 0., 0.)) * planeYZ),
        recBounds);

    auto atNegY = Surface::makeShared<PlaneSurface>(
        std::make_shared<Transform3D>(
            trf * Translation3D(Vector3D(0., -3, 0.)) * planeZX),
        recBounds);

    auto atPosY = Surface::makeShared<PlaneSurface>(
        std::make_shared<Transform3D>(
            trf * Translation3D(Vector3D(0., 3., 0.)) * planeZX),
        recBounds);

    if (not shared) {
      std::vector<const Surface*> surfaces = {atNegX.get(), atNegY.get(),
                                              atPosX.get(), atPosY.get()};

      return ProtoLayer(tgContext, surfaces);
    }
    std::vector<std::shared_ptr<const Surface>> sharedSurfaces = {
        atNegX, atNegY, atPosX, atPosY};
    return ProtoLayer(tgContext, sharedSurfaces);
  };

  // Test 0 - check constructor with surfaces and shared surfaces
  auto pLayerSf = createProtoLayer(Transform3D::Identity());
  auto pLayerSfShared = createProtoLayer(Transform3D::Identity());

  BOOST_CHECK(pLayerSf.extent.ranges == pLayerSfShared.extent.ranges);
  BOOST_CHECK(pLayerSf.envelope == pLayerSfShared.envelope);

  // CHECK That you have 4 surfaces
  BOOST_CHECK(pLayerSf.surfaces().size() == 4);
  // Add one surface
  auto rB = std::make_shared<RectangleBounds>(30., 60.);
  auto pSurface = Surface::makeShared<PlaneSurface>(
      std::make_shared<Transform3D>(Transform3D::Identity()), rB);
  pLayerSf.add(tgContext, *pSurface.get());
  // CHECK That if you now have 5 surfaces
  BOOST_CHECK(pLayerSf.surfaces().size() == 5);

  // That should invalidate the ranges
  BOOST_CHECK(pLayerSf.extent.ranges != pLayerSfShared.extent.ranges);

  // Test 1 - identity transform
  auto protoLayer = createProtoLayer(Transform3D::Identity());

  CHECK_CLOSE_ABS(protoLayer.range(binX), 12., 1e-8);
  CHECK_CLOSE_ABS(protoLayer.medium(binX), 0., 1e-8);
  CHECK_CLOSE_ABS(protoLayer.min(binX), -6., 1e-8);
  CHECK_CLOSE_ABS(protoLayer.max(binX), 6., 1e-8);
  CHECK_CLOSE_ABS(protoLayer.range(binY), 6., 1e-8);
  CHECK_CLOSE_ABS(protoLayer.medium(binY), 0., 1e-8);
  CHECK_CLOSE_ABS(protoLayer.min(binY), -3., 1e-8);
  CHECK_CLOSE_ABS(protoLayer.max(binY), 3., 1e-8);
  CHECK_CLOSE_ABS(protoLayer.range(binZ), 12., 1e-8);
  CHECK_CLOSE_ABS(protoLayer.medium(binZ), 0., 1e-8);
  CHECK_CLOSE_ABS(protoLayer.min(binZ), -6., 1e-8);
  CHECK_CLOSE_ABS(protoLayer.max(binZ), 6., 1e-8);
  CHECK_CLOSE_ABS(protoLayer.max(binR), std::sqrt(3 * 3 + 6 * 6), 1e-8);
  CHECK_CLOSE_ABS(protoLayer.min(binR), 3., 1e-8);

  // Test 1a

  // Test 2 - rotate around Z-Axis, should leave R, Z untouched,
  // only preserves medium values
  auto protoLayerRot = createProtoLayer(AngleAxis3D(-0.345, Vector3D::UnitZ()) *
                                        Transform3D::Identity());

  BOOST_CHECK_NE(protoLayer.min(binX), -6.);
  CHECK_CLOSE_ABS(protoLayerRot.medium(binX), 0., 1e-8);
  CHECK_CLOSE_ABS(protoLayerRot.medium(binY), 0., 1e-8);
  CHECK_CLOSE_ABS(protoLayerRot.range(binZ), 12., 1e-8);
  CHECK_CLOSE_ABS(protoLayerRot.medium(binZ), 0., 1e-8);
  CHECK_CLOSE_ABS(protoLayerRot.min(binZ), -6., 1e-8);
  CHECK_CLOSE_ABS(protoLayerRot.max(binZ), 6., 1e-8);
  CHECK_CLOSE_ABS(protoLayerRot.min(binR), 3., 1e-8);
  CHECK_CLOSE_ABS(protoLayerRot.max(binR), std::sqrt(3 * 3 + 6 * 6), 1e-8);
}

BOOST_AUTO_TEST_SUITE_END()
}  // namespace Layers
}  // namespace Test

}  // namespace Acts
