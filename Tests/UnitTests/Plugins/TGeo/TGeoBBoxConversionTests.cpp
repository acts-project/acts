// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Plugins/TGeo/TGeoSurfaceConverter.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Visualization/GeometryView.hpp"
#include "Acts/Visualization/ObjVisualization.hpp"
#include "TGeoBBox.h"
#include "TGeoManager.h"
#include "TGeoMaterial.h"
#include "TGeoMatrix.h"
#include "TGeoMedium.h"
#include "TGeoVolume.h"
#include "TView.h"

namespace bdata = boost::unit_test::data;
namespace tt = boost::test_tools;

namespace Acts {

namespace Test {

GeometryContext tgContext = GeometryContext();

ViewConfig red({200, 0, 0});
ViewConfig green({0, 200, 0});
ViewConfig blue({0, 0, 200});

/// @brief Unit test to convert a Bbox into a Plane
///
/// This test also tests:
/// * the "(x/X)(y/Y)(z/Z)" orientations
/// * the scaling functionality
BOOST_AUTO_TEST_CASE(TGeoBBox_to_PlaneSurface) {
  ObjVisualization objVis;

  double x = 10.;
  double y = 30.;
  double z = 1.;

  new TGeoManager("box", "poza1");
  TGeoMaterial *mat = new TGeoMaterial("Al", 26.98, 13, 2.7);
  TGeoMedium *med = new TGeoMedium("MED", 1, mat);
  TGeoVolume *top = gGeoManager->MakeBox("TOP", med, 100, 100, 100);
  gGeoManager->SetTopVolume(top);
  TGeoVolume *vol = gGeoManager->MakeBox("BOX", med, x, y, z);
  vol->SetLineWidth(2);
  top->AddNode(vol, 1);
  gGeoManager->CloseGeometry();

  // Upper case ---------------------------------
  auto plane_XYZ = TGeoSurfaceConverter::toSurface(*vol->GetShape(),
                                                   *gGeoIdentity, "XY*", 1);
  BOOST_TEST(plane_XYZ != nullptr);
  BOOST_TEST(plane_XYZ->type() == Surface::Plane);

  auto bounds_XYZ =
      dynamic_cast<const RectangleBounds *>(&(plane_XYZ->bounds()));
  BOOST_TEST(bounds_XYZ != nullptr);
  double maxX = bounds_XYZ->get(RectangleBounds::eMaxX);
  double minX = bounds_XYZ->get(RectangleBounds::eMinX);
  double maxY = bounds_XYZ->get(RectangleBounds::eMaxY);
  double minY = bounds_XYZ->get(RectangleBounds::eMinY);
  CHECK_CLOSE_ABS(maxX - minX, 2 * x, s_epsilon);
  CHECK_CLOSE_ABS(maxY - minY, 2 * y, s_epsilon);

  // Check if the surface is the (negative) identity
  auto transform_XYZ = plane_XYZ->transform(tgContext);
  auto rotation_XYZ = transform_XYZ.rotation();
  BOOST_TEST(transform_XYZ.isApprox(Transform3D::Identity()));

  const Vector3D offset_XYZ{-5.5 * x, 0., 0.};
  GeometryView::drawSurface(objVis, *plane_XYZ, tgContext,
                            Transform3D(Translation3D{offset_XYZ}));
  const Vector3D center_XYZ = plane_XYZ->center(tgContext) + offset_XYZ;
  GeometryView::drawArrowForward(
      objVis, center_XYZ,
      center_XYZ + 0.6 * (maxX - minX) * rotation_XYZ.col(0), 4., 2.5, red);
  GeometryView::drawArrowForward(
      objVis, center_XYZ,
      center_XYZ + 0.6 * (maxY - minY) * rotation_XYZ.col(1), 4., 2.5, green);
  GeometryView::drawArrowForward(
      objVis, center_XYZ, center_XYZ + 2 * rotation_XYZ.col(2), 4., 2.5, blue);

  // Lower case ---------------------------------
  auto plane_xyz = TGeoSurfaceConverter::toSurface(*vol->GetShape(),
                                                   *gGeoIdentity, "xy*", 1);
  BOOST_TEST(plane_xyz != nullptr);
  BOOST_TEST(plane_xyz->type() == Surface::Plane);

  auto bounds_xyz =
      dynamic_cast<const RectangleBounds *>(&(plane_XYZ->bounds()));
  BOOST_TEST(bounds_xyz != nullptr);
  BOOST_TEST(bounds_xyz == bounds_XYZ);
  auto transform_xyz = plane_xyz->transform(tgContext);
  auto rotation_xyz = transform_xyz.rotation();
  BOOST_TEST(rotation_xyz.col(0).isApprox(-1 * rotation_XYZ.col(0)));
  BOOST_TEST(rotation_xyz.col(1).isApprox(-1 * rotation_XYZ.col(1)));
  BOOST_TEST(rotation_xyz.col(2).isApprox(rotation_XYZ.col(2)));

  const Vector3D offset_xyz{-2 * x, 0., 0.};
  GeometryView::drawSurface(objVis, *plane_xyz, tgContext,
                            Transform3D(Translation3D{offset_xyz}));
  const Vector3D center_xyz = plane_xyz->center(tgContext) + offset_xyz;
  GeometryView::drawArrowForward(
      objVis, center_xyz,
      center_xyz + 0.6 * (maxX - minX) * rotation_xyz.col(0), 4., 2.5, red);
  GeometryView::drawArrowForward(
      objVis, center_xyz,
      center_xyz + 0.6 * (maxY - minY) * rotation_xyz.col(1), 4., 2.5, green);
  GeometryView::drawArrowForward(
      objVis, center_xyz, center_xyz + 2 * rotation_xyz.col(2), 4., 2.5, blue);

  // Mixed case ---------------------------------
  auto plane_xYz = TGeoSurfaceConverter::toSurface(*vol->GetShape(),
                                                   *gGeoIdentity, "xY*", 1);
  BOOST_TEST(plane_xYz != nullptr);
  BOOST_TEST(plane_xYz->type() == Surface::Plane);

  auto bounds_xYz =
      dynamic_cast<const RectangleBounds *>(&(plane_xYz->bounds()));
  BOOST_TEST(bounds_xYz != nullptr);
  BOOST_TEST(bounds_xYz == bounds_xYz);
  auto transform_xYz = plane_xYz->transform(tgContext);
  auto rotation_xYz = transform_xYz.rotation();
  BOOST_TEST(rotation_xYz.col(0).isApprox(-1 * rotation_XYZ.col(0)));
  BOOST_TEST(rotation_xYz.col(1).isApprox(rotation_XYZ.col(1)));
  BOOST_TEST(rotation_xYz.col(2).isApprox(-1. * rotation_XYZ.col(2)));

  const Vector3D offset_xYz{2 * x, 0., 0.};
  GeometryView::drawSurface(
      objVis, *plane_xYz, tgContext,
      Translation3D{offset_xYz} * Transform3D::Identity());
  const Vector3D center_xYz = plane_xYz->center(tgContext) + offset_xYz;
  GeometryView::drawArrowForward(
      objVis, center_xYz,
      center_xYz + 0.6 * (maxX - minX) * rotation_xYz.col(0), 4., 2.5, red);
  GeometryView::drawArrowForward(
      objVis, center_xYz,
      center_xYz + 0.6 * (maxY - minY) * rotation_xYz.col(1), 4., 2.5, green);
  GeometryView::drawArrowForward(
      objVis, center_xYz, center_xYz + 2 * rotation_xYz.col(2), 4., 2.5, blue);

  // Swap case --------------------------------- (x/y) here
  auto plane_YXz = TGeoSurfaceConverter::toSurface(*vol->GetShape(),
                                                   *gGeoIdentity, "YX*", 1);
  BOOST_TEST(plane_YXz != nullptr);
  BOOST_TEST(plane_YXz->type() == Surface::Plane);
  auto bounds_YXz =
      dynamic_cast<const RectangleBounds *>(&(plane_YXz->bounds()));
  maxX = bounds_YXz->get(RectangleBounds::eMaxX);
  minX = bounds_YXz->get(RectangleBounds::eMinX);
  maxY = bounds_YXz->get(RectangleBounds::eMaxY);
  minY = bounds_YXz->get(RectangleBounds::eMinY);
  CHECK_CLOSE_ABS(maxX - minX, 2 * y, s_epsilon);
  CHECK_CLOSE_ABS(maxY - minY, 2 * x, s_epsilon);

  auto transform_YXz = plane_YXz->transform(tgContext);
  auto rotation_YXz = transform_YXz.rotation();
  BOOST_TEST(rotation_YXz.col(0).isApprox(rotation_XYZ.col(1)));
  BOOST_TEST(rotation_YXz.col(1).isApprox(rotation_XYZ.col(0)));
  BOOST_TEST(rotation_YXz.col(2).isApprox(-1. * rotation_XYZ.col(2)));

  const Vector3D offset_YXz{5.5 * x, 0., 0.};
  GeometryView::drawSurface(objVis, *plane_YXz, tgContext,
                            Transform3D(Translation3D{offset_YXz}));
  const Vector3D center_YXz = plane_YXz->center(tgContext) + offset_YXz;
  GeometryView::drawArrowForward(
      objVis, center_YXz,
      center_YXz + 0.6 * (maxX - minX) * rotation_YXz.col(0), 4., 2.5, red);
  GeometryView::drawArrowForward(
      objVis, center_YXz,
      center_YXz + 0.6 * (maxY - minY) * rotation_YXz.col(1), 4., 2.5, green);
  GeometryView::drawArrowForward(
      objVis, center_YXz, center_YXz + 2 * rotation_YXz.col(2), 4., 2.5, blue);

  // Scaling test ---------------------------------
  auto plane_XYZ10 = TGeoSurfaceConverter::toSurface(*vol->GetShape(),
                                                     *gGeoIdentity, "xY*", 10);
  BOOST_TEST(plane_XYZ10 != nullptr);

  auto bounds_XYZ10 =
      dynamic_cast<const RectangleBounds *>(&(plane_XYZ10->bounds()));
  double maxX10 = bounds_XYZ10->get(RectangleBounds::eMaxX);
  double minX10 = bounds_XYZ10->get(RectangleBounds::eMinX);
  double maxY10 = bounds_XYZ10->get(RectangleBounds::eMaxY);
  double minY10 = bounds_XYZ10->get(RectangleBounds::eMinY);
  CHECK_CLOSE_ABS(maxX10 - minX10, 20 * x, s_epsilon);
  CHECK_CLOSE_ABS(maxY10 - minY10, 20 * y, s_epsilon);

  objVis.write("TGeoConversion_TGeoBBox_PlaneSurface");
}

}  // namespace Test

}  // namespace Acts