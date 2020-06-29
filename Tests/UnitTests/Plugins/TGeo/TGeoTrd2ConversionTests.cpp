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
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/TrapezoidBounds.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Visualization/GeometryView.hpp"
#include "Acts/Visualization/ObjVisualization.hpp"
#include "TGeoManager.h"
#include "TGeoMaterial.h"
#include "TGeoMatrix.h"
#include "TGeoMedium.h"
#include "TGeoTrd1.h"
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

/// @brief Unit test to convert a TGeoTrd2 into a Plane
///
/// * The TGeoTrd2 has x/z orientation
BOOST_AUTO_TEST_CASE(TGeoTrd2_xz_to_PlaneSurface) {
  ObjVisualization objVis;

  double hxmin = 10.;
  double hxmax = 30.;
  double t = 1.;
  double hy = 40.;

  new TGeoManager("trd1", "poza9");
  TGeoMaterial *mat = new TGeoMaterial("Al", 26.98, 13, 2.7);
  TGeoMedium *med = new TGeoMedium("MED", 1, mat);
  TGeoVolume *top = gGeoManager->MakeBox("TOP", med, 100, 100, 100);
  gGeoManager->SetTopVolume(top);
  TGeoVolume *vol = gGeoManager->MakeTrd2("Trd2", med, hxmin, hxmax, t, t, hy);
  gGeoManager->CloseGeometry();

  // Check the 4 possible ways
  std::vector<std::string> axesTypes = {"XZ*", "xZ*", "xz*", "Xz*"};

  size_t itrd = 0;
  for (const auto &axes : axesTypes) {
    auto plane = TGeoSurfaceConverter::toSurface(*vol->GetShape(),
                                                 *gGeoIdentity, axes, 1);
    BOOST_TEST(plane != nullptr);
    BOOST_TEST(plane->type() == Surface::Plane);

    auto bounds = dynamic_cast<const TrapezoidBounds *>(&(plane->bounds()));
    BOOST_TEST(bounds != nullptr);
    double hXminY = bounds->get(TrapezoidBounds::eHalfLengthXnegY);
    double hXmaxY = bounds->get(TrapezoidBounds::eHalfLengthXposY);
    double hY = bounds->get(TrapezoidBounds::eHalfLengthY);

    CHECK_CLOSE_ABS(hxmin, std::min(hXminY, hXmaxY), s_epsilon);
    CHECK_CLOSE_ABS(hxmax, std::max(hXminY, hXmaxY), s_epsilon);
    CHECK_CLOSE_ABS(hy, hY, s_epsilon);

    // Check if the surface is the (negative) identity
    auto transform = plane->transform(tgContext);
    auto rotation = transform.rotation();
    const Vector3D offset{(-5.5 + (itrd++) * 2.5) * hxmax, 0., 0.};
    GeometryView::drawSurface(objVis, *plane, tgContext,
                              Translation3D{offset} * Transform3D::Identity());
    const Vector3D center = plane->center(tgContext) + offset;
    GeometryView::drawArrowForward(
        objVis, center, center + 1.2 * (hXminY + hXmaxY) * rotation.col(0), 4.,
        2.5, red);
    GeometryView::drawArrowForward(
        objVis, center, center + 1.2 * hY * rotation.col(1), 4., 2.5, green);
    GeometryView::drawArrowForward(objVis, center, center + 2 * rotation.col(2),
                                   4., 2.5, blue);
  }
  objVis.write("TGeoConversion_TGeoTrd2_xz_PlaneSurface");

  // Check exceptions for not allowed axis definition
  std::vector<std::string> notAllowed = {"XY*", "xy*", "Xy*", "xY*"};
  for (const auto &naxis : notAllowed) {
    BOOST_CHECK_THROW(TGeoSurfaceConverter::toSurface(*vol->GetShape(),
                                                      *gGeoIdentity, naxis, 1),
                      std::invalid_argument);
  }
}

/// @brief Unit test to convert a TGeoTrd2 into a Plane
///
/// * The TGeoTrd2 has y/z orientation
BOOST_AUTO_TEST_CASE(TGeoTrd2_yz_to_PlaneSurface) {
  ObjVisualization objVis;

  double hxmin = 10.;
  double hxmax = 30.;
  double t = 1.;
  double hy = 40.;

  new TGeoManager("trd1", "poza9");
  TGeoMaterial *mat = new TGeoMaterial("Al", 26.98, 13, 2.7);
  TGeoMedium *med = new TGeoMedium("MED", 1, mat);
  TGeoVolume *top = gGeoManager->MakeBox("TOP", med, 100, 100, 100);
  gGeoManager->SetTopVolume(top);
  TGeoVolume *vol = gGeoManager->MakeTrd2("Trd2", med, t, t, hxmin, hxmax, hy);
  gGeoManager->CloseGeometry();

  // Check the 4 possible ways
  std::vector<std::string> axesTypes = {"YZ*", "yZ*", "yz*", "Yz*"};

  size_t itrd = 0;
  for (const auto &axes : axesTypes) {
    auto plane = TGeoSurfaceConverter::toSurface(*vol->GetShape(),
                                                 *gGeoIdentity, axes, 1);
    BOOST_TEST(plane != nullptr);
    BOOST_TEST(plane->type() == Surface::Plane);

    auto bounds = dynamic_cast<const TrapezoidBounds *>(&(plane->bounds()));
    BOOST_TEST(bounds != nullptr);
    double hXminY = bounds->get(TrapezoidBounds::eHalfLengthXnegY);
    double hXmaxY = bounds->get(TrapezoidBounds::eHalfLengthXposY);
    double hY = bounds->get(TrapezoidBounds::eHalfLengthY);

    CHECK_CLOSE_ABS(hxmin, std::min(hXminY, hXmaxY), s_epsilon);
    CHECK_CLOSE_ABS(hxmax, std::max(hXminY, hXmaxY), s_epsilon);
    CHECK_CLOSE_ABS(hy, hY, s_epsilon);

    // Check if the surface is the (negative) identity
    auto transform = plane->transform(tgContext);
    auto rotation = transform.rotation();
    const Vector3D offset{(-5.5 + (itrd++) * 2.5) * hxmax, 0., 0.};
    GeometryView::drawSurface(objVis, *plane, tgContext,
                              Translation3D{offset} * Transform3D::Identity());
    const Vector3D center = plane->center(tgContext) + offset;
    GeometryView::drawArrowForward(
        objVis, center, center + 1.2 * (hXminY + hXmaxY) * rotation.col(0), 4.,
        2.5, red);
    GeometryView::drawArrowForward(
        objVis, center, center + 1.2 * hY * rotation.col(1), 4., 2.5, green);
    GeometryView::drawArrowForward(objVis, center, center + 2 * rotation.col(2),
                                   4., 2.5, blue);
  }
  objVis.write("TGeoConversion_TGeoTrd2_yz_PlaneSurface");

  // Check exceptions for not allowed axis definition
  std::vector<std::string> notAllowed = {"YX*", "yx*", "yX*", "Yx*"};
  for (const auto &naxis : notAllowed) {
    BOOST_CHECK_THROW(TGeoSurfaceConverter::toSurface(*vol->GetShape(),
                                                      *gGeoIdentity, naxis, 1),
                      std::invalid_argument);
  }
}

}  // namespace Test

}  // namespace Acts