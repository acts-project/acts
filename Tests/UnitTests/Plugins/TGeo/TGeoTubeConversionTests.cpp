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
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Visualization/GeometryVisualization.hpp"
#include "Acts/Visualization/ObjVisualization.hpp"
#include "TGeoManager.h"
#include "TGeoMaterial.h"
#include "TGeoMatrix.h"
#include "TGeoMedium.h"
#include "TGeoTube.h"
#include "TGeoVolume.h"
#include "TView.h"

namespace Acts {

namespace Test {

GeometryContext tgContext = GeometryContext();

std::vector<std::string> allowedAxes = {"XY*", "Xy*", "xy*", "xY*",
                                        "YX*", "yx*", "yX*", "Yx*"};

std::vector<std::string> notAllowedAxes = {"YZ*", "ZX*", "ZY*"};

/// @brief Unit test to convert a TGeoTube into a CylinderSurface
///
/// * The TGeoTrd1 can only have (x/X)(y/Y) orientation
BOOST_AUTO_TEST_CASE(TGeoTube_to_CylinderSurface) {
  ObjVisualization objVis;

  double rmin = 10.;
  double rmax = 11;
  double hz = 40.;
  double phimin = -45.;
  double phimax = 45.;

  new TGeoManager("trd1", "poza9");
  TGeoMaterial *mat = new TGeoMaterial("Al", 26.98, 13, 2.7);
  TGeoMedium *med = new TGeoMedium("MED", 1, mat);
  TGeoVolume *top = gGeoManager->MakeBox("TOP", med, 100, 100, 100);
  gGeoManager->SetTopVolume(top);
  TGeoVolume *vol = gGeoManager->MakeTube("Tube", med, rmin, rmax, hz);
  TGeoVolume *vols =
      gGeoManager->MakeTubs("Tube", med, rmin, rmax, hz, phimin, phimax);
  gGeoManager->CloseGeometry();

  size_t icyl = 0;
  for (const auto &axes : allowedAxes) {
    auto cylinder = TGeoSurfaceConverter::toSurface(*vol->GetShape(),
                                                    *gGeoIdentity, axes, 1);
    BOOST_TEST(cylinder != nullptr);
    BOOST_TEST(cylinder->type() == Surface::Cylinder);

    auto bounds = dynamic_cast<const CylinderBounds *>(&(cylinder->bounds()));
    BOOST_TEST(bounds != nullptr);
    double bR = bounds->get(CylinderBounds::eR);
    double bhZ = bounds->get(CylinderBounds::eHalfLengthZ);

    CHECK_CLOSE_ABS(bR, 10.5, s_epsilon);
    CHECK_CLOSE_ABS(bhZ, hz, s_epsilon);

    auto transform = cylinder->transform(tgContext);
    auto rotation = transform.rotation();

    // Check if the surface is the (negative) identity
    GeometryVisualization::drawSurface(objVis, *cylinder, tgContext,
                                       Transform3D::Identity(), 80, false,
                                       {20, 20, 120});
    const Vector3D center = cylinder->center(tgContext);
    GeometryVisualization::drawArrowForward(objVis, center,
                                            center + 1.2 * bR * rotation.col(0),
                                            0.2, 4., 2.5, 72, {200, 0, 0});
    GeometryVisualization::drawArrowForward(objVis, center,
                                            center + 1.2 * bR * rotation.col(1),
                                            0.2, 4., 2.5, 72, {0, 200, 0});
    GeometryVisualization::drawArrowForward(
        objVis, center, center + 1.2 * bhZ * rotation.col(2), 0.2, 4., 2.5, 72,
        {0, 0, 200});

    objVis.write("TGeoConversion_TGeoTube_CylinderSurface_" +
                 std::to_string(icyl));
    objVis.clear();

    if (icyl < 2) {
      auto cylinderSegment = TGeoSurfaceConverter::toSurface(
          *vols->GetShape(), *gGeoIdentity, axes, 1);
      BOOST_TEST(cylinderSegment != nullptr);
      BOOST_TEST(cylinderSegment->type() == Surface::Cylinder);

      auto boundsSegment =
          dynamic_cast<const CylinderBounds *>(&(cylinderSegment->bounds()));
      BOOST_TEST(boundsSegment != nullptr);
      bR = boundsSegment->get(CylinderBounds::eR);
      bhZ = boundsSegment->get(CylinderBounds::eHalfLengthZ);
      double hphi = boundsSegment->get(CylinderBounds::eHalfPhiSector);
      double mphi = boundsSegment->get(CylinderBounds::eAveragePhi);
      CHECK_CLOSE_ABS(bR, 10.5, s_epsilon);
      CHECK_CLOSE_ABS(bhZ, hz, s_epsilon);
      CHECK_CLOSE_ABS(hphi, 0.25 * M_PI, s_epsilon);
      CHECK_CLOSE_ABS(mphi, 0., s_epsilon);
      GeometryVisualization::drawSurface(objVis, *cylinderSegment, tgContext,
                                         Transform3D::Identity(), 80, false,
                                         {20, 20, 120});
      GeometryVisualization::drawArrowForward(
          objVis, center, center + 1.2 * bR * rotation.col(0), 0.2, 4., 2.5, 72,
          {200, 0, 0});
      GeometryVisualization::drawArrowForward(
          objVis, center, center + 1.2 * bR * rotation.col(1), 0.2, 4., 2.5, 72,
          {0, 200, 0});
      GeometryVisualization::drawArrowForward(
          objVis, center, center + 1.2 * bhZ * rotation.col(2), 0.2, 4., 2.5,
          72, {0, 0, 200});
      objVis.write("TGeoConversion_TGeoTube_CylinderSegmentSurface_" +
                   std::to_string(icyl));
      objVis.clear();
    } else {
      BOOST_CHECK_THROW(TGeoSurfaceConverter::toSurface(*vols->GetShape(),
                                                        *gGeoIdentity, axes, 1),
                        std::invalid_argument);
    }
    ++icyl;
  }

  // Check exceptions for not allowed axis definition
  for (const auto &naxes : notAllowedAxes) {
    BOOST_CHECK_THROW(TGeoSurfaceConverter::toSurface(*vol->GetShape(),
                                                      *gGeoIdentity, naxes, 1),
                      std::invalid_argument);
  }
}

/// @brief Unit test to convert a TGeoTube into a DiscSurface
///
/// * The TGeoTrd1 can only have (x/X)(y/Y) orientation
BOOST_AUTO_TEST_CASE(TGeoTube_to_DiscSurface) {
  ObjVisualization objVis;

  double rmin = 5.;
  double rmax = 25;
  double hz = 2.;
  double phimin = -45.;
  double phimax = 45.;

  new TGeoManager("trd1", "poza9");
  TGeoMaterial *mat = new TGeoMaterial("Al", 26.98, 13, 2.7);
  TGeoMedium *med = new TGeoMedium("MED", 1, mat);
  TGeoVolume *top = gGeoManager->MakeBox("TOP", med, 100, 100, 100);
  gGeoManager->SetTopVolume(top);
  TGeoVolume *vol = gGeoManager->MakeTube("Tube", med, rmin, rmax, hz);
  vol->SetLineWidth(2);
  TGeoVolume *vols =
      gGeoManager->MakeTubs("Tube", med, rmin, rmax, hz, phimin, phimax);
  gGeoManager->CloseGeometry();

  size_t idisc = 0;
  for (const auto &axes : allowedAxes) {
    auto disc = TGeoSurfaceConverter::toSurface(*vol->GetShape(), *gGeoIdentity,
                                                axes, 1);
    BOOST_TEST(disc != nullptr);
    BOOST_TEST(disc->type() == Surface::Disc);

    auto bounds = dynamic_cast<const RadialBounds *>(&(disc->bounds()));
    BOOST_TEST(bounds != nullptr);
    double bminr = bounds->get(RadialBounds::eMinR);
    double bmaxr = bounds->get(RadialBounds::eMaxR);

    CHECK_CLOSE_ABS(bminr, rmin, s_epsilon);
    CHECK_CLOSE_ABS(bmaxr, rmax, s_epsilon);

    // Check if the surface is the (negative) identity
    GeometryVisualization::drawSurface(objVis, *disc, tgContext,
                                       Transform3D::Identity(), 80, false,
                                       {80, 80, 120});
    const Vector3D center = disc->center(tgContext);
    GeometryVisualization::drawArrowForward(
        objVis, center, center + 1.2 * rmax * Vector3D::UnitX(), 0.2, 4., 2.5,
        72, {200, 0, 0});
    GeometryVisualization::drawArrowForward(
        objVis, center, center + 1.2 * rmax * Vector3D::UnitY(), 0.2, 4., 2.5,
        72, {0, 200, 0});
    GeometryVisualization::drawArrowForward(
        objVis, center, center + 1.2 * hz * Vector3D::UnitZ(), 0.2, 4., 2.5, 72,
        {0, 0, 200});
    objVis.write("TGeoConversion_TGeoTube_DiscSurface_" +
                 std::to_string(idisc));
    objVis.clear();

    if (idisc < 2) {
      auto discSegment = TGeoSurfaceConverter::toSurface(
          *vols->GetShape(), *gGeoIdentity, axes, 1);
      BOOST_TEST(discSegment != nullptr);
      BOOST_TEST(discSegment->type() == Surface::Disc);

      auto boundsSegment =
          dynamic_cast<const RadialBounds *>(&(discSegment->bounds()));
      BOOST_TEST(boundsSegment != nullptr);
      bminr = boundsSegment->get(RadialBounds::eMinR);
      bmaxr = boundsSegment->get(RadialBounds::eMaxR);
      double hphi = boundsSegment->get(RadialBounds::eHalfPhiSector);
      double mphi = boundsSegment->get(RadialBounds::eAveragePhi);
      CHECK_CLOSE_ABS(bminr, rmin, s_epsilon);
      CHECK_CLOSE_ABS(bmaxr, rmax, s_epsilon);
      CHECK_CLOSE_ABS(hphi, 0.25 * M_PI, s_epsilon);
      CHECK_CLOSE_ABS(mphi, 0., s_epsilon);
      GeometryVisualization::drawSurface(objVis, *discSegment, tgContext,
                                         Transform3D::Identity(), 80, false,
                                         {80, 80, 120});
      GeometryVisualization::drawArrowForward(
          objVis, center, center + 1.2 * bmaxr * Vector3D::UnitX(), 0.2, 4.,
          2.5, 72, {200, 0, 0});
      GeometryVisualization::drawArrowForward(
          objVis, center, center + 1.2 * bmaxr * Vector3D::UnitY(), 0.2, 4.,
          2.5, 72, {0, 200, 0});
      GeometryVisualization::drawArrowForward(
          objVis, center, center + 1.2 * hz * Vector3D::UnitZ(), 0.2, 4., 2.5,
          72, {0, 0, 200});
      objVis.write("TGeoConversion_TGeoTube_DiscSegmentSurface_" +
                   std::to_string(idisc));
      objVis.clear();

    } else {
      BOOST_CHECK_THROW(TGeoSurfaceConverter::toSurface(*vols->GetShape(),
                                                        *gGeoIdentity, axes, 1),
                        std::invalid_argument);
    }
    ++idisc;
  }

  // Check exceptions for not allowed axis definition
  for (const auto &naxes : notAllowedAxes) {
    BOOST_CHECK_THROW(TGeoSurfaceConverter::toSurface(*vol->GetShape(),
                                                      *gGeoIdentity, naxes, 1),
                      std::invalid_argument);
  }
}

}  // namespace Test

}  // namespace Acts
