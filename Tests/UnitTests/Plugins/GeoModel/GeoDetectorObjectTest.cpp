// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Plugins/GeoModel/GeoModelDetectorObjectFactory.hpp"
#include "Acts/Plugins/GeoModel/GeoModelReader.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/DiamondBounds.hpp"
#include "Acts/Surfaces/LineBounds.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/StrawSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/TrapezoidBounds.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <typeinfo>

#include <GeoModelKernel/GeoBox.h>
#include <GeoModelKernel/GeoFullPhysVol.h>
#include <GeoModelKernel/GeoLogVol.h>
#include <GeoModelKernel/GeoMaterial.h>
#include <GeoModelKernel/GeoSimplePolygonBrep.h>
#include <GeoModelKernel/GeoTrd.h>
#include <GeoModelKernel/GeoTube.h>

BOOST_AUTO_TEST_SUITE(GeoModelDetObj)

BOOST_AUTO_TEST_CASE(GeoModelDetectorObjectFactory) {
  // define materials
  GeoIntrusivePtr<GeoMaterial> material(new GeoMaterial("Material", 1.0));
  GeoIntrusivePtr<GeoMaterial> al(new GeoMaterial("Aluminium", 1.0));

  // define dimensions
  double gmBoxHlx = 100, gmBoxHly = 200, gmBoxHlz = 2;
  double gmTubeRmin = 5, gmTubeRmax = 6, gmTubeHlz = 100;
  double gmRsurfHlx = 100, gmRsurfHly = 200, gmRsurfHlz = 2;
  double gmPolyZ = 2;
  std::vector<double> trapXVerts = {-103, 103, 183, -183};
  std::vector<double> trapYVerts = {-50, -50, 50, 50};
  std::vector<double> polyXVerts = {-60, 60, 153, 123, -123, -153};
  std::vector<double> polyYVerts = {-50, -50, 0, 50, 50, 0};

  // create shapes
  GeoIntrusivePtr<GeoBox> boxXY(new GeoBox(gmBoxHlx, gmBoxHly, gmBoxHlz));
  GeoIntrusivePtr<GeoTube> tube(new GeoTube(gmTubeRmin, gmTubeRmax, gmTubeHlz));
  GeoIntrusivePtr<GeoBox> ssurface(
      new GeoBox(gmRsurfHlx, gmRsurfHly, gmRsurfHlz));
  double halfX1 = fabs(trapXVerts[0] - trapXVerts[1]) / 2;
  double halfX2 = fabs(trapXVerts[2] - trapXVerts[3]) / 2;
  double halfY1 = fabs(trapYVerts[0] - trapYVerts[2]) / 2;
  GeoIntrusivePtr<GeoTrd> trd(new GeoTrd(1, 1, halfX1, halfX2, halfY1));
  GeoIntrusivePtr<GeoSimplePolygonBrep> trap(new GeoSimplePolygonBrep(gmPolyZ));
  for (long unsigned int i = 0; i < trapXVerts.size(); i++) {
    trap->addVertex(trapXVerts[i], trapYVerts[i]);
  }
  GeoIntrusivePtr<GeoSimplePolygonBrep> poly(new GeoSimplePolygonBrep(gmPolyZ));
  for (long unsigned int i = 0; i < polyXVerts.size(); i++) {
    poly->addVertex(polyXVerts[i], polyYVerts[i]);
  }

  // create logvols
  GeoIntrusivePtr<GeoLogVol> logXY(
      new GeoLogVol("LogVolumeXY", boxXY, material.get()));
  GeoIntrusivePtr<GeoLogVol> logTube(new GeoLogVol("LogTube", tube, al));
  GeoIntrusivePtr<GeoLogVol> logSurface(
      new GeoLogVol("LogSurface", ssurface, al));
  GeoIntrusivePtr<GeoLogVol> logTrap(new GeoLogVol("LogTrap", trap, al));
  GeoIntrusivePtr<GeoLogVol> logTrd(new GeoLogVol("LogTrd", trd, al));
  GeoIntrusivePtr<GeoLogVol> logPoly(new GeoLogVol("LogPoly", poly, al));

  // create physvols
  GeoIntrusivePtr<GeoFullPhysVol> fphysXY(new GeoFullPhysVol(logXY));
  GeoIntrusivePtr<GeoFullPhysVol> physTube(new GeoFullPhysVol(logTube));
  GeoIntrusivePtr<GeoFullPhysVol> physSurface(new GeoFullPhysVol(logSurface));
  GeoIntrusivePtr<GeoFullPhysVol> physTrap(new GeoFullPhysVol(logTrap));
  GeoIntrusivePtr<GeoFullPhysVol> physTrd(new GeoFullPhysVol(logTrd));
  GeoIntrusivePtr<GeoFullPhysVol> physPoly(new GeoFullPhysVol(logPoly));

  // build hierarchy
  fphysXY->add(physTube);
  fphysXY->add(physSurface);
  fphysXY->add(physTrap);
  fphysXY->add(physTrd);
  fphysXY->add(physPoly);

  PVConstLink physVol{fphysXY};
  auto rBounds = std::make_shared<Acts::RectangleBounds>(100, 200);

  // create pars for conversion
  Acts::GeoModelDetectorObjectFactory::Config gmConfig;
  gmConfig.convertBox = {"LogVolumeXY"};
  Acts::GeometryContext gContext;
  Acts::GeoModelDetectorObjectFactory::Cache gmCache;

  // create factory instance
  Acts::GeoModelDetectorObjectFactory factory(gmConfig);

  // convert GeoFullPhysVol
  factory.convertFpv("LogVolumeXY", fphysXY, gmCache, gContext);

  // checking the dimension of the converted bounding boxes
  for (const auto& box : gmCache.boundingBoxes) {
    const Acts::VolumeBounds& bounds = box->volumeBounds();
    BOOST_CHECK(gmBoxHlx == bounds.values()[0]);
    BOOST_CHECK(gmBoxHly == bounds.values()[1]);
    BOOST_CHECK(gmBoxHlz == bounds.values()[2]);
    std::vector<const Acts::Surface*> surfaces = box->surfaces();

    for (auto surface : surfaces) {
      const Acts::SurfaceBounds& sbounds = surface->bounds();
      // check straws
      if (surface->type() == Acts::Surface::SurfaceType::Straw) {
        const auto* lineBounds =
            dynamic_cast<const Acts::LineBounds*>(&sbounds);
        BOOST_CHECK(gmTubeRmax == lineBounds->get(Acts::LineBounds::eR));
        BOOST_CHECK(gmTubeHlz ==
                    lineBounds->get(Acts::LineBounds::eHalfLengthZ));
      }
      // rectangle Surface check corner position without trf
      if (sbounds.type() == Acts::SurfaceBounds::eRectangle) {
        double csxmin = sbounds.values()[0];
        double csymin = sbounds.values()[1];
        double csxmax = sbounds.values()[2];
        double csymax = sbounds.values()[3];
        BOOST_CHECK(gmRsurfHlx == -csxmin);
        BOOST_CHECK(gmRsurfHly == -csymin);
        BOOST_CHECK(gmRsurfHlx == csxmax);
        BOOST_CHECK(gmRsurfHly == csymax);
      }
      // trap Surface without trf
      if (sbounds.type() == Acts::SurfaceBounds::eTrapezoid) {
        const auto* trapBounds =
            dynamic_cast<const Acts::TrapezoidBounds*>(&sbounds);
        std::vector<Acts::Vector2> trapVerts = trapBounds->vertices();
        for (long unsigned int i = 0; i < trapVerts.size(); i++) {
          BOOST_CHECK(trapVerts[i][0] == trapXVerts[i]);
          BOOST_CHECK(trapVerts[i][1] == trapYVerts[i]);
        }
      }
      if (sbounds.type() == Acts::SurfaceBounds::eDiamond) {
        const auto* polyBounds =
            dynamic_cast<const Acts::DiamondBounds*>(&sbounds);
        std::vector<Acts::Vector2> polyVerts = polyBounds->vertices();
        for (long unsigned int i = 0; i < polyVerts.size(); i++) {
          BOOST_CHECK(polyVerts[i][0] == polyXVerts[i]);
          BOOST_CHECK(polyVerts[i][1] == polyYVerts[i]);
        }
      }
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()
