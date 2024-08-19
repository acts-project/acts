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
#include "Acts/Surfaces/LineBounds.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/TrapezoidBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <typeinfo>

#include <GeoModelKernel/GeoBox.h>
#include <GeoModelKernel/GeoTube.h>
#include <GeoModelKernel/GeoSimplePolygonBrep.h>
#include <GeoModelKernel/GeoFullPhysVol.h>
#include <GeoModelKernel/GeoLogVol.h>
#include <GeoModelKernel/GeoMaterial.h>
#include <GeoModelKernel/GeoTube.h>

BOOST_AUTO_TEST_SUITE(GeoModelDetObj)

BOOST_AUTO_TEST_CASE(GeoModelDetectorObjectFactory) {
  // define materials
  auto material = new GeoMaterial("Material", 1.0);
  auto al = new GeoMaterial("Aluminium", 1.0);

  //define dimensions
  double gm_box_hlx = 100, gm_box_hly = 200, gm_box_hlz = 2;
  double gm_tube_rmin = 5, gm_tube_rmax =5, gm_tube_hlz = 100;
  double gm_rsurf_hlx = 100, gm_rsurf_hly = 200, gm_rsurf_hlz = 2;
  double gm_poly_z = 2;
  std::vector<double> x_verts = {-153, 153,-123, 123};
  std::vector<double> y_verts = {50, 50, -50, -50};
  

  //create shapes
  auto boxXY = new GeoBox(gm_box_hlx, gm_box_hly, gm_box_hlz);
  auto tube = new GeoTube(gm_tube_rmin, gm_tube_rmax, gm_tube_hlz);
  auto ssurface = new GeoBox(gm_rsurf_hlx, gm_rsurf_hly, gm_rsurf_hlz);
  auto poly = new GeoSimplePolygonBrep(gm_poly_z);
  for(int i=0;i<x_verts.size();i++){
    poly->addVertex(x_verts[i], y_verts[i]);
  }


  // create logvols
  auto logXY = new GeoLogVol("LogVolumeXY", boxXY, material);
  auto logTube = new GeoLogVol("LogTube", tube, al);
  auto logSurface = new GeoLogVol("LogSurface", ssurface, al);
  auto logPoly = new GeoLogVol("LogPoly", poly, al);

  // create physvols
  auto fphysXY = new GeoFullPhysVol(logXY);
  auto physTube = new GeoFullPhysVol(logTube);
  auto physSurface = new GeoFullPhysVol(logSurface);
  auto physPoly = new GeoFullPhysVol(logPoly);

  // build hierarchy
  fphysXY->add(physTube);
  fphysXY->add(physSurface);
  fphysXY->add(physPoly);

  PVConstLink physVol{fphysXY};
  auto rBounds = std::make_shared<Acts::RectangleBounds>(100, 200);

  // create pars for conversion
  Acts::GeoModelDetectorObjectFactory::Config gmConfig;
  gmConfig.convertBox = {"LogVolumeXY"};
  Acts::GeometryContext gContext;
  Acts::GeoModelDetectorObjectFactory::Cache gmCache;

  // create factory instance
  Acts::GeoModelDetectorObjectFactory factory =
      Acts::GeoModelDetectorObjectFactory(gmConfig);

  // convert GeoFullPhysVol
  factory.convertFpv("LogVolumeXY", fphysXY, gmCache, gContext);

  // checking the dimension of the converted bounding boxes
  for (auto box : gmCache.boundingBoxes) {
    const Acts::VolumeBounds& bounds = box->volumeBounds();
    BOOST_CHECK(gm_box_hlx == bounds.values()[0]);
    BOOST_CHECK(gm_box_hly == bounds.values()[1]);
    BOOST_CHECK(gm_box_hlz == bounds.values()[2]);
    std::vector<const Acts::Surface*> surfaces = box->surfaces();

    for (auto surface : surfaces) {
      const Acts::SurfaceBounds& sbounds = surface->bounds();
      //check straws
      if (surface->type() == Acts::Surface::SurfaceType::Straw) {

      }
      // rectangle Surface check corner position without trf
      if (sbounds.type() == Acts::SurfaceBounds::eRectangle) {
        //std::cout << surface->type() << std::endl;
        double csxmin = sbounds.values()[0];
        double csymin = sbounds.values()[1];
        double csxmax = sbounds.values()[2];
        double csymax = sbounds.values()[3];
        BOOST_CHECK(gm_rsurf_hlx == -csxmin);
        BOOST_CHECK(gm_rsurf_hly == -csymin);
        BOOST_CHECK(gm_rsurf_hlx == csxmax);
        BOOST_CHECK(gm_rsurf_hly == csymax);
      }
      // trap Surface without trf
      if (sbounds.type() == Acts::SurfaceBounds::eTrapezoid) {
        const auto* trapBounds = dynamic_cast<const Acts::TrapezoidBounds*>(&sbounds);
        std::vector<Acts::Vector2>verts = trapBounds->vertices();
        BOOST_CHECK(verts[3][0]==x_verts[0]);
        BOOST_CHECK(verts[3][1]==y_verts[0]);
        BOOST_CHECK(verts[2][0]==x_verts[1]);
        BOOST_CHECK(verts[2][1]==y_verts[1]);
        BOOST_CHECK(verts[0][0]==x_verts[2]);
        BOOST_CHECK(verts[0][1]==y_verts[2]);
        BOOST_CHECK(verts[1][0]==x_verts[3]);
        BOOST_CHECK(verts[1][1]==y_verts[3]);
      }
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()
