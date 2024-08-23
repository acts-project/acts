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
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <typeinfo>

#include <GeoModelKernel/GeoBox.h>
#include <GeoModelKernel/GeoFullPhysVol.h>
#include <GeoModelKernel/GeoLogVol.h>
#include <GeoModelKernel/GeoMaterial.h>
#include <GeoModelKernel/GeoTube.h>

BOOST_AUTO_TEST_SUITE(GeoModelDetObj)

BOOST_AUTO_TEST_CASE(GeoModelDetectorObjectFactory) {
  // define materials
  auto material = new GeoMaterial("Material", 1.0);
  auto al = new GeoMaterial("Aluminium", 1.0);

  // define dimensions
  double gmhlx = 100, gmhly = 200, gmhlz = 2;
  double gmrmin = 5, gmrmax = 5, gmhlzt = 100;
  double gmhlxs = 100, gmhlys = 200, gmhlzs = 2;

  // create shapes
  auto boxXY = new GeoBox(gmhlx, gmhly, gmhlz);
  auto tube = new GeoTube(gmrmin, gmrmax, gmhlzt);
  auto ssurface = new GeoBox(gmhlxs, gmhlys, gmhlzs);

  // create logvols
  auto logXY = new GeoLogVol("LogVolumeXY", boxXY, material);
  auto logTube = new GeoLogVol("LogTube", tube, al);
  auto logSurface = new GeoLogVol("LogSurface", ssurface, al);

  // create physvols
  auto fphysXY = new GeoFullPhysVol(logXY);
  auto physTube = new GeoFullPhysVol(logTube);
  auto physSurface = new GeoFullPhysVol(logSurface);

  // build hierarchy
  fphysXY->add(physTube);
  fphysXY->add(physSurface);

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
    BOOST_CHECK(gmhlx == bounds.values()[0]);
    BOOST_CHECK(gmhly == bounds.values()[1]);
    BOOST_CHECK(gmhlz == bounds.values()[2]);
    std::vector<const Acts::Surface*> surfaces = box->surfaces();

    for (auto surface : surfaces) {
      const Acts::SurfaceBounds& sbounds = surface->bounds();
      // Straw check outer radius and length without trf
      if (surface->type() == Acts::Surface::SurfaceType::Straw) {
        BOOST_CHECK(sbounds.values()[0] == gmrmax);
        BOOST_CHECK(sbounds.values()[1] == gmhlzt);
      }

      // plane Surface check corner position without trf
      if (surface->type() == Acts::Surface::SurfaceType::Plane) {
        double csxmin = sbounds.values()[0];
        double csymin = sbounds.values()[1];
        double csxmax = sbounds.values()[2];
        double csymax = sbounds.values()[3];
        BOOST_CHECK(gmhlxs == -csxmin);
        BOOST_CHECK(gmhlys == -csymin);
        BOOST_CHECK(gmhlxs == csxmax);
        BOOST_CHECK(gmhlys == csymax);
      }
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()
