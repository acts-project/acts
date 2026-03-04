// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/DiamondBounds.hpp"
#include "Acts/Surfaces/LineBounds.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/StrawSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/TrapezoidBounds.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsPlugins/GeoModel/GeoModelDetectorObjectFactory.hpp"
#include "ActsPlugins/GeoModel/GeoModelReader.hpp"

#include <typeinfo>

#include <GeoModelKernel/GeoBox.h>
#include <GeoModelKernel/GeoFullPhysVol.h>
#include <GeoModelKernel/GeoLogVol.h>
#include <GeoModelKernel/GeoMaterial.h>
#include <GeoModelKernel/GeoSimplePolygonBrep.h>
#include <GeoModelKernel/GeoTrd.h>
#include <GeoModelKernel/GeoTube.h>

using namespace Acts;
using namespace ActsPlugins;

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(GeoModelSuite)

BOOST_AUTO_TEST_CASE(GeoModelDetectorObjectFactory) {
  auto al = make_intrusive<GeoMaterial>("Aluminium", 1.0);

  std::vector<std::vector<double>> trapVerts = {
      {-103, -50}, {103, -50}, {183, 50}, {-183, 50}};
  std::vector<std::vector<double>> polyVerts = {
      {-60, -50}, {60, -50}, {153, 0}, {123, 50}, {-123, 50}, {-153, 0}};
  std::vector<std::vector<double>> errVerts = {
      {60, -50}, {153, 0}, {123, 50}, {-123, 50}, {-153, 0}};
  double poly_z = 2;

  auto trap = make_intrusive<GeoSimplePolygonBrep>(poly_z);
  for (const auto& tVert : trapVerts) {
    trap->addVertex(tVert[0], tVert[1]);
  }
  auto poly = make_intrusive<GeoSimplePolygonBrep>(poly_z);
  for (const auto& pVert : polyVerts) {
    poly->addVertex(pVert[0], pVert[1]);
  }
  auto err = make_intrusive<GeoSimplePolygonBrep>(poly_z);
  for (const auto& eVert : errVerts) {
    err->addVertex(eVert[0], eVert[1]);
  }
  auto logTrap = make_intrusive<GeoLogVol>("LogTrap", trap, al);
  auto logPoly = make_intrusive<GeoLogVol>("LogPoly", poly, al);
  auto logErr = make_intrusive<GeoLogVol>("LogErr", err, al);

  auto physTrap = make_intrusive<GeoFullPhysVol>(logTrap);
  auto physPoly = make_intrusive<GeoFullPhysVol>(logPoly);
  auto physErr = make_intrusive<GeoFullPhysVol>(logErr);
  // create pars for conversion
  ActsPlugins::GeoModelDetectorObjectFactory::Config gmConfig;
  auto gContext = GeometryContext::dangerouslyDefaultConstruct();
  ActsPlugins::GeoModelDetectorObjectFactory::Cache trapCache;
  ActsPlugins::GeoModelDetectorObjectFactory::Cache polyCache;
  ActsPlugins::GeoModelDetectorObjectFactory::Cache errCache;

  // create factory instance
  ActsPlugins::GeoModelDetectorObjectFactory factory(gmConfig);

  // convert GeoFullPhysVol (to surfaces)
  factory.convertFpv("Trap", physTrap, trapCache, gContext);
  factory.convertFpv("Poly", physPoly, polyCache, gContext);
  BOOST_CHECK_THROW(factory.convertFpv("Error", physErr, errCache, gContext),
                    std::runtime_error);

  ActsPlugins::GeoModelSensitiveSurface trapSensSurface =
      trapCache.sensitiveSurfaces[0];
  ActsPlugins::GeoModelSensitiveSurface polySensSurface =
      polyCache.sensitiveSurfaces[0];
  std::shared_ptr<Surface> polySurface = std::get<1>(polySensSurface);
  std::shared_ptr<Surface> trapSurface = std::get<1>(trapSensSurface);

  const auto* polyBounds =
      dynamic_cast<const DiamondBounds*>(&polySurface->bounds());
  std::vector<Vector2> convPolyVerts = polyBounds->vertices();
  for (std::size_t i = 0; i < polyVerts.size(); i++) {
    BOOST_CHECK(polyVerts[i][0] == convPolyVerts[i][0]);
    BOOST_CHECK(polyVerts[i][1] == convPolyVerts[i][1]);
  }

  const auto* trapBounds =
      dynamic_cast<const TrapezoidBounds*>(&trapSurface->bounds());
  std::vector<Vector2> convTrapVerts = trapBounds->vertices();
  for (std::size_t i = 0; i < trapVerts.size(); i++) {
    BOOST_CHECK(trapVerts[i][0] == convTrapVerts[i][0]);
    BOOST_CHECK(trapVerts[i][1] == convTrapVerts[i][1]);
  }
}
BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
