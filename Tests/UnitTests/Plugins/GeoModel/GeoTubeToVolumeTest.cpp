// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

// In order to avoid conflicts with declarations in Geomodel that is fixed in
// v3.5
//  clang-formal off
#include "ActsPlugins/GeoModel/GeoModelDetectorObjectFactory.hpp"
// clang-formal on
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/SurfaceBounds.hpp"
#include "Acts/Surfaces/TrapezoidBounds.hpp"
#include "ActsPlugins/GeoModel/GeoModelConverters.hpp"
#include "ActsTests/CommonHelpers/FloatComparisons.hpp"

#include <GeoModelKernel/GeoFullPhysVol.h>
#include <GeoModelKernel/GeoLogVol.h>
#include <GeoModelKernel/GeoMaterial.h>
#include <GeoModelKernel/GeoTrap.h>
#include <GeoModelKernel/GeoTube.h>
#include <GeoModelKernel/GeoVPhysVol.h>

using namespace Acts;
using namespace ActsPlugins;

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(GeoModelPlugin)

// GeoBox conversion test case
BOOST_AUTO_TEST_CASE(GeoBoxToSensitiveConversion) {
  auto material = make_intrusive<GeoMaterial>("Material", 1.0);
  // Let's create a GeoFullPhysVol object

  // (BOX object) - XY
  std::vector<double> dims = {5, 6, 50};
  auto tube = make_intrusive<GeoTube>(dims[0], dims[1], dims[2]);
  auto logTube = make_intrusive<GeoLogVol>("Tube", tube, material);
  auto physTube = make_intrusive<GeoFullPhysVol>(logTube);

  // create pars for conversion
  GeoModelDetectorObjectFactory::Config gmConfig;
  gmConfig.convertBox = {"Tube"};
  auto gContext = GeometryContext::dangerouslyDefaultConstruct();
  GeoModelDetectorObjectFactory::Cache gmCache;

  // create factory instance
  GeoModelDetectorObjectFactory factory(gmConfig);

  factory.convertFpv("Tube", physTube, gmCache, gContext);
  BOOST_CHECK(!gmCache.volumeBoxFPVs.empty());
  const auto& volumeTube = gmCache.volumeBoxFPVs[0].volume;
  const auto* bounds =
      dynamic_cast<const CylinderVolumeBounds*>(&volumeTube->volumeBounds());
  std::vector<double> convDims = bounds->values();
  for (std::size_t i = 0; i < dims.size(); i++) {
    BOOST_CHECK(dims[i] == convDims[i]);
  }
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
