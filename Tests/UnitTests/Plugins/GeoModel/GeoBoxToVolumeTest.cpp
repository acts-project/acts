// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

// switching format off to avoid conflicting declaration in GeoModel
// needed until Acts GeoModel bumps to 6.5
//clang-format off
#include "ActsPlugins/GeoModel/GeoModelDetectorObjectFactory.hpp"
//clang-format on
#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/SurfaceBounds.hpp"
#include "Acts/Surfaces/TrapezoidBounds.hpp"
#include "ActsPlugins/GeoModel/GeoModelConverters.hpp"
#include "ActsTests/CommonHelpers/FloatComparisons.hpp"

#include <GeoModelKernel/GeoBox.h>
#include <GeoModelKernel/GeoFullPhysVol.h>
#include <GeoModelKernel/GeoLogVol.h>
#include <GeoModelKernel/GeoMaterial.h>
#include <GeoModelKernel/GeoTrap.h>
#include <GeoModelKernel/GeoTrd.h>
#include <GeoModelKernel/GeoVPhysVol.h>

using namespace Acts;
using namespace ActsPlugins;

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(GeoModelSuite)

// GeoBox conversion test case
BOOST_AUTO_TEST_CASE(GeoBoxToSensitiveConversion) {
  auto material = make_intrusive<GeoMaterial>("Material", 1.0);
  // Let's create a GeoFullPhysVol object

  // (BOX object) - XY
  std::vector<double> hls = {100, 200, 2};
  auto box = make_intrusive<GeoBox>(hls[0], hls[1], hls[2]);
  auto logBox = make_intrusive<GeoLogVol>("Box", box, material);
  auto physBox = make_intrusive<GeoFullPhysVol>(logBox);

  // create pars for conversion
  GeoModelDetectorObjectFactory::Config gmConfig;
  gmConfig.convertBox = {"Box"};
  auto gContext = GeometryContext::dangerouslyDefaultConstruct();
  GeoModelDetectorObjectFactory::Cache gmCache;

  // create factory instance
  GeoModelDetectorObjectFactory factory(gmConfig);

  factory.convertFpv("Box", physBox, gmCache, gContext);
  BOOST_CHECK(!gmCache.volumeBoxFPVs.empty());
  const auto& volumeBox = gmCache.volumeBoxFPVs[0].volume;
  const auto* bounds =
      dynamic_cast<const CuboidVolumeBounds*>(&volumeBox->volumeBounds());
  std::vector<double> convHls = bounds->values();
  for (std::size_t i = 0; i < hls.size(); i++) {
    BOOST_CHECK(hls[i] == convHls[i]);
  }
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
