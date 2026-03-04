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
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/TrapezoidVolumeBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/SurfaceBounds.hpp"
#include "Acts/Surfaces/TrapezoidBounds.hpp"
#include "ActsPlugins/GeoModel/GeoModelConverters.hpp"
#include "ActsTests/CommonHelpers/FloatComparisons.hpp"

#include <GeoModelKernel/GeoFullPhysVol.h>
#include <GeoModelKernel/GeoLogVol.h>
#include <GeoModelKernel/GeoMaterial.h>
#include <GeoModelKernel/GeoTrd.h>

using namespace Acts;
using namespace ActsPlugins;

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(GeoModelSuite)

// GeoBox conversion test case
BOOST_AUTO_TEST_CASE(GeoTrdToVolumeConversion) {
  auto material = make_intrusive<GeoMaterial>("Material", 1.0);
  // Let's create a GeoFullPhysVol object
  double geoHlX1 = 2, geoHlX2 = 2, geoHlY1 = 50, geoHlY2 = 80, geoHlZ = 60;
  auto trd = make_intrusive<GeoTrd>(geoHlX1, geoHlX1, geoHlY1, geoHlY2, geoHlZ);
  auto logTrd = make_intrusive<GeoLogVol>("Trd", trd, material);
  auto physTrd = make_intrusive<GeoFullPhysVol>(logTrd);

  // this should produce an error while converting
  auto errTrd = make_intrusive<GeoTrd>(2, 3, 25, 40, 30);
  auto errLogTrd = make_intrusive<GeoLogVol>("Trd", errTrd, material);
  auto errPhysTrd = make_intrusive<GeoFullPhysVol>(errLogTrd);

  // create pars for conversion
  GeoModelDetectorObjectFactory::Config gmConfig;
  gmConfig.convertBox = {"Trd"};
  auto gContext = GeometryContext::dangerouslyDefaultConstruct();
  GeoModelDetectorObjectFactory::Cache gmCache;
  GeoModelDetectorObjectFactory::Cache errCache;

  // create factory instance
  GeoModelDetectorObjectFactory factory(gmConfig);

  // test error case
  BOOST_CHECK_THROW(factory.convertFpv("Trd", errPhysTrd, errCache, gContext),
                    std::invalid_argument);
  factory.convertFpv("Trd", physTrd, gmCache, gContext);

  BOOST_CHECK(!gmCache.volumeBoxFPVs.empty());
  const auto& volumeTrd = gmCache.volumeBoxFPVs[0].volume;
  const auto* bounds =
      dynamic_cast<const TrapezoidVolumeBounds*>(&volumeTrd->volumeBounds());
  std::vector<double> convHls = bounds->values();
  // note: GeoTrd and Acts use different coordinates
  BOOST_CHECK(geoHlX1 == convHls[3]);
  BOOST_CHECK(geoHlX2 == convHls[3]);
  BOOST_CHECK(geoHlY1 == convHls[0]);
  BOOST_CHECK(geoHlY2 == convHls[1]);
  BOOST_CHECK(geoHlZ == convHls[2]);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
