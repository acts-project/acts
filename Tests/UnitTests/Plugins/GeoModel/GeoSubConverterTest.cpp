// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

//clang-format off
#include "ActsPlugins/GeoModel/GeoModelDetectorObjectFactory.hpp"
//clang-format on
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/SurfaceBounds.hpp"
#include "Acts/Surfaces/TrapezoidBounds.hpp"
#include "ActsPlugins/GeoModel/GeoModelConverters.hpp"
#include "ActsTests/CommonHelpers/FloatComparisons.hpp"

#include <GeoModelKernel/GeoFullPhysVol.h>
#include <GeoModelKernel/GeoLogVol.h>
#include <GeoModelKernel/GeoMaterial.h>
#include <GeoModelKernel/GeoShapeSubtraction.h>
#include <GeoModelKernel/GeoTrap.h>
#include <GeoModelKernel/GeoTrd.h>
#include <GeoModelKernel/GeoVPhysVol.h>

using namespace Acts;
using namespace ActsPlugins;

auto tContext = GeometryContext::dangerouslyDefaultConstruct();
RotationMatrix3 idRotation = RotationMatrix3::Identity();
Transform3 idTransform = Transform3::Identity();

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(GeoModelSuite)

// GeoBox conversion test case
BOOST_AUTO_TEST_CASE(GeoSubToSensitiveConversion) {
  auto material = make_intrusive<GeoMaterial>("Material", 1.0);

  // BOX object
  double hlX = 200, hlY = 100, hlZ = 2;
  auto shapeA = make_intrusive<GeoBox>(hlX, hlY, hlZ);
  // Trapezoid object
  auto shapeB = make_intrusive<GeoTrd>(2, 2, 50, 80, 60);

  // create subtraction
  auto geoSub = make_intrusive<GeoShapeSubtraction>(shapeA, shapeB);
  auto logSub = make_intrusive<GeoLogVol>("LogVolume", geoSub, material);
  auto fphysSub = make_intrusive<GeoFullPhysVol>(logSub);

  // create pars for conversion
  GeoModelDetectorObjectFactory::Config gmConfig;
  auto gContext = GeometryContext::dangerouslyDefaultConstruct();
  GeoModelDetectorObjectFactory::Cache subCache;

  // create factory instance
  GeoModelDetectorObjectFactory factory(gmConfig);

  // convert GeoFullPhysVol (to surfaces)
  factory.convertFpv("Sub", fphysSub, subCache, gContext);

  GeoModelSensitiveSurface subSensSurface = subCache.sensitiveSurfaces[0];
  std::shared_ptr<Surface> subSurface = std::get<1>(subSensSurface);
  const auto* subBounds =
      dynamic_cast<const RectangleBounds*>(&subSurface->bounds());
  BOOST_CHECK(subBounds->halfLengthX() == hlX);
  BOOST_CHECK(subBounds->halfLengthY() == hlY);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
