// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "ActsPlugins/GeoModel/GeoModelDetectorElement.hpp"

#include <GeoModelKernel/GeoBox.h>
#include <GeoModelKernel/GeoFullPhysVol.h>
#include <GeoModelKernel/GeoLogVol.h>
#include <GeoModelKernel/GeoMaterial.h>

using namespace Acts;
using namespace ActsPlugins;

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(GeoModelSuite)

BOOST_AUTO_TEST_CASE(GeoModelDetectorElementConstruction) {
  auto material = make_intrusive<GeoMaterial>("Material", 1.0);
  // Let's create a GeoFullPhysVol object

  // (BOX object)
  auto boxXY = make_intrusive<GeoBox>(100, 200, 2);
  auto logXY = make_intrusive<GeoLogVol>("LogVolumeXY", boxXY, material);
  auto fphysXY = make_intrusive<GeoFullPhysVol>(logXY);
  auto rBounds = std::make_shared<RectangleBounds>(100, 200);

  PVConstLink physXY{fphysXY};
  auto elementXY =
      GeoModelDetectorElement::createDetectorElement<PlaneSurface,
                                                     RectangleBounds>(
          physXY, rBounds, Transform3::Identity(), 2.0);

  const Surface& surface = elementXY->surface();
  BOOST_CHECK(surface.type() == Surface::SurfaceType::Plane);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
