// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "ActsTests/CommonHelpers/CubicTrackingGeometry.hpp"
#include "ActsTests/CommonHelpers/CylindricalTrackingGeometry.hpp"

#include <optional>

using namespace Acts;

namespace ActsTests {

// Create a test context
GeometryContext tgContext = GeometryContext::dangerouslyDefaultConstruct();

BOOST_AUTO_TEST_SUITE(GeometrySuite)

BOOST_AUTO_TEST_CASE(CylindricalTrackingGeometryTest) {
  CylindricalTrackingGeometry cGeometry(tgContext);
  auto tGeometry = cGeometry();
  BOOST_CHECK_NE(tGeometry, nullptr);

  BOOST_CHECK(tGeometry->geometryVersion() ==
              TrackingGeometry::GeometryVersion::Gen1);
}

BOOST_AUTO_TEST_CASE(CubicTrackingGeometryTest) {
  CubicTrackingGeometry cGeometry(tgContext);
  auto tGeometry = cGeometry();
  BOOST_CHECK_NE(tGeometry, nullptr);
  BOOST_CHECK(tGeometry->geometryVersion() ==
              TrackingGeometry::GeometryVersion::Gen1);
}

BOOST_AUTO_TEST_CASE(SurfaceLookupTracksNonSensitiveSurfaces) {
  CubicTrackingGeometry builder(tgContext);
  auto geometry = builder();
  BOOST_REQUIRE_NE(geometry, nullptr);

  std::optional<GeometryIdentifier> nonSensitiveSurfaceId;
  geometry->visitSurfaces(
      [&](const Surface* surface) {
        if (surface == nullptr) {
          return;
        }
        auto geoId = surface->geometryId();
        if (geoId == GeometryIdentifier{}) {
          return;
        }
        if (geoId.sensitive() == 0u) {
          nonSensitiveSurfaceId = geoId;
        }
      },
      false);

  BOOST_REQUIRE(nonSensitiveSurfaceId.has_value());

  BOOST_CHECK_MESSAGE(
      geometry->findSurface(*nonSensitiveSurfaceId) != nullptr,
      "TrackingGeometry::findSurface should return the surface with ID "
          << *nonSensitiveSurfaceId << " even when it is not sensitive");
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
