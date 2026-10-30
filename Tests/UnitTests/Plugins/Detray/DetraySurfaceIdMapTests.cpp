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
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsPlugins/Detray/DetrayGeometryConverter.hpp"
#include "ActsPlugins/Detray/DetrayPayloadConverter.hpp"
#include "ActsPlugins/Detray/DetraySurfaceIdMap.hpp"
#include "ActsTests/CommonHelpers/CylindricalTrackingGeometry.hpp"

#include <cstddef>
#include <memory>
#include <stdexcept>

#include <detray/core/detector.hpp>
#include <detray/detectors/default_metadata.hpp>
#include <vecmem/memory/host_memory_resource.hpp>

using namespace Acts;
using namespace ActsPlugins;

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(DetrayConversion)

BOOST_AUTO_TEST_CASE(DetraySurfaceIdMapInsertAndLookup) {
  DetraySurfaceIdMap map;
  BOOST_CHECK(map.empty());

  const auto id1 = GeometryIdentifier().withVolume(1).withSensitive(1);
  const auto id2 = GeometryIdentifier().withVolume(1).withSensitive(2);
  map.insert(10, id1);
  map.insert(11, id2);

  BOOST_CHECK_EQUAL(map.size(), 2u);
  BOOST_CHECK_EQUAL(map.toDetray(id1).value(), 10u);
  BOOST_CHECK_EQUAL(map.toDetray(id2).value(), 11u);
  BOOST_CHECK_EQUAL(map.fromDetray(10).value(), id1);
  BOOST_CHECK_EQUAL(map.fromDetray(11).value(), id2);

  // Unknown surfaces
  BOOST_CHECK(!map.toDetray(GeometryIdentifier().withVolume(2)).has_value());
  BOOST_CHECK(!map.fromDetray(12).has_value());

  // Duplicates on either side are rejected
  BOOST_CHECK_THROW(map.insert(10, GeometryIdentifier().withVolume(3)),
                    std::invalid_argument);
  BOOST_CHECK_THROW(map.insert(12, id1), std::invalid_argument);
  BOOST_CHECK_EQUAL(map.size(), 2u);
}

BOOST_AUTO_TEST_CASE(DetraySurfaceIdMapFromConversion) {
  using detector_t =
      detray::host::detector<detray::default_metadata<detray::array<double>>>;

  auto gctx = GeometryContext::dangerouslyDefaultConstruct();

  CylindricalTrackingGeometry cGeometry(gctx, true);
  std::shared_ptr<const TrackingGeometry> tGeometry =
      cGeometry(*getDefaultLogger("Geo", Logging::INFO));

  DetrayPayloadConverter::Config payloadCfg;
  tGeometry->apply([&payloadCfg](const TrackingVolume& volume) {
    if (volume.volumeName() == "Beampipe") {
      payloadCfg.beampipeVolume = &volume;
    }
  });

  DetrayGeometryConverter::Config cfg;
  cfg.payloadConverter = std::make_shared<DetrayPayloadConverter>(
      payloadCfg, getDefaultLogger("PayloadCnv", Logging::INFO));
  DetrayGeometryConverter converter(cfg,
                                    getDefaultLogger("GeoCnv", Logging::INFO));

  vecmem::host_memory_resource mr;
  auto geometry = converter.convert<detector_t::metadata>(mr, gctx, tGeometry);
  const detector_t& detector = *geometry.detector;
  const DetraySurfaceIdMap& map = geometry.surfaceIds;

  // Every sensitive detray surface is mapped to the ACTS surface it was
  // converted from, and back
  std::size_t nSensitive = 0;
  for (const auto& surface : detector.surfaces()) {
    if (!surface.is_sensitive()) {
      BOOST_CHECK(!map.fromDetray(surface.identifier().value()).has_value());
      continue;
    }
    ++nSensitive;
    const GeometryIdentifier id(surface.source);
    BOOST_CHECK_EQUAL(map.fromDetray(surface.identifier().value()).value(), id);
    BOOST_CHECK_EQUAL(map.toDetray(id).value(), surface.identifier().value());
  }
  BOOST_CHECK_GT(nSensitive, 0u);
  BOOST_CHECK_EQUAL(map.size(), nSensitive);

  // ... and every sensitive ACTS surface is in the map
  std::size_t nActsSensitive = 0;
  tGeometry->visitSurfaces([&](const Surface* surface) {
    ++nActsSensitive;
    BOOST_CHECK(map.toDetray(surface->geometryId()).has_value());
  });
  BOOST_CHECK_EQUAL(nActsSensitive, map.size());

  // Building the map from the detector directly gives the same result
  const auto rebuilt = DetraySurfaceIdMap::fromDetector(detector);
  BOOST_CHECK_EQUAL(rebuilt.size(), map.size());
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
