// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsPlugins/Detray/DetrayGeometryConverter.hpp"
#include "ActsPlugins/Detray/DetrayPayloadConverter.hpp"
#include "ActsPlugins/Detray/DetraySurfaceIdMap.hpp"
#include "ActsTests/CommonHelpers/CylindricalTrackingGeometry.hpp"

#include <cstddef>
#include <memory>
#include <stdexcept>
#include <vector>

#include <detray/core/detector.hpp>
#include <detray/detectors/default_metadata.hpp>
#include <detray/geometry/identifier.hpp>
#include <vecmem/memory/host_memory_resource.hpp>

using namespace Acts;
using namespace ActsPlugins;

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(DetrayConversion)

BOOST_AUTO_TEST_CASE(DetraySurfaceIdMapFromEntries) {
  using detray::geometry::identifier;

  const DetraySurfaceIdMap emptyMap;
  BOOST_CHECK(emptyMap.empty());

  auto makeSurface = [](GeometryIdentifier id) {
    auto surface = Surface::makeShared<PlaneSurface>(
        Transform3::Identity(), std::make_shared<RectangleBounds>(1., 1.));
    surface->assignGeometryId(id);
    return surface;
  };

  const auto id1 = GeometryIdentifier().withVolume(1).withSensitive(1);
  const auto id2 = GeometryIdentifier().withVolume(1).withSensitive(2);
  const auto surface1 = makeSurface(id1);
  const auto surface2 = makeSurface(id2);

  const std::vector<DetraySurfaceIdMap::Entry> entries = {
      {identifier{10}, surface1.get()}, {identifier{11}, surface2.get()}};
  const DetraySurfaceIdMap map(entries);

  BOOST_CHECK_EQUAL(map.size(), 2u);
  BOOST_CHECK(map.toDetray(id1).value() == identifier{10});
  BOOST_CHECK(map.toDetray(id2).value() == identifier{11});
  BOOST_CHECK_EQUAL(map.fromDetray(identifier{10}).value(), id1);
  BOOST_CHECK_EQUAL(map.fromDetray(identifier{11}).value(), id2);
  BOOST_CHECK_EQUAL(map.surface(identifier{10}), surface1.get());
  BOOST_CHECK_EQUAL(map.surface(identifier{11}), surface2.get());

  // Unknown surfaces
  BOOST_CHECK(!map.toDetray(GeometryIdentifier().withVolume(2)).has_value());
  BOOST_CHECK(!map.fromDetray(identifier{12}).has_value());
  BOOST_CHECK_EQUAL(map.surface(identifier{12}), nullptr);

  // Duplicates on either side and null surfaces are rejected
  const auto surface3 = makeSurface(GeometryIdentifier().withVolume(3));
  const std::vector<DetraySurfaceIdMap::Entry> duplicateDetray = {
      {identifier{10}, surface1.get()}, {identifier{10}, surface3.get()}};
  BOOST_CHECK_THROW(DetraySurfaceIdMap{duplicateDetray}, std::invalid_argument);
  const std::vector<DetraySurfaceIdMap::Entry> duplicateActs = {
      {identifier{10}, surface1.get()}, {identifier{12}, surface1.get()}};
  BOOST_CHECK_THROW(DetraySurfaceIdMap{duplicateActs}, std::invalid_argument);
  const std::vector<DetraySurfaceIdMap::Entry> nullSurface = {
      {identifier{10}, nullptr}};
  BOOST_CHECK_THROW(DetraySurfaceIdMap{nullSurface}, std::invalid_argument);
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
      BOOST_CHECK(!map.fromDetray(surface.identifier()).has_value());
      BOOST_CHECK_EQUAL(map.surface(surface.identifier()), nullptr);
      continue;
    }
    ++nSensitive;
    const GeometryIdentifier id(surface.source);
    BOOST_CHECK_EQUAL(map.fromDetray(surface.identifier()).value(), id);
    BOOST_CHECK(map.toDetray(id).value() == surface.identifier());
    BOOST_CHECK_EQUAL(map.surface(surface.identifier()),
                      tGeometry->findSurface(id));
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
  const DetraySurfaceIdMap rebuilt(detector, *tGeometry);
  BOOST_CHECK_EQUAL(rebuilt.size(), map.size());
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
