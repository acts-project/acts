// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsPlugins/Detray/DetrayConversionUtils.hpp"
#include "ActsPlugins/Detray/DetrayPayloadConverter.hpp"

#include <memory>
#include <string>
#include <unordered_map>

#include <detray/builders/detector_builder.hpp>
#include <detray/geometry/identifier.hpp>
#include <detray/io/backend/geometry_reader.hpp>
#include <detray/io/backend/homogeneous_material_reader.hpp>
#include <detray/io/backend/material_map_reader.hpp>
#include <detray/io/backend/surface_grid_reader.hpp>
#include <vecmem/memory/memory_resource.hpp>

namespace ActsPlugins {

/// @ingroup detray_plugin
/// @brief Converter from an ACTS TrackingGeometry to a detray detector
///
/// The geometry conversion is a two-step process: first the ACTS geometry is
/// converted into detray payloads by a @c DetrayPayloadConverter, then those
/// payloads are used to build the actual detray detector. This class drives
/// the second step and lets the call site fully customize the first step by
/// supplying its own configured @c DetrayPayloadConverter instance through the
/// @c Config.
class DetrayGeometryConverter {
 public:
  /// @brief Configuration for the geometry converter
  struct Config {
    /// The payload converter used to turn the ACTS geometry into detray
    /// payloads. Supplying it here allows the call site to fully customize the
    /// payload conversion (e.g. the beampipe volume, sensitive surface
    /// strategy or the navigation/material dispatchers).
    std::shared_ptr<const DetrayPayloadConverter> payloadConverter;

    /// Whether to convert material information from ACTS to detray
    bool convertMaterial = true;

    /// Whether to convert surface grid information from ACTS to detray
    bool convertSurfaceGrids = true;
  };

  /// @brief Combined result of a geometry conversion
  ///
  /// Bundles the built detray detector together with the bookkeeping that
  /// relates it back to the tracking geometry: the detray volume/surface name
  /// map and the map from detray surface identifiers to the source surfaces.
  /// @tparam metadata_t the detector metadata type
  template <typename metadata_t>
  struct DetrayGeometry {
    /// The built detray detector
    std::shared_ptr<detray::detector<metadata_t>> detector;

    /// The detray volume and surface name map
    detray::name_map names;

    /// Map from detray surface identifiers to the source surfaces for all
    /// sensitive surfaces in the detector. The GeometryIdentifier can be
    /// obtained from the surface via @c Surface::geometryId.
    std::unordered_map<detray::geometry::identifier, const Acts::Surface*>
        detrayToSurfaceMap;

    /// The tracking geometry the detector was converted from, retained to keep
    /// the surfaces referenced by @c detrayToSurfaceMap alive
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry;
  };

  /// Constructor
  /// @param config Configuration object
  /// @param logger Logger instance
  explicit DetrayGeometryConverter(
      Config config,
      std::unique_ptr<const Acts::Logger> logger = Acts::getDefaultLogger(
          "DetrayGeometryConverter", Acts::Logging::INFO))
      : m_cfg(std::move(config)), m_logger(std::move(logger)) {
    if (m_cfg.payloadConverter == nullptr) {
      throw std::invalid_argument(
          "DetrayGeometryConverter: payloadConverter must be set");
    }
  }

  /// @brief Convert an ACTS TrackingGeometry into a detray detector
  /// @tparam metadata_t the detector metadata type to build
  ///
  /// @param mr the memory resource to use for the detray detector construction
  /// @param gctx the geometry context
  /// @param trackingGeometry the ACTS tracking geometry to convert
  /// @param detectorName the name to set for the detray detector (optional, if
  ///     not set, it will be taken from the payloads or defaulted to empty)
  ///
  /// This method performs the following steps:
  /// 1. It converts the ACTS tracking geometry into detray payloads using the
  ///    configured DetrayPayloadConverter.
  /// 2. It builds a detray detector from the converted payloads using the
  ///    detray::detector_builder.
  /// 3. It builds the map from detray surface identifiers back to the source
  ///    surfaces.
  ///
  /// @return The built detray detector together with its name map, the
  ///     detray->surface map and the source tracking geometry.
  template <typename metadata_t>
  DetrayGeometry<metadata_t> convert(
      vecmem::memory_resource& mr, const Acts::GeometryContext& gctx,
      std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry,
      const std::string& detectorName = "") const {
    using detector_t = detray::detector<metadata_t>;

    if (trackingGeometry == nullptr) {
      throw std::invalid_argument(
          "DetrayGeometryConverter: trackingGeometry must not be null");
    }

    // ── Convert TrackingGeometry → detray payloads ────────────────────────
    auto payloads = m_cfg.payloadConverter->convertTrackingGeometry(
        gctx, *trackingGeometry);

    // ── Build detray detector from payloads ───────────────────────────────
    detray::detector_builder<metadata_t> detectorBuilder{};

    detray::io::geometry_reader::from_payload<detector_t>(detectorBuilder,
                                                          *payloads.detector);

    if (m_cfg.convertMaterial) {
      detray::io::homogeneous_material_reader::from_payload<detector_t>(
          detectorBuilder, *payloads.homogeneousMaterial);

      detray::io::material_map_reader<std::integral_constant<std::size_t, 2>>::
          from_payload<detector_t>(detectorBuilder,
                                   std::move(*payloads.materialGrids));
    }

    if (m_cfg.convertSurfaceGrids) {
      detray::io::surface_grid_reader<typename detector_t::surface_type,
                                      std::integral_constant<std::size_t, 0>,
                                      std::integral_constant<std::size_t, 2>>::
          template from_payload<detector_t>(detectorBuilder,
                                            *payloads.surfaceGrids);
    }

    if (!detectorName.empty()) {
      detectorBuilder.set_name(detectorName);
    } else if (payloads.names.contains(0)) {
      detectorBuilder.set_name(payloads.names.at(0));
    }

    DetrayGeometry<metadata_t> result{};
    result.detector =
        std::make_shared<detector_t>(detectorBuilder.build(mr, result.names));
    result.detrayToSurfaceMap =
        buildDetrayToSurfaceMap(*result.detector, *trackingGeometry);
    result.trackingGeometry = std::move(trackingGeometry);

    return result;
  }

 private:
  /// Build a mapping from detray surface identifiers to the source surfaces
  /// for all sensitive surfaces in the detray detector.
  /// @tparam detector_t the type of the detray detector
  /// @param detrayDetector the detray detector to build the mapping for
  /// @param trackingGeometry the geometry used to resolve the surfaces
  /// @return a map from detray surface identifiers to the source surfaces
  template <typename detector_t>
  std::unordered_map<detray::geometry::identifier, const Acts::Surface*>
  buildDetrayToSurfaceMap(
      const detector_t& detrayDetector,
      const Acts::TrackingGeometry& trackingGeometry) const {
    std::unordered_map<detray::geometry::identifier, const Acts::Surface*>
        detrayToSurfaceMap;

    for (const auto& surface : detrayDetector.surfaces()) {
      // surface.source is the GeometryIdentifier encoded as uint64
      const Acts::GeometryIdentifier geometryId(surface.source);
      const Acts::Surface* surfacePtr =
          trackingGeometry.findSurface(geometryId);
      if (surfacePtr == nullptr || !surfacePtr->isSensitive()) {
        continue;  // skip portals and passives
      }
      detrayToSurfaceMap[surface.identifier()] = surfacePtr;
    }

    ACTS_INFO("Built detray→surface map with " << detrayToSurfaceMap.size()
                                               << " sensitive surfaces");
    return detrayToSurfaceMap;
  }

  Config m_cfg;

  const Acts::Logger& logger() const { return *m_logger; }
  std::unique_ptr<const Acts::Logger> m_logger;
};

}  // namespace ActsPlugins
