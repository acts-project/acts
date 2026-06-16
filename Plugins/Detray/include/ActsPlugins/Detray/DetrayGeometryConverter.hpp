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
#include "ActsPlugins/Detray/DetrayConversionUtils.hpp"
#include "ActsPlugins/Detray/DetrayPayloadConverter.hpp"

#include <filesystem>
#include <fstream>

#include <detray/builders/detector_builder.hpp>
#include <detray/io/backend/geometry_reader.hpp>
#include <detray/io/backend/homogeneous_material_reader.hpp>
#include <detray/io/backend/material_map_reader.hpp>
#include <detray/io/backend/surface_grid_reader.hpp>
#include <detray/io/frontend/detector_writer.hpp>
#include <detray/io/frontend/detector_writer_config.hpp>
#include <detray/io/json/json.hpp>
#include <detray/utils/consistency_checker.hpp>
#include <vecmem/memory/host_memory_resource.hpp>

namespace ActsPlugins::DetrayGeometryConverter {

/// @brief conversion method from ACTS TrackingGeometry to detray detector
/// @tparam metadata_t the detector metadata type
///
/// @param mr the memory resource to use for the detray detector construction
/// @param gctx the geometry context
/// @param trackingGeometry the ACTS tracking geometry to convert
/// @param beampipeVolumeName the beampipe volume name
/// @param detectorName the name to set for the detray detector (optional,
///     if not set, it will be taken from the payloads or defaulted to empty)
/// @param logLevel the logging level to use for the conversion process
/// @param convertMaterial whether to convert material information from ACTS to detray
/// @param convertSurfaceGrids whether to convert surface grid information from ACTS to detray
///
/// This method performs the following steps:
/// 1. It searches for the beampipe volume in the ACTS tracking geometry using
///    the provided beampipeVolumeName. If found, it sets this volume in the
///    DetrayPayloadConverter configuration. If not found, it logs a warning.
/// 2. It creates a DetrayPayloadConverter instance with the configuration and
///    converts the ACTS tracking geometry into Detray payloads.
/// 3. It builds a detray detector from the converted payloads using the
/// detray::detector_builder.
///
/// @return A pair of the built detray detector and its volume name map.
template <typename metadata_t>
std::pair<std::shared_ptr<detray::detector<metadata_t>>, detray::name_map>
toDetray(vecmem::memory_resource& mr, const Acts::GeometryContext& gctx,
         const Acts::TrackingGeometry& trackingGeometry,
         const std::string& beampipeVolumeName,
         const std::string& detectorName = "",
         Acts::Logging::Level logLevel = Acts::Logging::INFO,
         bool convertMaterial = true, bool convertSurfaceGrids = true) {
  auto localLogger =
      Acts::getDefaultLogger("DetrayGeometryConverter", logLevel);
  auto payloadLogger = localLogger->clone("DetrayPayloadConverter");
  ACTS_LOCAL_LOGGER(std::move(localLogger));

  // ── Convert TrackingGeometry → detray payloads ──────────────────────────
  DetrayPayloadConverter::Config convCfg;

  ACTS_INFO("Looking for beampipe volume: " << beampipeVolumeName);

  // Find beampipe volume
  trackingGeometry.apply(
      [&beampipeVolumeName, &convCfg](const Acts::TrackingVolume& volume) {
        if (volume.volumeName() == beampipeVolumeName) {
          convCfg.beampipeVolume = &volume;
        }
      });

  if (convCfg.beampipeVolume == nullptr) {
    ACTS_WARNING("DetrayGeometryProvider: beampipe volume '"
                 << beampipeVolumeName << "' not found");
  }

  DetrayPayloadConverter converter(convCfg, std::move(payloadLogger));

  auto payloads = converter.convertTrackingGeometry(gctx, trackingGeometry);

  // ── Build detray detector from payloads ──────────────────────────────────
  using detector_t = detray::detector<metadata_t>;

  detray::detector_builder<metadata_t> detectorBuilder{};

  detray::io::geometry_reader::from_payload<detector_t>(detectorBuilder,
                                                        *payloads.detector);

  if (convertMaterial) {
    detray::io::homogeneous_material_reader::from_payload<detector_t>(
        detectorBuilder, *payloads.homogeneousMaterial);

    detray::io::material_map_reader<std::integral_constant<std::size_t, 2>>::
        from_payload<detector_t>(detectorBuilder,
                                 std::move(*payloads.materialGrids));
  }

  if (convertSurfaceGrids) {
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

  detray::name_map names{};
  auto det = std::make_shared<detector_t>(detectorBuilder.build(mr, names));

  return {std::move(det), std::move(names)};
}

/// Build a mapping from detray surface identifiers to ACTS
/// GeometryIdentifiers for all sensitive surfaces in the detray detector.
/// @tparam detector_t the type of the detray detector
/// @param detrayDetector the detray detector to build the mapping for
/// @param logLevel the logging level to use for the mapping process
///
/// @return an unordered map mapping detray surface identifiers to ACTS GeometryIdentifiers
template <typename detector_t>
std::unordered_map<detray::geometry::identifier, Acts::GeometryIdentifier>
buildDetrayToActsMap(const detector_t& detrayDetector,
                     Acts::Logging::Level logLevel = Acts::Logging::INFO) {
  auto localLogger =
      Acts::getDefaultLogger("DetrayGeometryConverter", logLevel);
  ACTS_LOCAL_LOGGER(std::move(localLogger));
  // ── Build detray→Acts geometry ID map
  // ───────────────────────────────────
  std::unordered_map<detray::geometry::identifier, Acts::GeometryIdentifier>
      detrayToActsMap;

  for (const auto& surface : detrayDetector.surfaces()) {
    // surface.source is the Acts GeometryIdentifier encoded as uint64
    const Acts::GeometryIdentifier actsId(surface.source);
    if (actsId.sensitive() == 0) {
      continue;  // skip portals and passives
    }
    detrayToActsMap[surface.identifier()] = actsId;
  }

  ACTS_INFO("DetrayGeometryProvider: built detray→Acts map with "
            << detrayToActsMap.size() << " sensitive surfaces");
  return detrayToActsMap;
}

}  // namespace ActsPlugins::DetrayGeometryConverter
