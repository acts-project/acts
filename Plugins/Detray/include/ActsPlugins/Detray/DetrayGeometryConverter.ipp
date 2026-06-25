// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Implementation of the DetrayGeometryConverter::convert member template.
//
// Include this header (instead of DetrayGeometryConverter.hpp) when you need to
// instantiate `convert` for a metadata type that is not part of the closed set
// declared in DetrayMetadata.hpp.

#include "ActsPlugins/Detray/DetrayGeometryConverter.hpp"

#include <utility>

#include <detray/builders/detector_builder.hpp>
#include <detray/io/backend/geometry_reader.hpp>
#include <detray/io/backend/homogeneous_material_reader.hpp>
#include <detray/io/backend/material_map_reader.hpp>
#include <detray/io/backend/surface_grid_reader.hpp>

namespace ActsPlugins {

template <typename metadata_t>
DetrayGeometryConverter::DetrayGeometry<metadata_t>
DetrayGeometryConverter::convert(
    vecmem::memory_resource& mr, const Acts::GeometryContext& gctx,
    const std::shared_ptr<const Acts::TrackingGeometry>& trackingGeometry,
    const std::string& detectorName) const {
  using detector_t = detray::detector<metadata_t>;

  if (trackingGeometry == nullptr) {
    throw std::invalid_argument(
        "DetrayGeometryConverter: trackingGeometry must not be null");
  }

  // ── Convert TrackingGeometry → detray payloads ──────────────────────────
  auto payloads =
      m_cfg.payloadConverter->convertTrackingGeometry(gctx, *trackingGeometry);

  // ── Build detray detector from payloads ─────────────────────────────────
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

  return result;
}

}  // namespace ActsPlugins
