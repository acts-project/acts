// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Json/TrackingGeometryMaterialJsonWriter.hpp"

#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Material/ProtoSurfaceMaterial.hpp"
#include "Acts/Material/detail/MaterialSurfaceRegistry.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/BinAdjustment.hpp"

#include <algorithm>
#include <stdexcept>
#include <vector>

namespace ActsExamples {
TrackingGeometryMaterialJsonWriter::TrackingGeometryMaterialJsonWriter(
    const Config& config, Acts::Logging::Level /*level*/)
    : m_config(config) {
  if (m_config.filePath.empty()) {
    throw std::invalid_argument(
        "TrackingGeometryMaterialJsonWriter needs an output file path");
  }
}

void TrackingGeometryMaterialJsonWriter::writeMaterial(
    const Acts::TrackingGeometryMaterial& material) {
  Acts::TrackingGeometryMaterialJsonConverter{}.toFile(
      material, m_config.filePath, m_config.options);
}

void TrackingGeometryMaterialJsonWriter::write(
    const Acts::TrackingGeometry& geometry) {
  std::vector<const Acts::Surface*> surfaces;
  Acts::SurfaceMaterialMaps materials;
  geometry.visitSurfaces(
      [&](const Acts::Surface* surface) {
        auto payload = surface->surfaceMaterialSharedPtr();
        // Legacy DD4hep proto materials defer their ranges using zero-width
        // axes. Resolve these with the same helper as material mapping.
        if (const auto* proto =
                dynamic_cast<const Acts::ProtoSurfaceMaterial*>(payload.get());
            proto != nullptr &&
            std::ranges::any_of(
                proto->binning().binningData(),
                [](const auto& axis) { return axis.min == axis.max; })) {
          payload = std::make_shared<Acts::ProtoSurfaceMaterial>(
              Acts::adjustBinUtility(proto->binning(), *surface,
                                     m_config.geoContext),
              proto->mappingType(), proto->materialKey());
        }
        if (!payload && m_config.includeNonMaterial) {
          payload = std::make_shared<Acts::ProtoGridSurfaceMaterial>(
              Acts::MultiAxisSpec2D({Acts::AxisSpec::DeferredEquidistant(1),
                                     Acts::AxisSpec::DeferredEquidistant(1)}));
        }
        if (payload) {
          surfaces.push_back(surface);
          materials.emplace(surface->geometryId(), std::move(payload));
        }
      },
      false);
  auto material = Acts::detail::MaterialSurfaceRegistry(surfaces).materialMaps(
      std::move(materials));
  geometry.visitVolumes([](const Acts::TrackingVolume* volume) {
    if (volume->volumeMaterial() != nullptr) {
      throw std::invalid_argument(
          "TrackingGeometryMaterialJsonWriter supports surface material only; "
          "volume material "
          "cannot be exported");
    }
  });
  writeMaterial(material);
}
}  // namespace ActsExamples
