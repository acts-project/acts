// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Detray/DetrayConverter.hpp"

Acts::DetrayConverter::DetrayConverter(
    std::unique_ptr<const Acts::Logger> logger)
    : m_logger(std::move(logger)) {}

void Acts::DetrayConverter::writeToJson(
    const DetrayDetector& dDetector,
    const typename DetrayDetector::name_map& names,
    detray::io::detector_writer_config writer_cfg) {
  writer_cfg.format(detray::io::format::json);
  detray::io::write_detector(dDetector, names, writer_cfg);
}

Acts::DetrayDetector Acts::DetrayConverter::convert(
    const GeometryContext& gctx, const Detector& detector,
    vecmem::memory_resource& mr,
    [[maybe_unused]] const DetrayConversionUtils::Options& options) {
  // The building cache object
  DetrayConversionUtils::GeometryIdCache geoIdCache;

  DetrayDetector::name_map names = {{0u, detector.name()}};

  // build detector
  detray::detector_builder<DetrayDetector::metadata> detectorBuilder{};
  // (1) geometry
  detray::io::detector_payload detectorPayload =
      DetrayGeometryConverter::convertDetector(geoIdCache, gctx, detector,
                                               logger());
  detray::io::geometry_reader::convert<DetrayDetector>(detectorBuilder, names,
                                                       detectorPayload);
  // (2) material
  if constexpr (detray::detail::has_material_grids_v<DetrayDetector>) {
    if (options.convertMaterial) {
      detray::io::detector_grids_payload<detray::io::material_slab_payload,
                                         detray::io::material_id>
          materialPayload =
              DetrayMaterialConverter::convertSurfaceMaterialGrids(
                  geoIdCache, detector, logger());
      detray::io::material_map_reader<>::convert<DetrayDetector>(
          detectorBuilder, names, std::move(materialPayload));
    }
  }

  DetrayDetector detrayDetector(detectorBuilder.build(mr));

  // Checks and print
  detray::detail::check_consistency(detrayDetector);

  // If configured, write the detector to json
  if (options.writeToJson) {
    // Create a writer configuration and write it out
    detray::io::detector_writer_config writerConfig{};
    writerConfig.m_write_material = options.convertMaterial;
    writerConfig.m_write_grids = options.convertSurfaceGrids;
    writeToJson(detrayDetector, names, writerConfig);
  }

  return detrayDetector;
}
