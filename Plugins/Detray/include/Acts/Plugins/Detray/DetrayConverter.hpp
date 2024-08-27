// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Detector/Detector.hpp"
#include "Acts/Plugins/Detray/DetrayConversionUtils.hpp"
#include "Acts/Plugins/Detray/DetrayGeometryConverter.hpp"
#include "Acts/Plugins/Detray/DetrayMaterialConverter.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <memory>

#include <detray/io/common/geometry_reader.hpp>
#include <detray/io/common/material_map_reader.hpp>
#include <detray/io/frontend/detector_writer_config.hpp>

namespace Acts {

using namespace Experimental;

class DetrayConverter {
 public:
  /// Constructor with logger
  DetrayConverter(std::unique_ptr<const Logger> logger =
                      getDefaultLogger("DetrayConverter", Logging::INFO));

  /// Convert an Acts::Experimental::Detector to a detray::detector object
  ///
  /// @param gctx the geometry context
  /// @param detector the detector to be converted
  /// @param mr the memory resource to be used
  /// @param options the conversion options
  ///
  /// @returns a detector of requested return type
  template <typename detector_t = DetrayDetector>
  detector_t convert(
      const GeometryContext& gctx, const Detector& detector,
      vecmem::memory_resource& mr,
      [[maybe_unused]] const DetrayConversionUtils::Options& options = {}) {
    // The building cache object
    DetrayConversionUtils::GeometryIdCache geoIdCache;

    typename detector_t::name_map names = {{0u, detector.name()}};

    // build detector
    detray::detector_builder<typename detector_t::metadata> detectorBuilder{};
    // (1) geometry
    detray::io::detector_payload detectorPayload =
        DetrayGeometryConverter::convertDetector(geoIdCache, gctx, detector,
                                                 logger());
    detray::io::geometry_reader::convert<detector_t>(detectorBuilder, names,
                                                     detectorPayload);
    // (2) material
    if constexpr (detray::detail::has_material_grids_v<detector_t>) {
      if (options.convertMaterial) {
        detray::io::detector_grids_payload<detray::io::material_slab_payload,
                                           detray::io::material_id>
            materialPayload =
                DetrayMaterialConverter::convertSurfaceMaterialGrids(
                    geoIdCache, detector, logger());
        detray::io::material_map_reader<>::convert<detector_t>(
            detectorBuilder, names, materialPayload);
      }
    }

    detector_t detrayDetector(detectorBuilder.build(mr));

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

  /// Write the detector to json output
  ///
  /// @param dDetector is the detray detector (converted)
  /// @param names a name map for the detector volumes
  /// @param writer_cfg the writer configuration
  static void writeToJson(const DetrayDetector& dDetector,
                          const typename DetrayDetector::name_map& names = {},
                          detray::io::detector_writer_config writer_cfg = {});

 private:
  /// The logger instance
  std::unique_ptr<const Logger> m_logger = nullptr;

  // Return the logging instance
  const Acts::Logger& logger() const { return *m_logger; }
};

}  // namespace Acts
