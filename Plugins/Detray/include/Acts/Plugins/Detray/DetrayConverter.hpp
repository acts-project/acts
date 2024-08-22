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
  DetrayDetector convert(
      const GeometryContext& gctx, const Detector& detector,
      vecmem::memory_resource& mr,
      [[maybe_unused]] const DetrayConversionUtils::Options& options = {});

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
