// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Material/IMaterialDecorator.hpp"
#include "Acts/Material/ISurfaceMaterial.hpp"
#include "Acts/Material/IVolumeMaterial.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "ActsPlugins/Json/MaterialMapJsonConverter.hpp"

// Convenience shorthand

namespace Acts {

/// @addtogroup json_plugin
/// @{

/// @brief Material decorator from Json format
///
/// This reads in material maps for surfaces and volumes
/// from a json file
class JsonMaterialDecorator : public IMaterialDecorator {
 public:
  /// Constructor with configuration
  /// @param rConfig the configuration for the material map reader
  /// @param jFileName the json file name to read
  /// @param level the logging level
  JsonMaterialDecorator(const MaterialMapJsonConverter::Config& rConfig,
                        const std::string& jFileName,
                        Acts::Logging::Level level);

  /// Parsed, format-independent maps for applying to a completed geometry.
  /// @return Material maps read from the input file
  const TrackingGeometryMaterial& materialMaps() const {
    return m_materialMaps;
  }

  /// Decorate a surface
  ///
  /// @param surface the non-cost surface that is decorated
  void decorate(Surface& surface) const final;

  /// Decorate a TrackingVolume
  ///
  /// @param volume the non-cost volume that is decorated
  void decorate(TrackingVolume& volume) const final;

 private:
  MaterialMapJsonConverter::Config m_readerConfig;
  TrackingGeometryMaterial m_materialMaps;

  std::unique_ptr<const Logger> m_logger;

  const Logger& logger() const { return *m_logger; }
};

/// @}
}  // namespace Acts
