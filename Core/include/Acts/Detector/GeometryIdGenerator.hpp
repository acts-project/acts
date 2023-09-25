// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Detector/interface/IGeometryIdGenerator.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <any>

namespace Acts {

class Surface;

namespace Experimental {

class Portal;
class DetectorVolume;

/// @brief This is the default implementation of the geometry id generator
///
/// It is a simple counter based generator that assigns ids to the geometry
/// and increments the counter for each object type.
///
/// Sub counters, i.e. for the sensitive surfaces, are reset after the volume
/// call, such that a new volume or layer volume would start from counter 0
/// again.
///
/// If the generator is configured to override existing ids, it will do so
/// and not respect previously assigned ids.
///
/// If the generator is configured in container mode, it will increase the
/// layer id for each layer volume with confined surfaces.
///
class GeometryIdGenerator final : public IGeometryIdGenerator {
 public:
  /// @brief  Nested config struct
  struct Config {
    /// Container mode
    bool containerMode = false;
    /// Container id (if container mode), will not be incremented
    unsigned int containerId = 0u;
    /// Resetting mode
    bool resetSubCounters = true;
    /// Force override existing ids
    bool overrideExistingIds = false;
  };

  /// @brief Nested cache struct
  struct Cache {
    /// Cache count of the volume, for non-container mode
    unsigned int volumeCount = 0u;
    /// Cache count of the layer volume, for container mode
    unsigned int layerCount = 0u;
    /// Cache count of the portal surfaces
    unsigned int portalCount = 0u;
    /// Cache count of passive surfaces
    unsigned int passiveCount = 0u;
    /// Cache count of sensitive surfaces
    unsigned int sensitiveCount = 0u;
  };

  /// @brief Constructor with config
  ///
  /// @param cfg is the geometry configuration object
  /// @param mlogger is the logging instance
  GeometryIdGenerator(const Config& cfg,
                      std::unique_ptr<const Logger> mlogger = getDefaultLogger(
                          "GeometryIdGenerator", Logging::INFO))
      : m_cfg(cfg), m_logger(std::move(mlogger)) {}

  ~GeometryIdGenerator() override = default;

  /// @brief Interface method to generata a geometry id cache
  /// @return a geometry id cache decorated in a std::any object
  IGeometryIdGenerator::GeoIdCache generateCache() const final;

  /// @brief Method for assigning a geometry id to a detector volume
  ///
  /// @param cache is the cache object for e.g. object counting
  /// @param dVolume the detector volume to assign the geometry id to
  void assignGeometryId(IGeometryIdGenerator::GeoIdCache& cache,
                        DetectorVolume& dVolume) const final;

  /// @brief Method for assigning a geometry id to a portal
  ///
  /// @param cache is the cache object for e.g. object counting
  /// @param portal the portal to assign the geometry id to
  void assignGeometryId(IGeometryIdGenerator::GeoIdCache& cache,
                        Portal& portal) const final;

  /// @brief Method for assigning a geometry id to a surface
  ///
  /// @param cache is the cache object for e.g. object counting
  /// @param surface the surface to assign the geometry id to
  void assignGeometryId(IGeometryIdGenerator::GeoIdCache& cache,
                        Surface& surface) const final;

 private:
  /// @brief Helper method to get the volume id from the cache
  ///
  /// @param cache the provided cache
  /// @param incrementLayer if true, the layer counter is incremented
  ///
  /// @return a valid geometry identifier
  GeometryIdentifier volumeId(Cache& cache, bool incrementLayer = true) const;

  /// Configuration object
  Config m_cfg;

  /// Private access method to the logger
  const Logger& logger() const { return *m_logger; }

  /// logging instance
  std::unique_ptr<const Logger> m_logger;
};

}  // namespace Experimental
}  // namespace Acts
