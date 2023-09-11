// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Detector/interface/IGeometryIdGenerator.hpp"

#include <any>

namespace Acts {

class Surface;

namespace Experimental {

class Portal;
class DetectorVolume;

/// @brief This is a straight forward version of the geometry id generator
/// It assigns the geometry ids in a sequential order, loopoing over volumes
/// portals and surfaces.
///
class SequentialGeometryIdGenerator : public IGeometryIdGenerator {
 public:
  /// @brief  Nested config struct
  struct Config {
    bool overrideExistingIds = false;
  };

  /// @brief Nested cache struct
  struct Cache {
    unsigned int volumeCount = 0u;
    unsigned int portalCount = 0u;
    unsigned int passiveCount = 0u;
    unsigned int sensitiveCount = 0u;
  };

  /// @brief Constructor with config
  ///
  /// @param cfg is the geometry configuration object
  SequentialGeometryIdGenerator(const Config& cfg) : m_cfg(cfg) {}

  virtual ~SequentialGeometryIdGenerator() = default;

  /// @brief  Interface method to generata a geometry id cache
  /// @return a geometry id cache decorated in a std::any object
  IGeometryIdGenerator::GeoIdCache generateCache() const final override;

  /// @brief Method for assigning a geometry id to a detector volume
  ///
  /// @param cache is the cahce object for e.g. object counting
  /// @param dVolume the detector volume to assign the geometry id to
  void assignGeometryId(IGeometryIdGenerator::GeoIdCache& cache,
                        DetectorVolume& dVolume) const final override;

  /// @brief Method for assigning a geometry id to a portal
  ///
  /// @param cache is the cahce object for e.g. object counting
  /// @param portal the portal to assign the geometry id to
  void assignGeometryId(IGeometryIdGenerator::GeoIdCache& cache,
                        Portal& portal) const final override;

  /// @brief Method for assigning a geometry id to a surface
  ///
  /// @param cache is the cahce object for e.g. object counting
  /// @param surface the surface to assign the geometry id to
  void assignGeometryId(IGeometryIdGenerator::GeoIdCache& cache,
                        Surface& surface) const final override;

 private:
  Config m_cfg;
};

}  // namespace Experimental
}  // namespace Acts
