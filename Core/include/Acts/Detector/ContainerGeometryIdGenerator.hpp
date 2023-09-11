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
#include <optional>

namespace Acts {

class Surface;

namespace Experimental {

class Portal;
class DetectorVolume;

/// @brief This is a geometry ID generator for a layered container volume
///
/// It can optionally reverse the labelling and override existing ids
/// It also allows to optionally set a volume id by configuration
class ContainerGeometryIdGenerator : public IGeometryIdGenerator {
 public:
  /// @brief  Nested config struct
  struct Config {
    /// Allow overriding of already existing sub ids
    bool overrideExistingIds = false;
    /// Reverse labelling
    bool reverse = false;
    /// Optional volume id
    std::optional<unsigned int> volumeId = std::nullopt;
  };

  /// @brief Constructor with config
  ///
  /// @param cfg is the geometry configuration object
  ContainerGeometryIdGenerator(const Config& cfg) : m_cfg(cfg) {}

  virtual ~ContainerGeometryIdGenerator() = default;

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
