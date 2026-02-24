// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Utilities/TypeDispatcher.hpp"
#include "ActsPlugins/Json/JsonKindDispatcher.hpp"

#include <cstddef>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>

#include <nlohmann/json.hpp>

namespace Acts {

class PortalLinkBase;
class TrackingGeometry;
class TrackingVolume;
class VolumeBounds;

/// @addtogroup json_plugin
/// @{

/// Converter for tracking geometry JSON payloads focused on volumes,
/// volume bounds and portals.
class TrackingGeometryJsonConverter {
 public:
  /// JSON serialization options for tracking geometry conversion.
  struct Options {};

  /// Lookup structure from volume pointer to serialized volume ID.
  struct VolumeIdLookup {
    bool emplace(const TrackingVolume& volume, std::size_t volumeId) {
      return m_volumeIds.emplace(&volume, volumeId).second;
    }

    std::size_t at(const TrackingVolume& volume) const {
      auto it = m_volumeIds.find(&volume);
      if (it == m_volumeIds.end()) {
        throw std::invalid_argument(
            "Volume lookup failed: volume is outside serialized hierarchy");
      }
      return it->second;
    }

   private:
    std::unordered_map<const TrackingVolume*, std::size_t> m_volumeIds;
  };

  /// Lookup structure from serialized volume ID to reconstructed volume
  /// pointer.
  struct VolumePointerLookup {
    bool emplace(std::size_t volumeId, TrackingVolume& volume) {
      return m_volumes.emplace(volumeId, &volume).second;
    }

    TrackingVolume* find(std::size_t volumeId) const {
      auto it = m_volumes.find(volumeId);
      return it == m_volumes.end() ? nullptr : it->second;
    }

    TrackingVolume& at(std::size_t volumeId) const {
      auto* volume = find(volumeId);
      if (volume == nullptr) {
        throw std::invalid_argument(
            "Volume pointer lookup failed: unknown serialized volume ID");
      }
      return *volume;
    }

   private:
    std::unordered_map<std::size_t, TrackingVolume*> m_volumes;
  };

  using VolumeBoundsEncoder =
      TypeDispatcher<VolumeBounds, nlohmann::json()>;

  using PortalLinkEncoder = TypeDispatcher<
      PortalLinkBase,
      nlohmann::json(const GeometryContext&,
                     const TrackingGeometryJsonConverter&,
                     const VolumeIdLookup&)>;

  using VolumeBoundsDecoder =
      JsonKindDispatcher<std::unique_ptr<VolumeBounds>>;

  using PortalLinkDecoder =
      JsonKindDispatcher<std::unique_ptr<PortalLinkBase>,
                         const GeometryContext&,
                         const TrackingGeometryJsonConverter&,
                         const VolumePointerLookup&>;

  /// Configuration for the tracking geometry JSON converter.
  struct Config {
    /// Dispatcher for volume bounds serialization.
    VolumeBoundsEncoder encodeVolumeBounds{};

    /// Dispatcher for portal link serialization.
    PortalLinkEncoder encodePortalLink{};

    /// Decoder dispatcher for volume bounds by kind tag.
    VolumeBoundsDecoder decodeVolumeBounds{"kind", "volume bounds"};

    /// Decoder dispatcher for portal links by kind tag.
    PortalLinkDecoder decodePortalLink{"kind", "portal link"};

    /// Construct default config with all supported converters registered.
    static Config defaultConfig();
  };

  explicit TrackingGeometryJsonConverter(Config config = Config::defaultConfig());

  /// Convert a tracking geometry to JSON.
  nlohmann::json toJson(const GeometryContext& gctx,
                        const TrackingGeometry& geometry,
                        const Options& options = Options{}) const;

  /// Convert a tracking volume hierarchy to JSON.
  nlohmann::json toJson(const GeometryContext& gctx,
                        const TrackingVolume& world,
                        const Options& options = Options{}) const;

  /// Reconstruct a tracking volume hierarchy from JSON.
  std::shared_ptr<TrackingVolume> trackingVolumeFromJson(
      const GeometryContext& gctx, const nlohmann::json& encoded,
      const Options& options = Options{}) const;

  /// Reconstruct a tracking geometry from JSON.
  std::shared_ptr<TrackingGeometry> trackingGeometryFromJson(
      const GeometryContext& gctx, const nlohmann::json& encoded,
      const Options& options = Options{}) const;

  /// Serialize one portal link using the configured dispatcher.
  nlohmann::json portalLinkToJson(const GeometryContext& gctx,
                                  const PortalLinkBase& link,
                                  const VolumeIdLookup& volumeIds) const;

  /// Deserialize one portal link using configured decoders.
  std::unique_ptr<PortalLinkBase> portalLinkFromJson(
      const GeometryContext& gctx, const nlohmann::json& encoded,
      const VolumePointerLookup& volumes) const;

 private:
  Config m_cfg;
};

/// @}
}  // namespace Acts
