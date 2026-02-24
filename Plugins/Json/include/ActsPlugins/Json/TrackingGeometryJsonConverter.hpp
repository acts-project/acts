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
class Portal;
class TrackingGeometry;
class TrackingVolume;
class VolumeBounds;

/// @addtogroup json_plugin
/// @{

/// Converter for tracking geometry JSON payloads focused on volumes, volume
/// bounds and portals.
///
/// High-level conversion overview:
/// - Serialization:
///   - traverse the `TrackingVolume::volumes()` tree in depth-first order
///   - assign stable in-file volume IDs
///   - collect unique portals and assign stable in-file portal IDs
///   - write each volume transform, bounds payload, children IDs, and portal IDs
///   - write all unique portals once in a top-level portal table
///   - encode portal links by concrete kind via registered dispatchers
/// - Deserialization:
///   - validate schema header and collect all volume records
///   - instantiate all volumes first and build ID->pointer lookup
///   - attach child volumes to reconstruct the tree
///   - decode unique portals by kind, then attach shared portal pointers to
///     volumes via portal IDs
///   - return a reconstructed world `TrackingVolume` (or `TrackingGeometry`)
class TrackingGeometryJsonConverter {
 public:
  /// JSON serialization options for tracking geometry conversion.
  struct Options {};

  /// Lookup structure from volume pointer to serialized volume ID.
  struct VolumeIdLookup {
    /// Insert a new volume to ID mapping.
    ///
    /// @param volume is the source volume pointer key
    /// @param volumeId is the serialized ID to assign
    ///
    /// @return true if insertion happened, false if the volume was already
    ///         present
    bool emplace(const TrackingVolume& volume, std::size_t volumeId) {
      return m_volumeIds.emplace(&volume, volumeId).second;
    }

    /// Resolve a serialized volume ID from a volume reference.
    ///
    /// @param volume is the source volume key
    ///
    /// @return associated serialized volume ID
    ///
    /// @throw std::invalid_argument if the volume is not in the lookup
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
    /// Insert a new serialized ID to volume pointer mapping.
    ///
    /// @param volumeId is the serialized ID key
    /// @param volume is the target volume object
    ///
    /// @return true if insertion happened, false if the ID was already present
    bool emplace(std::size_t volumeId, TrackingVolume& volume) {
      return m_volumes.emplace(volumeId, &volume).second;
    }

    /// Try to find a mapped volume pointer by serialized ID.
    ///
    /// @param volumeId is the serialized ID key
    ///
    /// @return raw pointer to the mapped volume, or nullptr if not found
    TrackingVolume* find(std::size_t volumeId) const {
      auto it = m_volumes.find(volumeId);
      return it == m_volumes.end() ? nullptr : it->second;
    }

    /// Resolve a mapped volume reference by serialized ID.
    ///
    /// @param volumeId is the serialized ID key
    ///
    /// @return reference to mapped volume
    ///
    /// @throw std::invalid_argument if the ID is not mapped
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

  /// Lookup structure from portal pointer to serialized portal ID.
  struct PortalIdLookup {
    /// Insert a new portal to ID mapping.
    ///
    /// @param portal is the source portal pointer key
    /// @param portalId is the serialized ID to assign
    ///
    /// @return true if insertion happened, false if the portal was already
    ///         present
    bool emplace(const Portal& portal, std::size_t portalId) {
      return m_portalIds.emplace(&portal, portalId).second;
    }

    /// Resolve a serialized portal ID from a portal reference.
    ///
    /// @param portal is the source portal key
    ///
    /// @return associated serialized portal ID
    ///
    /// @throw std::invalid_argument if the portal is not in the lookup
    std::size_t at(const Portal& portal) const {
      auto it = m_portalIds.find(&portal);
      if (it == m_portalIds.end()) {
        throw std::invalid_argument(
            "Portal lookup failed: portal is outside serialized hierarchy");
      }
      return it->second;
    }

   private:
    std::unordered_map<const Portal*, std::size_t> m_portalIds;
  };

  /// Lookup structure from serialized portal ID to reconstructed shared portal
  /// pointer.
  struct PortalPointerLookup {
    /// Insert a new serialized ID to portal pointer mapping.
    ///
    /// @param portalId is the serialized ID key
    /// @param portal is the target shared portal object
    ///
    /// @return true if insertion happened, false if the ID was already present
    bool emplace(std::size_t portalId, std::shared_ptr<Portal> portal) {
      return m_portals.emplace(portalId, std::move(portal)).second;
    }

    /// Try to find a mapped portal pointer by serialized ID.
    ///
    /// @param portalId is the serialized ID key
    ///
    /// @return mapped shared portal pointer, or nullptr if not found
    std::shared_ptr<Portal> find(std::size_t portalId) const {
      auto it = m_portals.find(portalId);
      return it == m_portals.end() ? nullptr : it->second;
    }

    /// Resolve a mapped portal pointer by serialized ID.
    ///
    /// @param portalId is the serialized ID key
    ///
    /// @return mapped shared portal pointer
    ///
    /// @throw std::invalid_argument if the ID is not mapped
    const std::shared_ptr<Portal>& at(std::size_t portalId) const {
      auto it = m_portals.find(portalId);
      if (it == m_portals.end()) {
        throw std::invalid_argument(
            "Portal pointer lookup failed: unknown serialized portal ID");
      }
      return it->second;
    }

   private:
    std::unordered_map<std::size_t, std::shared_ptr<Portal>> m_portals;
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

  /// Construct converter with custom or default dispatch configuration.
  ///
  /// @param config is the conversion dispatch configuration
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
