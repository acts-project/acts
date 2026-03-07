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

  /// Generic lookup from object pointer identity to serialized object ID.
  template <typename object_t, const char* kContext>
  struct PointerToIdLookup {
    /// Insert a new object to ID mapping.
    ///
    /// @param object is the source object pointer key
    /// @param objectId is the serialized ID to assign
    ///
    /// @return true if insertion happened, false if the object was already
    ///         present
    bool emplace(const object_t& object, std::size_t objectId) {
      return m_objectIds.emplace(&object, objectId).second;
    }

    /// Resolve a serialized object ID from an object reference.
    ///
    /// @param object is the source object key
    ///
    /// @return associated serialized object ID
    ///
    /// @throw std::invalid_argument if the object is not in the lookup
    std::size_t at(const object_t& object) const {
      auto it = m_objectIds.find(&object);
      if (it == m_objectIds.end()) {
        throw std::invalid_argument(
            "Pointer-to-ID lookup failed for " + std::string{kContext} +
            ": object is outside serialized hierarchy");
      }
      return it->second;
    }

   private:
    std::unordered_map<const object_t*, std::size_t> m_objectIds;
  };

  /// Generic lookup from serialized ID to pointer-like object holder.
  ///
  /// `pointer_t` can be a raw pointer (`object_t*`) or an owning pointer-like
  /// type such as `std::shared_ptr<object_t>`.
  template <typename object_t, typename pointer_t, const char* kContext>
  struct IdToPointerLikeLookup {
    /// Insert a new serialized ID to object mapping.
    ///
    /// @param objectId is the serialized ID key
    /// @param object is the target pointer-like object
    ///
    /// @return true if insertion happened, false if the ID was already present
    bool emplace(std::size_t objectId, pointer_t object) {
      return m_objects.emplace(objectId, std::move(object)).second;
    }

    /// Try to find a mapped pointer-like object by serialized ID.
    ///
    /// @param objectId is the serialized ID key
    ///
    /// @return mapped pointer-like object, or null-equivalent if not found
    pointer_t find(std::size_t objectId) const {
      auto it = m_objects.find(objectId);
      return it == m_objects.end() ? pointer_t{} : it->second;
    }

    /// Resolve a mapped pointer-like object by serialized ID.
    ///
    /// @param objectId is the serialized ID key
    ///
    /// @return mapped pointer-like object reference
    ///
    /// @throw std::invalid_argument if the ID is not mapped
    const pointer_t& at(std::size_t objectId) const {
      auto it = m_objects.find(objectId);
      if (it == m_objects.end()) {
        throw std::invalid_argument(
            "ID-to-pointer lookup failed for " + std::string{kContext} +
            ": unknown serialized object ID");
      }
      return it->second;
    }

   private:
    std::unordered_map<std::size_t, pointer_t> m_objects;
  };

  static inline constexpr char kVolumeLookupContext[] = "volume";
  static inline constexpr char kPortalLookupContext[] = "portal";

  using VolumeIdLookup =
      PointerToIdLookup<TrackingVolume, kVolumeLookupContext>;
  using VolumePointerLookup = IdToPointerLikeLookup<
      TrackingVolume, TrackingVolume*, kVolumeLookupContext>;
  using PortalIdLookup = PointerToIdLookup<Portal, kPortalLookupContext>;
  using PortalPointerLookup = IdToPointerLikeLookup<
      Portal, std::shared_ptr<Portal>, kPortalLookupContext>;

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
