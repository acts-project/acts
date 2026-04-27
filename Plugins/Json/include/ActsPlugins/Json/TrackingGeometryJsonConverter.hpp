// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Navigation/INavigationPolicy.hpp"
#include "Acts/Navigation/SurfaceArrayNavigationPolicy.hpp"
#include "Acts/Surfaces/RegularSurface.hpp"
#include "Acts/Utilities/TypeDispatcher.hpp"
#include "ActsPlugins/Json/JsonKindDispatcher.hpp"

#include <memory>
#include <string>

#include <nlohmann/json.hpp>

namespace Acts {

class PortalLinkBase;
class Portal;
class Surface;
class TrackingGeometry;
class TrackingVolume;
class VolumeBounds;

/// @addtogroup json_plugin
/// @{

/// @brief Converter for tracking geometry JSON payloads
///
/// High-level conversion overview:
/// - Serialization:
///   - traverse the `TrackingVolume::volumes()` tree in depth-first order
///   - collect unique instances of surfaces, portals, volumes and assign stable
///   in-file IDs
///   - serialize the instances into their independent top-level tables
///   - encode object-to-object relationships through the assigned IDs
/// - Deserialization:
///   - validate schema header and collect all volume records
///   - instantiate volumes, portals, surfaces and build ID->pointer lookup
///   - attach child volumes to reconstruct the tree
///   - attach surfaces to portals and volumes via ID lookup
///   - portals to volumes via ID lookup
///   - return deserialized geometry
class TrackingGeometryJsonConverter {
 public:
  /// JSON serialization options for tracking geometry conversion.
  struct Options {};

  /// Generic lookup from object pointer identity to serialized object ID.
  template <typename object_t, const char* kContext>
  struct PointerToIdLookup;

  /// Generic lookup from serialized ID to pointer-like object holder.
  ///
  /// `pointer_t` can be a raw pointer (`object_t*`) or an owning pointer-like
  /// type such as `std::shared_ptr<object_t>`.
  template <typename object_t, typename pointer_t, const char* kContext>
  struct IdToPointerLikeLookup;

  /// Exception context for surfaces
  static inline constexpr char kSurfaceLookupContext[] = "surface";
  /// Exception context for portals
  static inline constexpr char kPortalLookupContext[] = "portal";
  /// Exception context for volumes
  static inline constexpr char kVolumeLookupContext[] = "volume";

  /// Surface map to its JSON ID
  using SurfaceIdLookup = PointerToIdLookup<Surface, kSurfaceLookupContext>;
  /// JSON ID map to its surface
  using SurfacePointerLookup =
      IdToPointerLikeLookup<RegularSurface, std::shared_ptr<RegularSurface>,
                            kSurfaceLookupContext>;

  /// Portal map to its JSON ID
  using PortalIdLookup = PointerToIdLookup<Portal, kPortalLookupContext>;
  /// JSON ID map to its portal
  using PortalPointerLookup =
      IdToPointerLikeLookup<Portal, std::shared_ptr<Portal>,
                            kPortalLookupContext>;

  /// Tracking volume map to its JSON ID
  using VolumeIdLookup =
      PointerToIdLookup<TrackingVolume, kVolumeLookupContext>;
  /// JSON ID map to its tracking volume
  using VolumePointerLookup =
      IdToPointerLikeLookup<TrackingVolume, TrackingVolume*,
                            kVolumeLookupContext>;

  /// Portal link encoder
  using PortalLinkEncoder =
      TypeDispatcher<PortalLinkBase,
                     nlohmann::json(const GeometryContext&,
                                    const TrackingGeometryJsonConverter&,
                                    const SurfaceIdLookup&,
                                    const VolumeIdLookup&)>;
  /// Portal link decoder
  using PortalLinkDecoder = JsonKindDispatcher<
      std::unique_ptr<PortalLinkBase>, const TrackingGeometryJsonConverter&,
      const SurfacePointerLookup&, const VolumePointerLookup&>;

  /// Volume bounds encoder
  using VolumeBoundsEncoder = TypeDispatcher<VolumeBounds, nlohmann::json()>;
  /// Volume bounds decoder
  using VolumeBoundsDecoder = JsonKindDispatcher<std::unique_ptr<VolumeBounds>>;

  /// Navigation policy encoder
  using NavigationPolicyEncoder =
      TypeDispatcher<INavigationPolicy,
                     nlohmann::json(
                         const Acts::TrackingGeometryJsonConverter&)>;
  /// Navigation policy decoder
  using NavigationPolicyDecoder =
      JsonKindDispatcher<std::unique_ptr<INavigationPolicy>,
                         const GeometryContext&,
                         const TrackingGeometryJsonConverter&,
                         const TrackingVolume&, const Acts::Logger&>;

  /// Configuration for the tracking geometry JSON converter.
  struct Config {
    /// Dispatcher for portal link serialization.
    PortalLinkEncoder encodePortalLink{};
    /// Dispatcher for volume bounds serialization.
    VolumeBoundsEncoder encodeVolumeBounds{};
    /// Dispatcher for navigation policy serialization.
    NavigationPolicyEncoder encodeNavigationPolicy{};
    /// Decoder dispatcher for portal links by kind tag.
    PortalLinkDecoder decodePortalLink{"kind", "portal link"};
    /// Decoder dispatcher for volume bounds by kind tag.
    VolumeBoundsDecoder decodeVolumeBounds{"kind", "volume bounds"};
    /// Decoder dispatcher for portal links by kind tag.
    NavigationPolicyDecoder decodeNavigationPolicy{"kind", "navigation policy"};

    /// Construct default config with all supported converters registered.
    ///
    /// @return configuration instance
    static Config defaultConfig();
  };

  /// @brief Construct converter with custom or default dispatch configuration.
  ///
  /// @param config is the conversion dispatch configuration
  explicit TrackingGeometryJsonConverter(
      Config config = Config::defaultConfig());

  /// @brief Convert a tracking geometry to JSON.
  ///
  /// @param gctx geometry context
  /// @param geometry tracking geometry to convert
  /// @param options options for the conversion
  ///
  /// @return serialized tracking geometry
  nlohmann::json toJson(const GeometryContext& gctx,
                        const TrackingGeometry& geometry,
                        const Options& options = Options{}) const;

  /// @brief Reconstruct a tracking geometry from JSON.
  ///
  /// @param gctx geometry context
  /// @param encoded serialized tracking geometry
  /// @param options options for the conversion
  ///
  /// @return pointer to deserialized geometry
  std::shared_ptr<TrackingGeometry> fromJson(
      const GeometryContext& gctx, const nlohmann::json& encoded,
      const Options& options = Options{}) const;

  /// @brief Convert a tracking volume hierarchy to JSON.
  ///
  /// @param gctx geometry context
  /// @param world top tracking volume in the hierarchy
  /// @param options options for the conversion
  ///
  /// @return serialized tracking volume hierarchy
  ///
  /// @note the geometry context is applied to the transformations
  /// during the serialization
  nlohmann::json trackingVolumeToJson(const GeometryContext& gctx,
                                      const TrackingVolume& world,
                                      const Options& options = Options{}) const;

  /// @brief Reconstruct a tracking volume hierarchy from JSON.
  ///
  /// @param gctx geometry context
  /// @param encoded serialized tracking volume hierarchy
  /// @param options options for the conversion
  ///
  /// @return pointer to deserialized tracking volume hierarchy
  ///
  /// @note currently the geometry context is only propagated to the
  /// Portal construction and the NavigationPolicy assignment
  std::shared_ptr<TrackingVolume> trackingVolumeFromJson(
      const GeometryContext& gctx, const nlohmann::json& encoded,
      const Options& options = Options{}) const;

  /// @brief Serialize one portal link using the configured dispatcher.
  ///
  /// @param gctx geometry context
  /// @param link portal link to serialize
  /// @param surfaceIds surface-to-id map for internal lookup
  /// @param volumeIds volume-to-id map for internal lookup
  ///
  /// @return serialized portal link
  nlohmann::json portalLinkToJson(const GeometryContext& gctx,
                                  const PortalLinkBase& link,
                                  const SurfaceIdLookup& surfaceIds,
                                  const VolumeIdLookup& volumeIds) const;

  /// @brief Deserialize one portal link using configured decoders.
  ///
  /// @param encoded serialized portal link
  /// @param surfaces id-to-surface map for internal lookup
  /// @param volumes id-to-volume map for internal lookup
  ///
  /// @return pointer to deserialized portal link
  std::unique_ptr<PortalLinkBase> portalLinkFromJson(
      const nlohmann::json& encoded, const SurfacePointerLookup& surfaces,
      const VolumePointerLookup& volumes) const;

  /// @brief Serialize navigation policy using the configured dispatcher.
  ///
  /// @param policy navigation policy to serialize
  ///
  /// @return serialized navigation policy
  nlohmann::json navigationPolicyToJson(const INavigationPolicy& policy) const;

  /// @brief Deserialize navigation policy using configured decoders.
  ///
  /// @param gctx geometry context
  /// @param encoded serialized navigation policy
  /// @param volume tracking volume to assign navigation policy to
  /// @param logger logging instance
  ///
  /// @return pointer to deserialized navigation policy
  std::unique_ptr<Acts::INavigationPolicy> navigationPolicyFromJson(
      const Acts::GeometryContext& gctx, const nlohmann::json& encoded,
      const Acts::TrackingVolume& volume, const Acts::Logger& logger) const;

 private:
  Config m_cfg;
};

NLOHMANN_JSON_SERIALIZE_ENUM(
    SurfaceArrayNavigationPolicy::LayerType,
    {{SurfaceArrayNavigationPolicy::LayerType::Cylinder, "Cylinder"},
     {SurfaceArrayNavigationPolicy::LayerType::Disc, "Disc"},
     {SurfaceArrayNavigationPolicy::LayerType::Plane, "Plane"}})

/// @}
}  // namespace Acts
