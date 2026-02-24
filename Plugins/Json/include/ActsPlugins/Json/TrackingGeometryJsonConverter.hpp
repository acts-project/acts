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

#include <memory>
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

  /// Mapping from volume pointer to serialized volume ID.
  using VolumeIdMap = std::unordered_map<const TrackingVolume*, std::size_t>;

  /// Mapping from serialized volume ID to reconstructed volume pointer.
  using VolumePointerMap = std::unordered_map<std::size_t, TrackingVolume*>;

  using VolumeBoundsEncoder =
      TypeDispatcher<VolumeBounds, nlohmann::json()>;

  using PortalLinkEncoder = TypeDispatcher<
      PortalLinkBase,
      nlohmann::json(const GeometryContext&,
                     const TrackingGeometryJsonConverter&, const VolumeIdMap&)>;

  using VolumeBoundsDecoder =
      JsonKindDispatcher<std::unique_ptr<VolumeBounds>>;

  using PortalLinkDecoder =
      JsonKindDispatcher<std::unique_ptr<PortalLinkBase>,
                         const GeometryContext&,
                         const TrackingGeometryJsonConverter&,
                         const VolumePointerMap&>;

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
                                  const VolumeIdMap& volumeIds) const;

  /// Deserialize one portal link using configured decoders.
  std::unique_ptr<PortalLinkBase> portalLinkFromJson(
      const GeometryContext& gctx, const nlohmann::json& encoded,
      const VolumePointerMap& volumes) const;

 private:
  Config m_cfg;
};

/// @}
}  // namespace Acts
