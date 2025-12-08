// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Material/BinnedSurfaceMaterial.hpp"
#include "Acts/Material/GridSurfaceMaterial.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Material/ProtoSurfaceMaterial.hpp"
#include "Acts/Navigation/CylinderNavigationPolicy.hpp"
#include "Acts/Navigation/INavigationPolicy.hpp"
#include "Acts/Navigation/MultiLayerNavigationPolicy.hpp"
#include "Acts/Navigation/MultiNavigationPolicy.hpp"
#include "Acts/Navigation/SurfaceArrayNavigationPolicy.hpp"
#include "Acts/Navigation/TryAllNavigationPolicy.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/TypeDispatcher.hpp"

#include <map>
#include <memory>

/// @cond
namespace detray {

using dindex = unsigned int;

namespace io {

struct detector_payload;
struct transform_payload;
struct mask_payload;
struct surface_payload;
struct volume_payload;
struct material_slab_payload;
struct material_volume_payload;
struct detector_homogeneous_material_payload;

template <typename, typename>
struct grid_payload;
enum class material_id : unsigned int;

template <typename, typename>
struct detector_grids_payload;

enum class accel_id : unsigned int;

}  // namespace io
}  // namespace detray
/// @endcond

namespace Acts {
class GeometryContext;
class TrackingGeometry;
class SurfaceBounds;
class Surface;
class Portal;
class TrackingVolume;
class PortalLinkBase;
class MaterialSlab;
class ISurfaceMaterial;
}  // namespace Acts

namespace ActsPlugins {

/// @ingroup detray_plugin
class DetrayPayloadConverter {
 public:
  using DetraySurfaceMaterial =
      std::variant<detray::io::grid_payload<detray::io::material_slab_payload,
                                            detray::io::material_id>,
                   detray::io::material_slab_payload>;

  using DetraySurfaceGrid =
      detray::io::grid_payload<std::size_t, detray::io::accel_id>;

  /// Function type for looking up surface indices in detray conversion
  using SurfaceLookupFunction =
      std::function<std::size_t(const Acts::Surface*)>;

  static std::optional<DetraySurfaceMaterial> convertHomogeneousSurfaceMaterial(
      const Acts::HomogeneousSurfaceMaterial& material);

  static std::optional<DetraySurfaceMaterial> convertGridSurfaceMaterial(
      const Acts::IGridSurfaceMaterialBase& material);

  static std::optional<DetraySurfaceMaterial> convertBinnedSurfaceMaterial(
      const Acts::BinnedSurfaceMaterial& material);

  static std::optional<DetraySurfaceMaterial>
  convertProtoSurfaceMaterialBinUtility(
      const Acts::ProtoSurfaceMaterialT<Acts::BinUtility>& material);

  static std::optional<DetraySurfaceMaterial>
  convertProtoSurfaceMaterialProtoAxes(
      const Acts::ProtoSurfaceMaterialT<std::vector<Acts::DirectedProtoAxis>>&
          material);

  static std::optional<DetraySurfaceGrid> convertSurfaceArray(
      const Acts::SurfaceArrayNavigationPolicy& policy,
      const Acts::GeometryContext& gctx,
      const SurfaceLookupFunction& surfaceLookup, const Acts::Logger& logger);

  static std::optional<DetraySurfaceGrid> convertTryAllNavigationPolicy(
      const Acts::TryAllNavigationPolicy& policy,
      const Acts::GeometryContext& gctx,
      const SurfaceLookupFunction& surfaceLookup, const Acts::Logger& logger);

  static std::optional<DetraySurfaceGrid> convertCylinderNavigationPolicy(
      const Acts::CylinderNavigationPolicy& policy,
      const Acts::GeometryContext& gctx,
      const SurfaceLookupFunction& surfaceLookup, const Acts::Logger& logger);

  static std::optional<DetraySurfaceGrid> convertMultiLayerNavigationPolicy(
      const Acts::Experimental::MultiLayerNavigationPolicy& policy,
      const Acts::GeometryContext& gctx,
      const SurfaceLookupFunction& surfaceLookup, const Acts::Logger& logger);

  // This is a noop, the payload converter will actually traverse the children
  // via `visit`.
  static std::optional<DetraySurfaceGrid> convertMultiNavigationPolicy(
      const Acts::MultiNavigationPolicy& policy,
      const Acts::GeometryContext& gctx,
      const SurfaceLookupFunction& surfaceLookup, const Acts::Logger& logger);

  struct Config {
    Config() = default;

    Config(const Config&) = default;
    Config(Config&&) = default;
    Config& operator=(const Config&) = default;
    Config& operator=(Config&&) = default;
    enum class SensitiveStrategy {
      /// Checks if the sensitive component of the surface is set to check if
      /// it's a sensitive surface
      Identifier,
      /// Check if the surface is a sensitive surface by checking for an
      /// associated detector element
      DetectorElement
    };
    SensitiveStrategy sensitiveStrategy = SensitiveStrategy::Identifier;

    /// Detray MUST have beampipe volume at index 0
    const Acts::TrackingVolume* beampipeVolume = nullptr;

    Acts::TypeDispatcher<Acts::INavigationPolicy,
                         std::optional<DetraySurfaceGrid>(
                             const Acts::GeometryContext& gctx,
                             const SurfaceLookupFunction& surfaceLookup,
                             const Acts::Logger& logger)>
        convertNavigationPolicy{
            convertSurfaceArray, convertTryAllNavigationPolicy,
            convertCylinderNavigationPolicy, convertMultiLayerNavigationPolicy,
            convertMultiNavigationPolicy};

    Acts::TypeDispatcher<Acts::ISurfaceMaterial,
                         std::optional<DetraySurfaceMaterial>()>
        convertSurfaceMaterial{
            convertHomogeneousSurfaceMaterial, convertBinnedSurfaceMaterial,
            convertGridSurfaceMaterial, convertProtoSurfaceMaterialProtoAxes,
            convertProtoSurfaceMaterialBinUtility};
  };

  /// @param bounds the surface bounds to convert
  /// @param forPortal detray special cases the local parametrization for portals for performance reasons
  static detray::io::mask_payload convertMask(const Acts::SurfaceBounds& bounds,
                                              bool forPortal);

  detray::io::surface_payload convertSurface(const Acts::GeometryContext& gctx,
                                             const Acts::Surface& surface,
                                             bool portal = false) const;

  detray::io::volume_payload convertVolume(
      const Acts::TrackingVolume& volume) const;

  struct Payloads {
    // Unique pointers used to be able to forward declare the type
    std::unique_ptr<detray::io::detector_payload> detector;

    std::unique_ptr<detray::io::detector_homogeneous_material_payload>
        homogeneousMaterial;

    std::unique_ptr<detray::io::detector_grids_payload<
        detray::io::material_slab_payload, detray::io::material_id>>
        materialGrids;

    std::unique_ptr<
        detray::io::detector_grids_payload<std::size_t, detray::io::accel_id>>
        surfaceGrids;

    std::map<detray::dindex, std::string> names;
  };

  Payloads convertTrackingGeometry(
      const Acts::GeometryContext& gctx,
      const Acts::TrackingGeometry& geometry) const;

  explicit DetrayPayloadConverter(
      const Config& config,
      std::unique_ptr<const Acts::Logger> logger = Acts::getDefaultLogger(
          "DetrayPayloadConverter", Acts::Logging::INFO));

  std::pair<std::vector<detray::io::grid_payload<
                detray::io::material_slab_payload, detray::io::material_id>>,
            detray::io::material_volume_payload>
  convertMaterial(const Acts::TrackingVolume& volume,

                  const std::unordered_map<const Acts::Surface*, std::size_t>&
                      surfaceIndices,
                  detray::io::volume_payload& volPayload) const;

 private:
  void handlePortalLink(
      const Acts::GeometryContext& gctx, const Acts::TrackingVolume& volume,
      detray::io::volume_payload& volPayload,
      const std ::function<std::size_t(const Acts::TrackingVolume*)>&
          volumeLookup,
      std::unordered_map<const Acts::Surface*, std::size_t>& surfaceIndices,
      const Acts::PortalLinkBase& link) const;

  void makeEndOfWorld(
      const Acts::GeometryContext& gctx, detray::io::volume_payload& volPayload,
      std::unordered_map<const Acts::Surface*, std::size_t>& surfaceIndices,
      const Acts::Surface& surface) const;

  void handlePortal(
      const Acts::GeometryContext& gctx, const Acts::TrackingVolume& volume,
      detray::io::volume_payload& volPayload,
      const std::function<std::size_t(const Acts::TrackingVolume*)>&
          volumeLookup,
      std::unordered_map<const Acts::Surface*, std::size_t>& surfaceIndices,
      const Acts::Portal& portal) const;

  Config m_cfg;

  const Acts::Logger& logger() const { return *m_logger; }
  std::unique_ptr<const Acts::Logger> m_logger;
};

/// @ingroup detray_plugin
class DetrayUnsupportedMaterialException : public std::runtime_error {
 public:
  explicit DetrayUnsupportedMaterialException(std::string_view name);
};

}  // namespace ActsPlugins
