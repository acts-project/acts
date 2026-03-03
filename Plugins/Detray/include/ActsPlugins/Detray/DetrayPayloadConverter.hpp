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
/// @brief Converter for ACTS geometry to Detray payload format
class DetrayPayloadConverter {
 public:
  /// Detray surface material payload type
  using DetraySurfaceMaterial =
      std::variant<detray::io::grid_payload<detray::io::material_slab_payload,
                                            detray::io::material_id>,
                   detray::io::material_slab_payload>;

  /// Detray surface grid payload type
  using DetraySurfaceGrid =
      detray::io::grid_payload<std::size_t, detray::io::accel_id>;

  /// Function type for looking up surface indices in detray conversion
  using SurfaceLookupFunction =
      std::function<std::size_t(const Acts::Surface*)>;

  /// Convert homogeneous surface material
  /// @param material Homogeneous surface material
  /// @return Detray surface material payload
  static std::optional<DetraySurfaceMaterial> convertHomogeneousSurfaceMaterial(
      const Acts::HomogeneousSurfaceMaterial& material);

  /// Convert grid surface material
  /// @param material Grid surface material
  /// @return Detray surface material payload
  static std::optional<DetraySurfaceMaterial> convertGridSurfaceMaterial(
      const Acts::IGridSurfaceMaterialBase& material);

  /// Convert binned surface material
  /// @param material Binned surface material
  /// @return Detray surface material payload
  static std::optional<DetraySurfaceMaterial> convertBinnedSurfaceMaterial(
      const Acts::BinnedSurfaceMaterial& material);

  /// Convert proto surface material with bin utility
  /// @param material Proto surface material
  /// @return Detray surface material payload
  static std::optional<DetraySurfaceMaterial>
  convertProtoSurfaceMaterialBinUtility(
      const Acts::ProtoSurfaceMaterialT<Acts::BinUtility>& material);

  /// Convert proto surface material with proto axes
  /// @param material Proto surface material
  /// @return Detray surface material payload
  static std::optional<DetraySurfaceMaterial>
  convertProtoSurfaceMaterialProtoAxes(
      const Acts::ProtoSurfaceMaterialT<std::vector<Acts::DirectedProtoAxis>>&
          material);

  /// Convert surface array navigation policy
  /// @param policy Surface array navigation policy
  /// @param gctx Geometry context
  /// @param surfaceLookup Surface lookup function
  /// @param logger Logger instance
  /// @return Detray surface grid payload
  static std::optional<DetraySurfaceGrid> convertSurfaceArray(
      const Acts::SurfaceArrayNavigationPolicy& policy,
      const Acts::GeometryContext& gctx,
      const SurfaceLookupFunction& surfaceLookup, const Acts::Logger& logger);

  /// Convert try all navigation policy
  /// @param policy Try all navigation policy
  /// @param gctx Geometry context
  /// @param surfaceLookup Surface lookup function
  /// @param logger Logger instance
  /// @return Detray surface grid payload
  static std::optional<DetraySurfaceGrid> convertTryAllNavigationPolicy(
      const Acts::TryAllNavigationPolicy& policy,
      const Acts::GeometryContext& gctx,
      const SurfaceLookupFunction& surfaceLookup, const Acts::Logger& logger);

  /// Convert cylinder navigation policy
  /// @param policy Cylinder navigation policy
  /// @param gctx Geometry context
  /// @param surfaceLookup Surface lookup function
  /// @param logger Logger instance
  /// @return Detray surface grid payload
  static std::optional<DetraySurfaceGrid> convertCylinderNavigationPolicy(
      const Acts::CylinderNavigationPolicy& policy,
      const Acts::GeometryContext& gctx,
      const SurfaceLookupFunction& surfaceLookup, const Acts::Logger& logger);

  /// Convert multi layer navigation policy
  /// @param policy Multi layer navigation policy
  /// @param gctx Geometry context
  /// @param surfaceLookup Surface lookup function
  /// @param logger Logger instance
  /// @return Detray surface grid payload
  static std::optional<DetraySurfaceGrid> convertMultiLayerNavigationPolicy(
      const Acts::Experimental::MultiLayerNavigationPolicy& policy,
      const Acts::GeometryContext& gctx,
      const SurfaceLookupFunction& surfaceLookup, const Acts::Logger& logger);

  /// Convert multi navigation policy
  /// @param policy Multi navigation policy
  /// @param gctx Geometry context
  /// @param surfaceLookup Surface lookup function
  /// @param logger Logger instance
  /// @return Detray surface grid payload
  /// @note This is a noop, the payload converter will actually traverse the children via `visit`.
  static std::optional<DetraySurfaceGrid> convertMultiNavigationPolicy(
      const Acts::MultiNavigationPolicy& policy,
      const Acts::GeometryContext& gctx,
      const SurfaceLookupFunction& surfaceLookup, const Acts::Logger& logger);

  /// @brief Configuration for the Detray payload converter
  struct Config {
    Config() = default;

    /// Copy constructor
    Config(const Config&) = default;
    /// Move constructor
    Config(Config&&) = default;
    /// Copy assignment
    /// @return Reference to this object
    Config& operator=(const Config&) = default;
    /// Move assignment
    /// @return Reference to this object
    Config& operator=(Config&&) = default;
    /// Strategy for determining sensitive surfaces
    enum class SensitiveStrategy {
      /// Checks if the sensitive component of the surface is set to check if
      /// it's a sensitive surface
      Identifier,
      /// Check if the surface is a sensitive surface by checking for an
      /// associated detector element
      DetectorElement
    };
    /// Strategy to use for sensitive surface detection
    SensitiveStrategy sensitiveStrategy = SensitiveStrategy::Identifier;

    /// Detray MUST have beampipe volume at index 0
    const Acts::TrackingVolume* beampipeVolume = nullptr;

    /// Type dispatcher for converting navigation policies
    Acts::TypeDispatcher<Acts::INavigationPolicy,
                         std::optional<DetraySurfaceGrid>(
                             const Acts::GeometryContext& gctx,
                             const SurfaceLookupFunction& surfaceLookup,
                             const Acts::Logger& logger)>
        convertNavigationPolicy{
            convertSurfaceArray, convertTryAllNavigationPolicy,
            convertCylinderNavigationPolicy, convertMultiLayerNavigationPolicy,
            convertMultiNavigationPolicy};

    /// Type dispatcher for converting surface materials
    Acts::TypeDispatcher<Acts::ISurfaceMaterial,
                         std::optional<DetraySurfaceMaterial>()>
        convertSurfaceMaterial{
            convertHomogeneousSurfaceMaterial, convertBinnedSurfaceMaterial,
            convertGridSurfaceMaterial, convertProtoSurfaceMaterialProtoAxes,
            convertProtoSurfaceMaterialBinUtility};
  };

  /// Convert surface bounds to detray mask payload
  /// @param bounds the surface bounds to convert
  /// @param forPortal detray special cases the local parametrization for portals for performance reasons
  /// @return Detray mask payload
  static detray::io::mask_payload convertMask(const Acts::SurfaceBounds& bounds,
                                              bool forPortal);

  /// Convert surface
  /// @param gctx Geometry context
  /// @param surface Surface to convert
  /// @param portal Is portal surface
  /// @return Detray surface payload
  detray::io::surface_payload convertSurface(const Acts::GeometryContext& gctx,
                                             const Acts::Surface& surface,
                                             bool portal = false) const;

  /// Convert volume
  /// @param gctx Geometry context
  /// @param volume Volume to convert
  /// @return Detray volume payload
  detray::io::volume_payload convertVolume(
      const Acts::GeometryContext& gctx,
      const Acts::TrackingVolume& volume) const;

  /// @brief Container for all Detray payload outputs
  struct Payloads {
    /// Detector payload
    std::unique_ptr<detray::io::detector_payload> detector;

    /// Homogeneous material payload
    std::unique_ptr<detray::io::detector_homogeneous_material_payload>
        homogeneousMaterial;

    /// Material grids payload
    std::unique_ptr<detray::io::detector_grids_payload<
        detray::io::material_slab_payload, detray::io::material_id>>
        materialGrids;

    /// Surface grids payload
    std::unique_ptr<
        detray::io::detector_grids_payload<std::size_t, detray::io::accel_id>>
        surfaceGrids;

    /// Volume and surface names
    std::map<detray::dindex, std::string> names;
  };

  /// Convert tracking geometry
  /// @param gctx Geometry context
  /// @param geometry Tracking geometry to convert
  /// @return Detray payloads
  Payloads convertTrackingGeometry(
      const Acts::GeometryContext& gctx,
      const Acts::TrackingGeometry& geometry) const;

  /// Constructor
  /// @param config Configuration object
  /// @param logger Logger instance
  explicit DetrayPayloadConverter(
      const Config& config,
      std::unique_ptr<const Acts::Logger> logger = Acts::getDefaultLogger(
          "DetrayPayloadConverter", Acts::Logging::INFO));

  /// Convert material
  /// @param volume Tracking volume
  /// @param surfaceIndices Surface indices map
  /// @param volPayload Volume payload
  /// @return Material grids and volume material payloads
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
/// @brief Exception thrown when encountering unsupported material types in Detray conversion
class DetrayUnsupportedMaterialException : public std::runtime_error {
 public:
  /// Constructor
  /// @param name Name of the unsupported material type
  explicit DetrayUnsupportedMaterialException(std::string_view name);
};

}  // namespace ActsPlugins
