// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/DetrayFwd.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <memory>

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

class DetrayPayloadConverter {
 public:
  struct Config {
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
    const TrackingVolume* beampipeVolume = nullptr;
  };

  static detray::io::transform_payload convertTransform(
      const Transform3& transform);

  /// @param forPortal detray special cases the local parametrization for portals for performance reasons
  static detray::io::mask_payload convertMask(const Acts::SurfaceBounds& bounds,
                                              bool forPortal);

  detray::io::surface_payload convertSurface(const GeometryContext& gctx,
                                             const Surface& surface,
                                             bool portal = false) const;

  detray::io::volume_payload convertVolume(const TrackingVolume& volume) const;

  struct Payloads {
    // Unique pointers used to be able to forward declare the type
    std::unique_ptr<detray::io::detector_payload> detector;
    std::unique_ptr<detray::io::detector_homogeneous_material_payload>
        homogeneousMaterial;
  };

  Payloads convertTrackingGeometry(const GeometryContext& gctx,
                                   const TrackingGeometry& geometry) const;

  explicit DetrayPayloadConverter(const Config& config,
                                  std::unique_ptr<const Logger> logger =
                                      getDefaultLogger("DetrayPayloadConverter",
                                                       Logging::INFO));

  std::pair<std::vector<detray::io::grid_payload<
                detray::io::material_slab_payload, detray::io::material_id>>,
            detray::io::material_volume_payload>
  convertMaterial(const TrackingVolume& volume,
                  detray::io::volume_payload& volPayload) const;

 private:
  void handlePortalLink(
      const GeometryContext& gctx, const TrackingVolume& volume,
      detray::io::volume_payload& volPayload,
      const std ::function<std::size_t(const TrackingVolume*)>& volumeLookup,
      const PortalLinkBase& link) const;

  void makeEndOfWorld(const GeometryContext& gctx,
                      detray::io::volume_payload& volPayload,
                      const Surface& surface) const;

  void handlePortal(
      const GeometryContext& gctx, const TrackingVolume& volume,
      detray::io::volume_payload& volPayload,
      const std::function<std::size_t(const TrackingVolume*)>& volumeLookup,
      const Portal& portal) const;

  Config m_cfg;

  const Logger& logger() const { return *m_logger; }
  std::unique_ptr<const Logger> m_logger;
};
}  // namespace Acts
