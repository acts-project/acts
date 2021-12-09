// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Common.hpp"
#include "Acts/Experimental/DetectorEnvironment.hpp"
#include "Acts/Experimental/Portal.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/VolumeBounds.hpp"
#include "Acts/Surfaces/BoundaryCheck.hpp"

#include <memory>
#include <stdexcept>
#include <string>
#include <vector>
#include <array>

/// @note This file is foreseen for the `Geometry` module

namespace Acts {

class Surface;
class DetectorVolume;

using SurfaceLinks = std::function<std::vector<SurfaceIntersection>(
    const GeometryContext&, const DetectorVolume&, const Vector3&,
    const Vector3&, const BoundaryCheck&, const std::array<ActsScalar,2>&, bool)>;

/// A detector volume description which can be:
///
/// - a) a layer volume
/// - b) a container volume
/// - c) a gap voume
///
/// @note A detector volume holds non-const objects internally
/// that are allowed to be modified as long as the geometry
/// is not yet closed. Like this, material can be loaded,
/// and GeometryId can be set.
/// The access the objects is given only as const access.
///
/// @note Navigation is always done by plain pointers, while
/// object ownership is done by shared/unique pointers.
class DetectorVolume {
 public:
  /// Nested object store that holds the internal (non-const),
  /// reference counted objects and provides an external
  /// (const raw pointer) access
  ///
  /// @tparam Internal is the internal storage representation,
  /// has to comply with std::shared_ptr or std::unique_ptr
  /// semantics
  template <typename Internal>
  struct ObjectStore {
    /// The internal storage vector
    std::vector<Internal> internal = {};

    /// The external storage vector, const raw pointer
    std::vector<const typename Internal::element_type*> external = {};

    /// Store constructor
    ///
    /// @param objects are the ones copied into the internal store
    ObjectStore(const std::vector<Internal>& objects)
        : internal(objects) {
      external = unpack_ref_to_const_vector(internal);
    }

    ObjectStore() = default;
  };

 private:
  /// Create a detector volume - layer volume constructor
  ///
  /// @param transform the transform defining the volume position
  /// @param bounds the volume bounds
  /// @param surfaces the contained surfaces 
  /// @param volumeSurfaceLinks the atacched links to surface (from volume)
  /// @param portalSurfaceLinks the attached links to surfaces (from portal)
  /// @param name the volume name
  ///
  DetectorVolume(const Transform3& transform,
                 std::unique_ptr<VolumeBounds> bounds,
                 const std::vector<std::shared_ptr<Surface>>& surfaces,
                 SurfaceLinks&& volumeSurfaceLinks = VoidSurfaceLink{},
                 std::vector<SurfaceLinks>&& portalSurfaceLinks = {},
                 const std::string& name = "Unnamed");

  /// Create a detector volume - container volume constructor
  ///
  /// @param transform the transform defining the volume position
  /// @param bounds the volume bounds
  /// @param volumes the contained volumes
  /// @param volumeLink the links for finding the volumes
  /// @param recastPortals boolean flag that indicates to recast portals
  /// @param recastValue the recast direction for portals
  /// @param name the volume name
  ///
  DetectorVolume(const Transform3& transform,
                 std::unique_ptr<VolumeBounds> bounds,
                 const std::vector<std::shared_ptr<DetectorVolume>>& volumes,
                 VolumeLink&& volumeLink, bool recastPortals,
                 BinningValue recastValue, const std::string& name = "Unnamed");

  /// Create a detector volume - gap volume constructor
  ///
  /// @param transform the transform defining the volume position
  /// @param bounds the volume bounds
  /// @param name the volume name
  ///
  DetectorVolume(const Transform3& transform,
                 std::unique_ptr<VolumeBounds> bounds,
                 const std::string& name = "Unnamed");

 public:
  /// Factory for producing memory managed instances of Surface.
  /// Will forward all parameters and will attempt to find a suitable
  /// constructor.
  ///
  /// @tparam Args the arguments that will be forwarded
  template <typename... Args>
  static std::shared_ptr<DetectorVolume> makeShared(Args&&... args) {
    return std::shared_ptr<DetectorVolume>(
        new DetectorVolume(std::forward<Args>(args)...));
  }

  /// Const access to the transform
  ///
  /// @param gctx the geometry contect (@note currently ignored)
  ///
  /// @return const reference to the contextual transform
  const Transform3& transform(
      const GeometryContext& gctx = GeometryContext()) const;

  /// Const access to the volume bounds
  ///
  /// @return const reference to the volume bounds object
  const VolumeBounds& volumeBounds() const;

  /// Inside/outside method
  ///
  /// @param position the position for the inside check
  /// @param tolerance is the tolerance parameter
  ///
  /// @return a bool to indicate inside/outside
  bool inside(const Vector3& position, ActsScalar tolerance = 0.) const;

  /// Non-const access to the portals
  ///
  /// @return the portal shared pointer store
  std::vector<std::shared_ptr<Portal>>& portalPtrs();

  /// Update of a portal
  ///
  /// @note this will re-initialize the portal store
  /// @note it will throw an exception if something is wrong
  ///
  /// @param updatedPortal the new portal
  /// @param portal the portal position
  /// @param keepPortalLink is a directive to either disregard/keep a volume link
  /// @param nDir the navigation direction of the portal link
  ///             to be kept
  void updatePortalPtr(std::shared_ptr<Portal> updatedPortal, size_t portal,
                       bool keepPortalLink = false,
                       NavigationDirection nDir = forward) noexcept(false);

  /// Get a new detector environment
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param position The global position on surface
  /// @param direction The direction on the surface
  /// @param pathRange The possible path range
  /// @param bCheck is the boundary check for the surface search
  /// @param provideAll is the boolean switch for trial&error navigation
  ///
  /// @return a new detector environment with portals and surfaces
  DetectorEnvironment environment(const GeometryContext& gctx,
                                  const Vector3& position,
                                  const Vector3& direction,
                                  const std::array<ActsScalar,2>& pathRange,
                                  const BoundaryCheck& bCheck,
                                  bool provideAll = false) const;

  /// Lowest volume in at navigation level
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  ///
  /// @return the loweset volume
  const DetectorVolume* lowest(const GeometryContext& gctx,
                               const Vector3& position) const;

  /// Const access to the detector portals
  ///
  /// @note an empty vector indicates a container volume
  /// that has not been properly connected
  ///
  /// @return a vector to const Portal raw pointers
  const std::vector<const Portal*>& portals() const;

  /// Const acess to the surfaces
  ///
  /// @note an empty vector indicates either gap volume
  /// or container volume, a non-empty vector indicates
  /// a layer volume.
  ///
  /// @return a vector to const Surface raw pointers
  const std::vector<const Surface*>& surfaces() const;

  /// Const access to sub volumes
  ///
  /// @note and empty vector indicates this is either a
  /// gap volume or a layer volume, in any case it means
  /// the volume is on navigation level and the portals
  /// need to be connected
  ///
  /// @return a vector to const DetectorVolume raw pointers
  const std::vector<const DetectorVolume*>& volumes() const;

  /// Set the name ov the volume
  void setName(const std::string& name);

  /// @return the name of the volume
  const std::string& name() const;

 private:
  /// Private method to build portal surfaces
  ///
  /// @param surfaceLinks is the list of surface links to be set
  /// @param recastPortals is a boolean flag to recast
  /// @param recastValue is the BinningValue for the portal recasting
  ///
  /// @note if the vector is not empty but the length does not correspond to
  /// the number of portals an exception will be thrown
  void createPortals(std::vector<SurfaceLinks>&& surfaceLinks = {},
                     bool recastPortals = false,
                     BinningValue recastValue = binX) noexcept(false);

  /// Transform to place the bolume
  Transform3 m_transform = Transform3::Identity();

  /// Volume boundaries
  std::unique_ptr<VolumeBounds> m_bounds = nullptr;

  /// Portal store (internal/external)
  ObjectStore<std::shared_ptr<Portal>> m_portals;

  /// Surface store (internal/external)
  ObjectStore<std::shared_ptr<Surface>> m_surfaces;
  /// Volume internal access to the surfaces
  SurfaceLinks m_surfaceLinks = VoidSurfaceLink{};

  /// Volume store (internal/external)
  ObjectStore<std::shared_ptr<DetectorVolume>> m_volumes;
  /// Volume internal access to the volumes
  VolumeLink m_volumeLink = VoidVolumeLink{};

  /// Name of the volume
  std::string m_name = "Unnamed";
};

inline const Transform3& DetectorVolume::transform(
    const GeometryContext& /*gctx*/) const {
  return m_transform;
}

inline const VolumeBounds& DetectorVolume::volumeBounds() const {
  return (*m_bounds.get());
}

inline std::vector<std::shared_ptr<Portal>>& DetectorVolume::portalPtrs() {
  return m_portals.internal;
}

inline const std::vector<const Portal*>& DetectorVolume::portals() const {
  return m_portals.external;
}

inline const std::vector<const Surface*>& DetectorVolume::surfaces() const {
  return m_surfaces.external;
}

inline const std::vector<const DetectorVolume*>& DetectorVolume::volumes()
    const {
  return m_volumes.external;
}

inline void DetectorVolume::setName(const std::string& name) {
  m_name = name;
}


inline const std::string& DetectorVolume::name() const {
  return m_name;
}

}  // namespace Acts
