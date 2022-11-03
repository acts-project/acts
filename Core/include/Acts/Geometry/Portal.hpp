// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Common.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Geometry/NavigationDelegates.hpp"
#include "Acts/Geometry/NavigationState.hpp"
#include "Acts/Surfaces/BoundaryCheck.hpp"
#include "Acts/Surfaces/Surface.hpp"

#include <map>
#include <memory>
#include <optional>

namespace Acts {

class ISurfaceMaterial;

namespace Experimental {

/// A portal description between the detector volumes
///
/// It has a Surface representation for navigation and propagation
/// and guides from one volume to the next.
///
/// The surface can carry material to allow mapping onto
/// portal positions if required.
///
class Portal : public std::enable_shared_from_this<Portal> {
 protected:
  /// Constructor from surface w/o portal links
  ///
  /// @param surface is the representing surface
  Portal(std::shared_ptr<Surface> surface);

 public:
  /// The volume links forward/backward with respect to the surface normal
  using VolumeLinks = std::array<ManagedDetectorVolumeLink, 2u>;

  /// The vector of attached volumes forward/backward, this is useful in the
  /// geometry building
  using AttachedVolumes =
      std::array<std::vector<std::shared_ptr<DetectorVolume>>, 2u>;

  /// Declare the DetectorVolume friend for portal setting
  friend class DetectorVolume;

  /// Factory for producing memory managed instances of Portal.
  /// Will forward all parameters and will attempt to find a suitable
  /// constructor.
  template <typename... Args>
  static std::shared_ptr<Portal> makeShared(Args&&... args) {
    return std::shared_ptr<Portal>(new Portal(std::forward<Args>(args)...));
  }

  /// Retrieve a @c std::shared_ptr for this surface (non-const version)
  ///
  /// @note Will error if this was not created through the @c makeShared factory
  ///       since it needs access to the original reference. In C++14 this is
  ///       undefined behavior (but most likely implemented as a @c bad_weak_ptr
  ///       exception), in C++17 it is defined as that exception.
  /// @note Only call this if you need shared ownership of this object.
  ///
  /// @return The shared pointer
  std::shared_ptr<Portal> getSharedPtr();

  /// Retrieve a @c std::shared_ptr for this surface (const version)
  ///
  /// @note Will error if this was not created through the @c makeShared factory
  ///       since it needs access to the original reference. In C++14 this is
  ///       undefined behavior, but most likely implemented as a @c bad_weak_ptr
  ///       exception, in C++17 it is defined as that exception.
  /// @note Only call this if you need shared ownership of this object.
  ///
  /// @return The shared pointer
  std::shared_ptr<const Portal> getSharedPtr() const;

  Portal() = delete;
  virtual ~Portal() = default;

  /// Const access to the surface representation
  const Surface& surface() const;

  /// Non-const access to the surface
  std::shared_ptr<Surface> surfacePtr();

  /// Update switching to next volume on this portal
  ///
  /// @param gctx is the current geometry context
  /// @param position is the position at the query
  /// @param direction is the direction at the query
  ///
  /// @return a raw pointer to the next volume, a nullptr
  ///         if no next volume is found, e.g. end of world
  const DetectorVolume* nextVolume(const GeometryContext& gctx,
                                   const Vector3& position,
                                   const Vector3& direction) const;

  /// Set the geometry identifier (to the underlying surface)
  ///
  /// @param geometryId the geometry identifier to be assigned
  void assignGeometryId(const GeometryIdentifier& geometryId);

  /// Fuse with another portal, this one is kept
  ///
  /// @param other is the portal that will be fused
  ///
  /// @note this will move the portal links from the other
  /// into this volume, it will throw an exception if the
  /// portals are not fusable
  void fuse(Portal& other) noexcept(false);

  /// Update the volume link
  ///
  /// @param nDir the navigation direction for the link
  /// @param dVolumeLink is the mangaged link delegate
  /// @param attachedVolumes is the list of attached volumes for book keeping
  ///
  /// @note this overwrites the existing link
  void updateVolumeLink(
      NavigationDirection nDir, ManagedDetectorVolumeLink&& dVolumeLink,
      const std::vector<std::shared_ptr<DetectorVolume>>& attachedVolumes);

  // Access to the portal targets
  const VolumeLinks& volumeLinks() const;

  // Access to the attached volumes - non-const access
  AttachedVolumes& attachedVolumes();

 private:
  /// The surface representation of this portal
  std::shared_ptr<Surface> m_surface;

  /// The portal targets along/opposite the normal vector
  VolumeLinks m_volumeLinks;

  /// The portal attaches to the following volumes
  AttachedVolumes m_attachedVolumes;
};

inline const Surface& Portal::surface() const {
  return *(m_surface.get());
}

inline std::shared_ptr<Surface> Portal::surfacePtr() {
  return m_surface;
}

inline const Portal::VolumeLinks& Portal::volumeLinks() const {
  return m_volumeLinks;
}

inline Portal::AttachedVolumes& Portal::attachedVolumes() {
  return m_attachedVolumes;
}

}  // namespace Experimental
}  // namespace Acts
