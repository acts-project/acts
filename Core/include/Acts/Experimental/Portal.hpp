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
#include "Acts/Experimental/NavigationState.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Surfaces/BoundaryCheck.hpp"
#include "Acts/Surfaces/Surface.hpp"

#include <memory>
#include <optional>
#include <unordered_set>

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

  /// Access to the surface representation
  const Surface& surface() const;

  /// Update switching to next volume on this portal
  ///
  /// @param gctx is the current geometry context
  /// @param position is the position at the query
  /// @param direction is the direction at the query
  ///
  const DetectorVolume* nextVolume(const GeometryContext& gctx,
                                   const Vector3& position,
                                   const Vector3& direction) const;

  /// Assign the surface material description
  ///
  /// The material is usually derived in a complicated way and loaded from
  /// a framework given source. As various surfaces may share the same source
  /// this is provided by a shared pointer
  ///
  /// @param material Material description associated to this surface
  void assignSurfaceMaterial(std::shared_ptr<const ISurfaceMaterial> material);

  /// Set the geometry identifier (to the underlying surface)
  ///
  /// @param geometryId the geometry identifier to be assigned
  void assignGeometryId(const GeometryIdentifier& geometryId);

  /// Update the volume link
  ///
  /// @param nDir the navigation direction for the link
  /// @param dVolumeLink is the link delegate
  /// @param dVolumeLinkStore is the link implementation store
  /// @param crationVolume is a boolean to flag if this was the
  /// creation volume
  ///
  /// @note this overwrites the existing link
  void updateVolumeLink(NavigationDirection nDir,
                        const DetectorVolumeLink& dVolumeLink,
                        DetectorVolumeLinkStore dVolumeLinkStore = nullptr,
                        bool creationVolume = false);

  /// Assign an outside volume link, i.e. the none-creation volume
  ///
  /// @param dVolumeLink is the link delegate
  /// @param dVolumeLinkStore is the link implementation store
  /// @param overwrite steers whether any existing link should
  ///        be overwritten
  ///
  /// @note this relies on the m_creationVolumedir optional to be set
  void updateOutsideVolumeLink(
      const DetectorVolumeLink& dVolumeLink,
      DetectorVolumeLinkStore dVolumeLinkStore = nullptr,
      bool overwrite = false);

 private:
  /// The surface representation of this portal
  std::shared_ptr<Surface> m_surface;

  /// Navigation direction of the creation volume (optionally set)
  /// this can  help with the glueing process
  std::optional<NavigationDirection> m_creationVolumeDir = std::nullopt;

  /// The environment updators - along normal direction
  DetectorVolumeLink m_forwardLink;

  /// The environment updators - opposite normal direction
  DetectorVolumeLink m_backwardLink;

  /// The store for potential lifetime control of the links
  std::unordered_set<DetectorVolumeLinkStore> m_linkStore;
};

inline const Surface& Portal::surface() const {
  return *(m_surface.get());
}

}  // namespace Experimental
}  // namespace Acts
