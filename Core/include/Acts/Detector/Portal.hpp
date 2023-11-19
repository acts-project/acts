// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Direction.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Navigation/NavigationDelegates.hpp"
#include "Acts/Navigation/NavigationState.hpp"
#include "Acts/Surfaces/BoundaryCheck.hpp"
#include "Acts/Surfaces/RegularSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"

#include <array>
#include <map>
#include <memory>
#include <optional>
#include <vector>

namespace Acts {

class ISurfaceMaterial;
class Surface;

namespace Experimental {
class DetectorVolume;
struct NavigationState;

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
  Portal(std::shared_ptr<RegularSurface> surface);

 public:
  /// The volume links forward/backward with respect to the surface normal
  using DetectorVolumeUpdators = std::array<DetectorVolumeUpdator, 2u>;

  /// The vector of attached volumes forward/backward, this is useful in the
  /// geometry building
  using AttachedDetectorVolumes =
      std::array<std::vector<std::shared_ptr<DetectorVolume>>, 2u>;

  /// Declare the DetectorVolume friend for portal setting
  friend class DetectorVolume;

  /// Factory for producing memory managed instances of Portal.
  static std::shared_ptr<Portal> makeShared(
      std::shared_ptr<RegularSurface> surface);

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
  const RegularSurface& surface() const;

  /// Non-const access to the surface reference
  RegularSurface& surface();

  /// Update the current volume
  ///
  /// @param gctx is the Geometry context of this call
  /// @param nState [in,out] the navigation state for the volume to be updated
  ///
  void updateDetectorVolume(const GeometryContext& gctx,
                            NavigationState& nState) const noexcept(false);

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
  ///
  /// @note that other will be overwritten to point to this
  void fuse(std::shared_ptr<Portal>& other) noexcept(false);

  /// Update the volume link
  ///
  /// @param dir the direction of the link
  /// @param dVolumeUpdator is the mangaged volume updator delegate
  /// @param attachedVolumes is the list of attached volumes for book keeping
  ///
  /// @note this overwrites the existing link
  void assignDetectorVolumeUpdator(
      Direction dir, DetectorVolumeUpdator dVolumeUpdator,
      std::vector<std::shared_ptr<DetectorVolume>> attachedVolumes);

  /// Update the volume link, w/o directive, i.e. it relies that there's only
  /// one remaining link to be set, throws an exception if that's not the case
  ///
  /// @param dVolumeUpdator is the mangaged volume updator delegate
  /// @param attachedVolumes is the list of attached volumes for book keeping
  ///
  /// @note this overwrites the existing link
  void assignDetectorVolumeUpdator(DetectorVolumeUpdator dVolumeUpdator,
                                   std::vector<std::shared_ptr<DetectorVolume>>
                                       attachedVolumes) noexcept(false);

  // Access to the portal targets: opposite/along normal vector
  const DetectorVolumeUpdators& detectorVolumeUpdators() const;

  // Access to the attached volumes - non-const access
  AttachedDetectorVolumes& attachedDetectorVolumes();

 private:
  /// The surface representation of this portal
  std::shared_ptr<RegularSurface> m_surface;

  /// The portal targets along/opposite the normal vector
  DetectorVolumeUpdators m_volumeUpdators = {unconnectedUpdator(),
                                             unconnectedUpdator()};

  /// The portal attaches to the following volumes
  AttachedDetectorVolumes m_attachedVolumes;
};

}  // namespace Experimental
}  // namespace Acts
