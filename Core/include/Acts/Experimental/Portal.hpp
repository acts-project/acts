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
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Surfaces/BoundaryCheck.hpp"
#include "Acts/Utilities/Intersection.hpp"

#include <array>
#include <functional>
#include <memory>
#include <vector>

namespace Acts {

class ISurfaceMaterial;
class DetectorVolume;
class Portal;
class Surface;

/// Definition of a SurfaceLinks function
///
using SurfaceLinks = std::function<std::vector<SurfaceIntersection>(
    const GeometryContext&, const DetectorVolume&, const Vector3&,
    const Vector3&, const BoundaryCheck&, const std::array<ActsScalar, 2>&,
    bool)>;

/// Decleare a void surface link that
struct VoidSurfaceLink {
  /// This returns an emtpy vector
  ///
  /// @note the input parameters are ignored
  std::vector<SurfaceIntersection> operator()(
      const GeometryContext& /*gctx*/, const DetectorVolume& /*volume*/,
      const Vector3& /*position*/, const Vector3& /*direction*/,
      const BoundaryCheck& /*bCheck*/,
      const std::array<ActsScalar, 2>& /*pathRange*/, bool) const {
    return {};
  }
};

/// Definition of a VolumeLink function
///
using VolumeLink =
    std::function<unsigned int(const Vector3&)>;

/// Declare a void volume index
struct VoidVolumeLink {
  /// @note the parameters are ignored in this context
  unsigned int operator()(const Vector3&) const {
    return 0u;
  }
};

/// Definition of an entry link
///
/// @param gctx is the current geometry context
/// @param portal the portal at this request
/// @param position is the current position on a portal
/// @param direction is the current direction at that portal
/// @param provideAll is a flag to switch on test-all navigation
///
using PortalLink = std::function<DetectorEnvironment(
    const GeometryContext& gctx, const Portal& portal, const Vector3& position,
    const Vector3& direction, const BoundaryCheck& bCheck, bool provideAll)>;

/// A void entry link
struct VoidPortalLink {
  /// Fullfills the call std::function call structure
  ///
  /// @note parameters are ignored
  DetectorEnvironment operator()(const GeometryContext& /*gctx*/,
                                 const Portal& /*portal*/,
                                 const Vector3& /*position*/,
                                 const Vector3& /*direction*/,
                                 const BoundaryCheck& /*bCheck*/,
                                 bool /*provideAll = false*/) const { 
    return DetectorEnvironment{};
  }
};

/// A portal between the detector volumes
///
/// It has a Surface representation for navigation and propagation
/// and guides into the next volumes.
///
/// The surface can also carry material to allow mapping onto
/// portal positions.
///
class Portal {
 public:
  /// Declare the DetectorVolume friend for portal setting
  friend class DetectorVolume;

  /// Constructor with argument
  ///
  /// @param surface is the representing surface
  Portal(std::shared_ptr<Surface> surface);

  Portal() = delete;
  virtual ~Portal() = default;

  /// Access to the surface representation
  const Surface& surfaceRepresentation() const;

  /// Portal intersection
  ///
  /// @param gctx is the current geometry conbtext
  /// @param position is the position at the query
  /// @param direction is the direction at the query
  ///
  PortalIntersection intersect(const GeometryContext& gctx,
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

  /// Update the portal link - move semantics, this is
  /// with respect to the normal vector of the surface
  ///
  /// @param portalLink the volume link to be updated
  /// @param nDir the navigation direction
  void updatePortalLink(PortalLink&& portalLink, NavigationDirection nDir);

  /// Retrieve the portalLink given the navigation direction
  ///
  /// @param nDir the navigation direction
  ///
  /// @return the portal link as a const object
  const PortalLink& portalLink(NavigationDirection nDir) const;

  /// Get the next detector environment once you have reached a portal
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param position The global position on surface
  /// @param direction The direction on the surface
  /// @param bCheck is the boundary check for the surface search
  /// @param provideAll is a flag to switch on trial&error navigation
  ///
  /// @return The attached detector volume at that position
  DetectorEnvironment next(const GeometryContext& gctx, const Vector3& position,
                           const Vector3& direction,
                           const BoundaryCheck& bCheck,
                           bool provideAll = false) const;

  /// Set the geometry identifier (to the underlying surface)
  ///
  /// @param geometryId the geometry identifier to be assigned
  void assignGeometryId(const GeometryIdentifier& geometryId);

 private:
  /// The surface representation of this portal
  std::shared_ptr<Surface> m_surface;
  /// The entry link along the surface normal direction
  PortalLink m_alongNormal = VoidPortalLink{};
  /// The entry link opposite the surfacea normal direction
  PortalLink m_oppositeNormal = VoidPortalLink{};
};

inline const Surface& Portal::surfaceRepresentation() const {
  return *(m_surface.get());
}

inline void Portal::updatePortalLink(PortalLink&& portalLink,
                                     NavigationDirection nDir) {
  if (nDir == forward) {
    m_alongNormal = std::move(portalLink);
  } else {
    m_oppositeNormal = std::move(portalLink);
  }
}

inline const PortalLink& Portal::portalLink(NavigationDirection nDir) const {
  if (nDir == forward) {
    return m_alongNormal;
  }
  return m_oppositeNormal;
}

}  // namespace Acts
