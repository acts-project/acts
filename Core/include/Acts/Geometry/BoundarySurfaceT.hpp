// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// BoundarySurfaceT.h, Acts project
///////////////////////////////////////////////////////////////////

#pragma once
#include <memory>
#include "Acts/Geometry/BoundarySurfaceFace.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/Volume.hpp"
#include "Acts/Utilities/BinnedArray.hpp"
#include "Acts/Utilities/Definitions.hpp"

namespace Acts {

class Surface;

/// @class BoundarySurfaceT
///
/// The boundary surface class combines a Surface with the information of a
/// volume.
/// It's templated in the type of volume in order to allow for a return type tat
/// is usable
/// in the navigation stream.
///
/// @note inside/outside definitions are given by the normal vector of the
/// surface
///
/// @todo change to one schema with BinnedArray0D and along/opposite
/// nominclature

template <class T>
class BoundarySurfaceT {
  /// declare the TrackingVolume as friend
  friend T;

  using VolumePtr = std::shared_ptr<const T>;
  using VolumeArray = BinnedArray<VolumePtr>;

 public:
  /// Default Constructor
  BoundarySurfaceT()
      : m_surface(nullptr),
        m_insideVolume(nullptr),
        m_outsideVolume(nullptr),
        m_insideVolumeArray(nullptr),
        m_outsideVolumeArray(nullptr) {}

  /// Constructor for a Boundary with exact two Volumes attached to it
  /// - usually used in a volume constructor
  ///
  /// @param surface The unqiue surface the boundary represents
  /// @param inside The inside volume the bounday surface points to
  /// @param outside The outside volume the boundary surface points to
  BoundarySurfaceT(std::shared_ptr<const Surface> surface, const T* inside,
                   const T* outside)
      : m_surface(std::move(surface)),
        m_insideVolume(inside),
        m_outsideVolume(outside),
        m_insideVolumeArray(nullptr),
        m_outsideVolumeArray(nullptr) {}

  /// Constructor for a Boundary with exact two Volumes attached to it
  /// - usually used in a volume constructor
  ///
  /// @param surface The unqiue surface the boundary represents
  /// @param inside The inside volume the bounday surface points to
  /// @param outside The outside volume the boundary surface points to
  BoundarySurfaceT(std::shared_ptr<const Surface> surface, VolumePtr inside,
                   VolumePtr outside)
      : m_surface(std::move(surface)),
        m_insideVolume(inside.get()),
        m_outsideVolume(outside.get()),
        m_insideVolumeArray(nullptr),
        m_outsideVolumeArray(nullptr) {}

  /// Constructor for a Boundary with exact multiple Volumes attached to it
  /// - usually used in a volume constructor
  ///
  /// @param surface The unqiue surface the boundary represents
  /// @param insideArray The inside volume array the bounday surface points to
  /// @param outsideArray The outside volume array the boundary surface
  /// points to
  BoundarySurfaceT(std::shared_ptr<const Surface> surface,
                   std::shared_ptr<const VolumeArray> insideArray,
                   std::shared_ptr<const VolumeArray> outsideArray)
      : m_surface(std::move(surface)),
        m_insideVolume(nullptr),
        m_outsideVolume(nullptr),
        m_insideVolumeArray(insideArray),
        m_outsideVolumeArray(outsideArray) {}

  /// Get the next Volume depending on GlobalPosition, GlobalMomentum, dir on
  /// the TrackParameters and the requested direction
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param pos The global position on surface
  /// @param mom The direction on the surface
  /// @param dir is an aditional direction corrective
  ///
  /// @return The attached volume at that position
  virtual const T* attachedVolume(const GeometryContext& gctx,
                                  const Vector3D& pos, const Vector3D& mom,
                                  NavigationDirection pdir) const;

  /// templated onBoundary method
  ///
  /// @tparam pars are the parameters to be checked
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param pars The parameters used for this call
  template <class P>
  bool onBoundary(const GeometryContext& gctx, const P& pars) const {
    return surfaceRepresentation().isOnSurface(gctx, pars);
  }

  /// The Surface Representation of this
  virtual const Surface& surfaceRepresentation() const;

  /// Virtual Destructor
  virtual ~BoundarySurfaceT() = default;

protected:
  /// Helper method: attach a Volume to this BoundarySurfaceT
  /// this is done during the geometry construction and only called by
  /// the friend templated volume
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param volume The volume to be attached
  /// @param inout The boundary orientation @todo update to along/opposite
  void
  attachVolume(const T* volume, BoundaryOrientation inout);

  /// Helper method: attach a Volume to this BoundarySurfaceT
  /// this is done during the geometry construction and only called by
  /// the friend templated volume
  ///
  /// @param volumes The volume array to be attached
  /// @param inout The boundary orientation @todo update to along/opposite
  void attachVolumeArray(std::shared_ptr<const VolumeArray> volumes,
                         BoundaryOrientation inout);

 protected:
  /// the represented surface by this
  std::shared_ptr<const Surface> m_surface;
  /// the inside (w.r.t. normal vector) volume to point to if only one exists
  const T* m_insideVolume;
  /// the outside (w.r.t. normal vector) volume to point to if only one exists
  const T* m_outsideVolume;
  /// the inside (w.r.t. normal vector) volume array to point to
  std::shared_ptr<const VolumeArray> m_insideVolumeArray;
  /// the outside (w.r.t. normal vector) volume array to point to
  std::shared_ptr<const VolumeArray> m_outsideVolumeArray;
};

template <class T>
inline const Surface& BoundarySurfaceT<T>::surfaceRepresentation() const {
  return (*(m_surface.get()));
}

template <class T>
void BoundarySurfaceT<T>::attachVolume(const T* volume,
                                       BoundaryOrientation inout) {
  if (inout == insideVolume) {
    m_insideVolume = volume;
  } else {
    m_outsideVolume = volume;
  }
}

template <class T>
void BoundarySurfaceT<T>::attachVolumeArray(
    const std::shared_ptr<const VolumeArray> volumes,
    BoundaryOrientation inout) {
  if (inout == insideVolume) {
    m_insideVolumeArray = volumes;
  } else {
    m_outsideVolumeArray = volumes;
  }
}

template <class T>
const T* BoundarySurfaceT<T>::attachedVolume(const GeometryContext& gctx,
                                             const Vector3D& pos,
                                             const Vector3D& mom,
                                             NavigationDirection pdir) const {
  const T* attVolume = nullptr;
  // dot product with normal vector to distinguish inside/outside
  if ((surfaceRepresentation().normal(gctx, pos)).dot(pdir * mom) > 0.) {
    attVolume = m_outsideVolumeArray ? m_outsideVolumeArray->object(pos).get()
                                     : m_outsideVolume;
  } else {
    attVolume = m_insideVolumeArray ? m_insideVolumeArray->object(pos).get()
                                    : m_insideVolume;
  }
  return attVolume;
}
}  // namespace Acts
