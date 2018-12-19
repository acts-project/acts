// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// BoundarySurfaceT.h, Acts project
///////////////////////////////////////////////////////////////////

#pragma once
#include <memory>
#include "Acts/Utilities/BinnedArray.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Volumes/BoundarySurfaceFace.hpp"
#include "Acts/Volumes/Volume.hpp"

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
class BoundarySurfaceT
{
  /// declare the TrackingVolume as friend
  friend T;

  using VolumePtr   = std::shared_ptr<const T>;
  using VolumeArray = BinnedArray<VolumePtr>;

public:
  /// Default Constructor
  BoundarySurfaceT()
    : m_surface(nullptr)
    , m_insideVolume(nullptr)
    , m_outsideVolume(nullptr)
    , m_insideVolumeArray(nullptr)
    , m_outsideVolumeArray(nullptr)
  {
  }

  /// Constructor for a Boundary with exact two Volumes attached to it
  /// - usually used in a volume constructor
  ///
  /// @param surface is the unqiue surface the boundary represents
  /// @param inside is the inside volume the bounday surface points to
  /// @param outside is the outside volume the boundary surface points to
  BoundarySurfaceT(std::shared_ptr<const Surface> surface,
                   const T*                       inside,
                   const T*                       outside)
    : m_surface(std::move(surface))
    , m_insideVolume(inside)
    , m_outsideVolume(outside)
    , m_insideVolumeArray(nullptr)
    , m_outsideVolumeArray(nullptr)
  {
  }

  /// Constructor for a Boundary with exact two Volumes attached to it
  /// - usually used in a volume constructor
  ///
  /// @param surface is the unqiue surface the boundary represents
  /// @param inside is the inside volume the bounday surface points to
  /// @param outside is the outside volume the boundary surface points to
  BoundarySurfaceT(std::shared_ptr<const Surface> surface,
                   VolumePtr                      inside,
                   VolumePtr                      outside)
    : m_surface(std::move(surface))
    , m_insideVolume(inside.get())
    , m_outsideVolume(outside.get())
    , m_insideVolumeArray(nullptr)
    , m_outsideVolumeArray(nullptr)
  {
  }

  /// Constructor for a Boundary with exact multiple Volumes attached to it
  /// - usually used in a volume constructor
  ///
  /// @param surface is the unqiue surface the boundary represents
  /// @param insideArray is the inside volume array the bounday surface points
  /// to
  /// @param outsideArray is the outside volume array the boundary surface
  /// points to
  BoundarySurfaceT(std::shared_ptr<const Surface>     surface,
                   std::shared_ptr<const VolumeArray> insideArray,
                   std::shared_ptr<const VolumeArray> outsideArray)
    : m_surface(std::move(surface))
    , m_insideVolume(nullptr)
    , m_outsideVolume(nullptr)
    , m_insideVolumeArray(insideArray)
    , m_outsideVolumeArray(outsideArray)
  {
  }

  /// Get the next Volume depending on GlobalPosition, GlobalMomentum, dir on
  /// the TrackParameters and the requested direction
  ///
  /// @param pos is the global position on surface
  /// @param mom is the direction on the surface
  /// @param dir is an aditional direction corrective
  ///
  /// @return is the attached volume at that position
  virtual const T*
  attachedVolume(const Vector3D&     pos,
                 const Vector3D&     mom,
                 NavigationDirection pdir) const;

  /// templated onBoundary method
  ///
  /// @tparam pars are the parameters to be checked
  template <class P>
  bool
  onBoundary(const P& pars) const
  {
    return surfaceRepresentation().isOnSurface(pars);
  }

  /// The Surface Representation of this
  virtual const Surface&
  surfaceRepresentation() const;

  /// Virtual Destructor
  virtual ~BoundarySurfaceT() = default;

protected:
  /// Helper metho: attach a Volume to this BoundarySurfaceT
  /// this si done during the geometry construction and only called by
  /// the friend templated volume
  ///
  /// @param volume is the volume to be attached
  /// @param inout is the boundary orientation @todo update to along/opposite
  void
  attachVolume(VolumePtr volume, BoundaryOrientation inout);

  /// Helper metho: attach a Volume to this BoundarySurfaceT
  /// this si done during the geometry construction and only called by
  /// the friend templated volume
  ///
  /// @param volumes is the volume array to be attached
  /// @param inout is the boundary orientation @todo update to along/opposite
  void
  attachVolumeArray(std::shared_ptr<const VolumeArray> volumes,
                    BoundaryOrientation                inout);

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
inline const Surface&
BoundarySurfaceT<T>::surfaceRepresentation() const
{
  return (*(m_surface.get()));
}

template <class T>
void
BoundarySurfaceT<T>::attachVolume(VolumePtr volume, BoundaryOrientation inout)
{
  if (inout == insideVolume) {
    m_insideVolume = volume.get();
  } else {
    m_outsideVolume = volume.get();
  }
}

template <class T>
void
BoundarySurfaceT<T>::attachVolumeArray(
    const std::shared_ptr<const VolumeArray> volumes,
    BoundaryOrientation                      inout)
{
  if (inout == insideVolume) {
    m_insideVolumeArray = volumes;
  } else {
    m_outsideVolumeArray = volumes;
  }
}

template <class T>
const T*
BoundarySurfaceT<T>::attachedVolume(const Vector3D&     pos,
                                    const Vector3D&     mom,
                                    NavigationDirection pdir) const
{
  const T* attVolume = nullptr;
  // dot product with normal vector to distinguish inside/outside
  if ((surfaceRepresentation().normal(pos)).dot(pdir * mom) > 0.) {
    attVolume = m_outsideVolumeArray ? m_outsideVolumeArray->object(pos).get()
                                     : m_outsideVolume;
  } else {
    attVolume = m_insideVolumeArray ? m_insideVolumeArray->object(pos).get()
                                    : m_insideVolume;
  }
  return attVolume;
}
}  // namespace Acts
