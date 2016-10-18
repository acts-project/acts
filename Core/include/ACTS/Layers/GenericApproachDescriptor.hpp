// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// GenericApproachDescriptor.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_DETECTOR_GENERICAPPROPACHDESCRIPTOR_H
#define ACTS_DETECTOR_GENERICAPPROPACHDESCRIPTOR_H 1

#include "ACTS/Utilities/ApproachDescriptor.hpp"

namespace Acts {

/// @class GenericApproachDescriptor
///
/// Class to decide and return which approaching surface to be taken,
/// it's a generic descriptor for n surfaces
///
/// It is templated in order to allow for BoundarySurfaces from
/// representing volumes of layers to be re-used

template <class T>
class GenericApproachDescriptor : public ApproachDescriptor
{
public:
  /// A generic approach descriptor for new Acts::Surface objects
  /// passing ownership
  ///
  /// @param aSurfaces are the approach surfaces
  GenericApproachDescriptor(const std::vector<T*>& aSurfaces)
    : ApproachDescriptor(), m_surfaces(), m_surfacesCache(aSurfaces)
  {
    // create the surface container with memory control
    for (auto& sf : (aSurfaces)) m_surfaces.push_back(std::shared_ptr<T>(sf));
  }

  /// A generic approach descriptor with shared surfaces to test
  /// can not be sed with Acts::Surfaces obejcts
  ///
  /// @param aSurfaces are the approach surfaces
  GenericApproachDescriptor(std::vector<std::shared_ptr<T>> aSurfaces)
    : ApproachDescriptor(), m_surfaces(aSurfaces), m_surfacesCache()
  {
    m_surfacesCache.reserve(m_surfaces.size());
    // cache the surfaces
    for (auto& sf : (aSurfaces))
      m_surfacesCache.push_back(&(sf->surfaceRepresentation()));
  }

  /// A generic approach descriptor with n surfaces to test
  ~GenericApproachDescriptor() {}
  /// register the Layer to the surfaces
  ///
  /// @param lay is the layer to be registerd
  void
  registerLayer(const Layer& lay) override;

  /// get the compatible surfaces
  ///
  /// @param gpos is the global posoition to start the approach from
  /// @param dir is the direction in which you approach the layer
  /// @param bchk is the boundary check presrcition
  /// @param ice is a (future) compatibility estimater if needed
  ///
  /// @return : a boolean indicating if an actual intersection had been tried
  const SurfaceIntersection
  approachSurface(const Vector3D&                gpos,
                  const Vector3D&                dir,
                  const BoundaryCheck&           bchk,
                  const ICompatibilityEstimator* ice = nullptr) const override;

  /// return all containes surfaces of this approach descriptor
  const std::vector<const Surface*>&
  containedSurfaces() const override;

private:
  /// approach surfaces with ownership control
  std::vector<std::shared_ptr<T>> m_surfaces;
  /// the surface container cache
  std::vector<const Surface*> m_surfacesCache;
};

template <class T>
void
GenericApproachDescriptor<T>::registerLayer(const Layer& lay)
{
  // go through the surfaces
  for (auto& sf : (m_surfacesCache)) sf->associateLayer(lay);
}

template <class T>
const SurfaceIntersection
GenericApproachDescriptor<T>::approachSurface(
    const Vector3D&      pos,
    const Vector3D&      dir,
    const BoundaryCheck& bchk,
    const ICompatibilityEstimator*) const
{
  // prepare the return surface
  Intersection   sIntersection;
  const Surface* aSurface    = nullptr;
  double         aPathLength = 10e10;
  // get the surfaces
  for (auto& sfIter : m_surfacesCache) {
    // get the intersection with the surface
    sIntersection = sfIter->intersectionEstimate(pos, dir, true, bchk);
    // validatie if it's ok and take the closest
    if (sIntersection.valid && sIntersection.pathLength < aPathLength) {
      aPathLength = sIntersection.pathLength;
      aSurface    = sfIter;
    }
  }
  // return what you have
  return SurfaceIntersection(sIntersection, aSurface, alongMomentum);
}

template <class T>
const std::vector<const Surface*>&
GenericApproachDescriptor<T>::containedSurfaces() const
{
  return m_surfacesCache;
}
}

#endif
