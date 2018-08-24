// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// GenericApproachDescriptor.h, Acts project
///////////////////////////////////////////////////////////////////

#pragma once

#include <algorithm>
#include "Acts/Utilities/ApproachDescriptor.hpp"
#include "Acts/Utilities/Intersection.hpp"

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
  GenericApproachDescriptor(const std::vector<const T*>& aSurfaces)
    : ApproachDescriptor(), m_surfaces(), m_surfaceCache(aSurfaces)
  {
    // create the surface container with memory control
    for (const auto& sf : aSurfaces) {
      m_surfaces.push_back(std::shared_ptr<const T>(sf));
    }
  }

  /// A generic approach descriptor with shared surfaces to test
  /// can not be used with Acts::Surfaces objects
  ///
  /// @param aSurfaces are the approach surfaces
  GenericApproachDescriptor(std::vector<std::shared_ptr<const T>> aSurfaces)
    : ApproachDescriptor(), m_surfaces(aSurfaces), m_surfaceCache()
  {
    m_surfaceCache.reserve(m_surfaces.size());
    // cache the surfaces
    for (const auto& sf : aSurfaces) {
      m_surfaceCache.push_back(&(sf->surfaceRepresentation()));
    }
  }

  /// A generic approach descriptor with n surfaces to test
  ~GenericApproachDescriptor() override {}
  /// register the Layer to the surfaces
  ///
  /// @param lay is the layer to be registerd
  void
  registerLayer(const Layer& lay) override;

  /// get the compatible surfaces
  ///
  /// @param gpos is the global position to start the approach from
  /// @param mom is the momentum vector
  /// @param bcheck is the boundary check prescription
  /// @param corrfnc is an noption correction function
  ///
  /// @return : a boolean indicating if an actual intersection had been tried
  ObjectIntersection<Surface>
  approachSurface(const Vector3D&      gpos,
                  const Vector3D&      mom,
                  NavigationDirection  navDir,
                  const BoundaryCheck& bcheck,
                  CorrFnc              corrfnc = nullptr) const override;

  /// return all contained surfaces of this approach descriptor
  const std::vector<const Surface*>&
  containedSurfaces() const override;

  /// Non-const version
  std::vector<const Surface*>&
  containedSurfaces() override;

private:
  /// approach surfaces with ownership control
  std::vector<std::shared_ptr<const T>> m_surfaces;

  /// the surface container cache
  ///
  /// We will need to mutate those surfaces in registerLayer, but the C++ type
  /// system has no const-correct way of expressing this constraint.
  ///
  std::vector<const Surface*> m_surfaceCache;
};

template <class T>
void
GenericApproachDescriptor<T>::registerLayer(const Layer& lay)
{
  // go through the surfaces
  for (auto& sf : m_surfaceCache) {
    auto mutableSf = const_cast<Surface*>(sf);
    mutableSf->associateLayer(lay);
  }
}

template <class T>
ObjectIntersection<Surface>
GenericApproachDescriptor<T>::approachSurface(const Vector3D&      gpos,
                                              const Vector3D&      gdir,
                                              NavigationDirection  navDir,
                                              const BoundaryCheck& bcheck,
                                              CorrFnc corrfnc) const
{
  // the intersection estimates
  std::vector<ObjectIntersection<Surface>> sIntersections;
  sIntersections.reserve(m_surfaceCache.size());
  for (auto& sf : m_surfaceCache) {
    // intersect
    auto intersection
        = sf->intersectionEstimate(gpos, gdir, navDir, bcheck, corrfnc);
    sIntersections.push_back(
        ObjectIntersection<Surface>(intersection, sf, navDir));
  }
  // sort depedning on the navigation direction
  if (navDir == forward) {
    std::sort(sIntersections.begin(), sIntersections.end());
  } else {
    std::sort(sIntersections.begin(), sIntersections.end(), std::greater<>());
  }
  // return what you have
  return (*sIntersections.begin());
}

template <class T>
const std::vector<const Surface*>&
GenericApproachDescriptor<T>::containedSurfaces() const
{
  return m_surfaceCache;
}

template <class T>
std::vector<const Surface*>&
GenericApproachDescriptor<T>::containedSurfaces()
{
  return m_surfaceCache;
}

}  // namespace Acts