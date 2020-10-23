// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/GenericApproachDescriptor.hpp"

#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Intersection.hpp"

#include <algorithm>

void Acts::GenericApproachDescriptor::registerLayer(const Layer& lay) {
  // go through the surfaces
  for (auto& sf : m_surfaceCache) {
    auto mutableSf = const_cast<Surface*>(sf);
    mutableSf->associateLayer(lay);
  }
}

Acts::ObjectIntersection<Acts::Surface>
Acts::GenericApproachDescriptor::approachSurface(
    const GeometryContext& gctx, const Vector3D& position,
    const Vector3D& direction, const BoundaryCheck& bcheck) const {
  std::vector<ObjectIntersection<Surface>> sIntersections;
  sIntersections.reserve(m_surfaceCache.size());
  for (auto& sf : m_surfaceCache) {
    auto sfIntersection = sf->intersect(gctx, position, direction, bcheck);
    // Overstepping is not allowed for approach surfaces
    if (sfIntersection.intersection.pathLength < 0. and
        sfIntersection.alternative.pathLength > 0.) {
      std::swap(sfIntersection.intersection, sfIntersection.alternative);
    }
    sIntersections.push_back(sfIntersection);
  }
  // Sort them & return the closest
  std::sort(sIntersections.begin(), sIntersections.end());
  return (*sIntersections.begin());
}

const std::vector<const Acts::Surface*>&
Acts::GenericApproachDescriptor::containedSurfaces() const {
  return m_surfaceCache;
}

std::vector<const Acts::Surface*>&
Acts::GenericApproachDescriptor::containedSurfaces() {
  return m_surfaceCache;
}
