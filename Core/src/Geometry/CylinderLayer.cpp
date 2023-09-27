// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/CylinderLayer.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/AbstractVolume.hpp"
#include "Acts/Geometry/BoundarySurfaceFace.hpp"
#include "Acts/Geometry/BoundarySurfaceT.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/GenericApproachDescriptor.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Helpers.hpp"

#include <vector>

using Acts::VectorHelpers::phi;

Acts::CylinderLayer::CylinderLayer(
    const Transform3& transform,
    const std::shared_ptr<const CylinderBounds>& cBounds,
    std::unique_ptr<SurfaceArray> surfaceArray, double thickness,
    std::unique_ptr<ApproachDescriptor> ades, LayerType laytyp)
    : CylinderSurface(transform, cBounds),
      Layer(std::move(surfaceArray), thickness, std::move(ades), laytyp) {
  // create the representing volume
  auto cVolumeBounds = std::make_shared<const CylinderVolumeBounds>(
      *CylinderSurface::m_bounds, thickness);
  // @todo rotate around x for the avePhi if you have a sector
  m_representingVolume =
      std::make_unique<AbstractVolume>(m_transform, cVolumeBounds);

  // associate the layer to the surface
  CylinderSurface::associateLayer(*this);
  // an approach descriptor is automatically created if there's a surface array
  if (!m_approachDescriptor && m_surfaceArray) {
    buildApproachDescriptor();
  }
  // register the layer to the approach descriptor surfaces
  if (m_approachDescriptor) {
    approachDescriptor()->registerLayer(*this);
  }
}

const Acts::CylinderSurface& Acts::CylinderLayer::surfaceRepresentation()
    const {
  return (*this);
}

Acts::CylinderSurface& Acts::CylinderLayer::surfaceRepresentation() {
  return (*this);
}

void Acts::CylinderLayer::buildApproachDescriptor() {
  // delete and reset as you build a new one
  m_approachDescriptor.reset(nullptr);

  // take the boundary surfaces of the representving volume if they exist
  if (m_representingVolume != nullptr) {
    // get the boundary surfaces
    const std::vector<std::shared_ptr<const BoundarySurfaceT<AbstractVolume>>>&
        bSurfaces = m_representingVolume->boundarySurfaces();

    // fill in the surfaces into the vector
    std::vector<std::shared_ptr<const Surface>> aSurfaces;
    if (bSurfaces.size() > size_t(tubeInnerCover)) {
      aSurfaces.push_back(
          bSurfaces.at(tubeInnerCover)->surfaceRepresentation().getSharedPtr());
    }
    aSurfaces.push_back(
        bSurfaces.at(tubeOuterCover)->surfaceRepresentation().getSharedPtr());
    // create an ApproachDescriptor with Boundary surfaces
    m_approachDescriptor =
        std::make_unique<const GenericApproachDescriptor>(std::move(aSurfaces));
  }

  for (auto& sfPtr : (m_approachDescriptor->containedSurfaces())) {
    if (sfPtr != nullptr) {
      auto& mutableSf = *(const_cast<Surface*>(sfPtr));
      mutableSf.associateLayer(*this);
    }
  }
}
