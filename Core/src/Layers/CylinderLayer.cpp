// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// CylinderLayer.cpp, Acts project
///////////////////////////////////////////////////////////////////

#include "Acts/Layers/CylinderLayer.hpp"
#include "Acts/Layers/GenericApproachDescriptor.hpp"
#include "Acts/Material/SurfaceMaterial.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Utilities/BinUtility.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Volumes/AbstractVolume.hpp"
#include "Acts/Volumes/BoundarySurfaceFace.hpp"
#include "Acts/Volumes/CylinderVolumeBounds.hpp"

using Acts::VectorHelpers::phi;

Acts::CylinderLayer::CylinderLayer(
    const std::shared_ptr<const Transform3D>&    transform,
    const std::shared_ptr<const CylinderBounds>& cBounds,
    std::unique_ptr<SurfaceArray>                surfaceArray,
    double                                       thickness,
    std::unique_ptr<ApproachDescriptor>          ades,
    LayerType                                    laytyp)
  : CylinderSurface(transform, cBounds)
  , Layer(std::move(surfaceArray), thickness, std::move(ades), laytyp)
{
  // create the representing volume
  CylinderVolumeBounds* cvBounds
      = new CylinderVolumeBounds(cBounds->r() - 0.5 * thickness,
                                 cBounds->r() + 0.5 * thickness,
                                 cBounds->halflengthZ());
  m_representingVolume
      = new AbstractVolume(transform, VolumeBoundsPtr(cvBounds));
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

const Acts::CylinderSurface&
Acts::CylinderLayer::surfaceRepresentation() const
{
  return (*this);
}

Acts::CylinderSurface&
Acts::CylinderLayer::surfaceRepresentation()
{
  return (*this);
}

void
Acts::CylinderLayer::buildApproachDescriptor()
{
  // delete it
  m_approachDescriptor.reset(nullptr);
  // delete the surfaces
  // take the boundary surfaces of the representving volume if they exist
  if (m_representingVolume != nullptr) {
    // get the boundary surfaces
    const std::vector<std::shared_ptr<const BoundarySurfaceT<AbstractVolume>>>&
        bSurfaces
        = m_representingVolume->boundarySurfaces();

    // fill in the surfaces into the vector
    std::vector<std::shared_ptr<const Surface>> aSurfaces;
    if (bSurfaces.size() > size_t(tubeOuterCover)) {
      aSurfaces.push_back(
          bSurfaces.at(tubeInnerCover)->surfaceRepresentation().getSharedPtr());
    }
    aSurfaces.push_back(
        bSurfaces.at(tubeOuterCover)->surfaceRepresentation().getSharedPtr());
    // create an ApproachDescriptor with Boundary surfaces
    m_approachDescriptor = std::make_unique<const GenericApproachDescriptor>(
        std::move(aSurfaces));
  } else {
    // create the new surfaces
    std::vector<std::shared_ptr<const Acts::Surface>> aSurfaces;
    aSurfaces.push_back(
        Surface::makeShared<CylinderSurface>(m_transform,
                                             m_bounds->r() - 0.5 * thickness(),
                                             m_bounds->halflengthZ()));
    aSurfaces.push_back(
        Surface::makeShared<CylinderSurface>(m_transform,
                                             m_bounds->r() + 0.5 * thickness(),
                                             m_bounds->halflengthZ()));
    // create an ApproachDescriptor with standard surfaces surfaces - these will
    // be deleted by the approach descriptor
    m_approachDescriptor = std::make_unique<const GenericApproachDescriptor>(
        std::move(aSurfaces));
  }
  for (auto& sfPtr : (m_approachDescriptor->containedSurfaces())) {
    if (sfPtr != nullptr) {
      auto& mutableSf = *(const_cast<Surface*>(sfPtr));
      mutableSf.associateLayer(*this);
    }
  }
}
