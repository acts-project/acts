// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// CylinderLayer.cpp, ACTS project
///////////////////////////////////////////////////////////////////

#include "ACTS/Layers/CylinderLayer.hpp"
#include "ACTS/Layers/GenericApproachDescriptor.hpp"
#include "ACTS/Material/SurfaceMaterial.hpp"
#include "ACTS/Surfaces/CylinderBounds.hpp"
#include "ACTS/Utilities/BinUtility.hpp"
#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Volumes/AbstractVolume.hpp"
#include "ACTS/Volumes/BoundarySurfaceFace.hpp"
#include "ACTS/Volumes/CylinderVolumeBounds.hpp"

Acts::CylinderLayer::CylinderLayer(
    std::shared_ptr<Transform3D>          transform,
    std::shared_ptr<const CylinderBounds> cBounds,
    std::unique_ptr<SurfaceArray>         surfaceArray,
    double                                thickness,
    ApproachDescriptor*                   ades,
    LayerType                             laytyp)
  : CylinderSurface(transform, cBounds)
  , Layer(std::move(surfaceArray), thickness, ades, laytyp)
{
  // create the representing volume
  CylinderVolumeBounds* cvBounds
      = new CylinderVolumeBounds(cBounds->r() - 0.5 * thickness,
                                 cBounds->r() + 0.5 * thickness,
                                 cBounds->halflengthZ());
  Layer::m_representingVolume
      = new AbstractVolume(transform, VolumeBoundsPtr(cvBounds));
  // associate teh layer
  CylinderSurface::associateLayer(*this);
  // an approach descriptor is automatically created if there's a surface array
  if (!ades && Layer::m_surfaceArray) buildApproachDescriptor();
  // register the layer if the approach descriptor was provided
  if (ades) m_approachDescriptor->registerLayer(*this);
  // set the material surface if present
  // material can be on any approach surface or on the representing surface
  if (m_approachDescriptor) {
    // the approach surfaces
    const std::vector<const Surface*>& approachSurfaces
        = m_approachDescriptor->containedSurfaces();
    for (auto& aSurface : approachSurfaces)
      if (aSurface->associatedMaterial()) m_materialSurface = aSurface;
  }
  if (surfaceRepresentation().associatedMaterial())
    m_materialSurface = &surfaceRepresentation();
}

Acts::CylinderLayer::CylinderLayer(const CylinderLayer& clay,
                                   const Transform3D&   transf)
  : CylinderSurface(clay, transf), Layer(clay)
{
  if (m_surfaceArray) buildApproachDescriptor();
}

const Acts::CylinderSurface&
Acts::CylinderLayer::surfaceRepresentation() const
{
  return (*this);
}

void
Acts::CylinderLayer::buildApproachDescriptor() const
{
  // delete it
  delete m_approachDescriptor;
  // delete the surfaces
  // take the boundary surfaces of the representving volume if they exist
  if (m_representingVolume) {
    // get teh boundary surfaces
    const std::
        vector<std::shared_ptr<const BoundarySurfaceT<AbstractVolume>>>&
            bSurfaces  = m_representingVolume->boundarySurfaces();
    // fill in the surfaces into the vector
    std::vector<std::shared_ptr<const BoundarySurfaceT<AbstractVolume>>>
        aSurfaces;
    if (bSurfaces.size() > size_t(tubeOuterCover))
      aSurfaces.push_back(bSurfaces.at(tubeInnerCover));
    aSurfaces.push_back(bSurfaces.at(tubeOuterCover));
    // create an ApproachDescriptor with Boundary surfaces
    m_approachDescriptor = new 
        GenericApproachDescriptor<const BoundarySurfaceT<AbstractVolume>>(
            aSurfaces);
  } else {
    // create the new surfaces
    std::vector<const Acts::Surface*> aSurfaces;
    aSurfaces.push_back(
        new CylinderSurface(m_transform,
                            m_bounds->r() - 0.5 * thickness(),
                            m_bounds->halflengthZ()));
    aSurfaces.push_back(
        new CylinderSurface(m_transform,
                            m_bounds->r() + 0.5 * thickness(),
                            m_bounds->halflengthZ()));
    // create an ApproachDescriptor with standard surfaces surfaces - these will
    // be deleted by the approach descriptor
    m_approachDescriptor
        = new GenericApproachDescriptor<const Surface>(aSurfaces);
  }
  for (auto& sIter : (m_approachDescriptor->containedSurfaces())) {
    if (sIter) sIter->associateLayer(*this);
  }
}
