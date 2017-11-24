// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// DiscLayer.cpp, ACTS project
///////////////////////////////////////////////////////////////////

#include "ACTS/Layers/DiscLayer.hpp"
#include "ACTS/Layers/GenericApproachDescriptor.hpp"
#include "ACTS/Material/SurfaceMaterial.hpp"
#include "ACTS/Surfaces/DiscBounds.hpp"
#include "ACTS/Surfaces/RadialBounds.hpp"
#include "ACTS/Utilities/BinUtility.hpp"
#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Volumes/AbstractVolume.hpp"
#include "ACTS/Volumes/BoundarySurfaceFace.hpp"
#include "ACTS/Volumes/CylinderVolumeBounds.hpp"

Acts::DiscLayer::DiscLayer(std::shared_ptr<const Transform3D>  transform,
                           std::shared_ptr<const DiscBounds>   dbounds,
                           std::unique_ptr<SurfaceArray_old>       surfaceArray,
                           double                              thickness,
                           std::unique_ptr<ApproachDescriptor> ades,
                           LayerType                           laytyp)
  : DiscSurface(transform, dbounds)
  , Layer(std::move(surfaceArray), thickness, std::move(ades), laytyp)
{
  // create the representing volume
  const RadialBounds* rBounds
      = dynamic_cast<const RadialBounds*>(dbounds.get());
  if (rBounds) {
    // @todo make a trapezoidal volume when you have DiscTrapezoidalBounds
    CylinderVolumeBounds* cvBounds = new CylinderVolumeBounds(
        rBounds->rMin(), rBounds->rMax(), 0.5 * thickness);
    Layer::m_representingVolume
        = new AbstractVolume(transform, VolumeBoundsPtr(cvBounds));
  }
  // associate the layer to the layer surface itself
  DiscSurface::associateLayer(*this);
  // build an approach descriptor if none provided
  if (!m_approachDescriptor && Layer::m_surfaceArray) buildApproachDescriptor();
  // register the layer to the approach descriptor
  if (m_approachDescriptor) approachDescriptor()->registerLayer(*this);
}

const Acts::DiscSurface&
Acts::DiscLayer::surfaceRepresentation() const
{
  return (*this);
}

Acts::DiscSurface&
Acts::DiscLayer::surfaceRepresentation()
{
  return (*this);
}

void
Acts::DiscLayer::buildApproachDescriptor()
{
  // delete it
  m_approachDescriptor = nullptr;
  // take the boundary surfaces of the representving volume if they exist
  if (m_representingVolume) {
    // get teh boundary surfaces
    const std::vector<std::shared_ptr<const BoundarySurfaceT<AbstractVolume>>>&
        bSurfaces
        = m_representingVolume->boundarySurfaces();
    // fill in the surfaces into the vector
    std::vector<std::shared_ptr<const BoundarySurfaceT<AbstractVolume>>>
        aSurfaces;
    aSurfaces.push_back(bSurfaces.at(negativeFaceXY));
    aSurfaces.push_back(bSurfaces.at(positiveFaceXY));
    // create an ApproachDescriptor with Boundary surfaces
    m_approachDescriptor = std::
        make_unique<const GenericApproachDescriptor<BoundarySurfaceT<AbstractVolume>>>(
            aSurfaces);
  } else {
    // create the new surfaces - positions first
    Vector3D aspPosition(center() + 0.5 * thickness() * normal());
    Vector3D asnPosition(center() - 0.5 * thickness() * normal());
    auto     asnTransform
        = std::make_shared<const Transform3D>(Translation3D(asnPosition));
    auto aspTransform
        = std::make_shared<const Transform3D>(Translation3D(aspPosition));
    // create the vector
    std::vector<const Surface*> aSurfaces;
    aSurfaces.push_back(new DiscSurface(asnTransform, m_bounds));
    aSurfaces.push_back(new DiscSurface(aspTransform, m_bounds));
    // create an ApproachDescriptor with standard surfaces surfaces - these will
    // be deleted by the approach descriptor
    m_approachDescriptor
        = std::make_unique<const GenericApproachDescriptor<Surface>>(aSurfaces);
  }
  // @todo check if we can give the layer at curface creation
  for (auto& sfPtr : (m_approachDescriptor->containedSurfaces())) {
    if (sfPtr) {
      auto& mutableSf = *(const_cast<Surface*>(sfPtr));
      mutableSf.associateLayer(*this);
    }
  }
}
