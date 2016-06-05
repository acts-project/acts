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

// Geometry module
#include "ACTS/Layers/DiscLayer.hpp"

#include "ACTS/Detector/GenericApproachDescriptor.hpp"
#include "ACTS/Detector/GenericOverlapDescriptor.hpp"
#include "ACTS/Material/SurfaceMaterial.hpp"
#include "ACTS/Surfaces/DiscBounds.hpp"
#include "ACTS/Surfaces/RadialBounds.hpp"
#include "ACTS/Utilities/BinUtility.hpp"
#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Volumes/AbstractVolume.hpp"
#include "ACTS/Volumes/BoundarySurfaceFace.hpp"
#include "ACTS/Volumes/CylinderVolumeBounds.hpp"
// Core module

Acts::DiscLayer::DiscLayer(std::shared_ptr<Acts::Transform3D>      transform,
                           std::shared_ptr<const Acts::DiscBounds> dbounds,
                           std::unique_ptr<SurfaceArray>           surfaceArray,
                           double                                  thickness,
                           Acts::OverlapDescriptor*                olap,
                           Acts::ApproachDescriptor*               ades,
                           int                                     laytyp)
  : DiscSurface(transform, dbounds)
  , Layer(std::move(surfaceArray), thickness, olap, ades, laytyp)
{
  // just create a generic overlap descriptor if none is there
  if (!Layer::m_overlapDescriptor)
    Layer::m_overlapDescriptor = new GenericOverlapDescriptor();
  // create the representing volume
  const Acts::RadialBounds* rBounds
      = dynamic_cast<const Acts::RadialBounds*>(dbounds.get());
  if (rBounds) {
    // @TODO make a trapezoidal volume when you have DiscTrapezoidalBounds
    CylinderVolumeBounds* cvBounds = new CylinderVolumeBounds(
        rBounds->rMin(), rBounds->rMax(), 0.5 * thickness);
    Layer::m_representingVolume
        = new AbstractVolume(transform, VolumeBoundsPtr(cvBounds));
  }
  // associate teh layer to this
  DiscSurface::associateLayer(*this);
  if (!ades && Layer::m_surfaceArray) buildApproachDescriptor();
  // register the layer
  if (ades) m_approachDescriptor->registerLayer(*this);
}

Acts::DiscLayer::DiscLayer(const Acts::DiscLayer&   dlay,
                           const Acts::Transform3D& transf)
  : DiscSurface(dlay, transf), Layer(dlay)
{
  DiscSurface::associateLayer(*this);
  if (m_surfaceArray) buildApproachDescriptor();
}

const Acts::DiscSurface&
Acts::DiscLayer::surfaceRepresentation() const
{
  return (*this);
}

/** build approach surfaces */
void
Acts::DiscLayer::buildApproachDescriptor() const
{
  // delete it
  delete m_approachDescriptor;
  // delete the surfaces
  // take the boundary surfaces of the representving volume if they exist
  if (m_representingVolume) {
    // get teh boundary surfaces
    const std::
        vector<std::shared_ptr<const Acts::
                                   BoundarySurface<Acts::AbstractVolume>>>&
            bSurfaces
        = m_representingVolume->boundarySurfaces();
    // fill in the surfaces into the vector
    std::vector<std::shared_ptr<const Acts::
                                    BoundarySurface<Acts::AbstractVolume>>>
        aSurfaces;
    aSurfaces.push_back(bSurfaces.at(negativeFaceXY));
    aSurfaces.push_back(bSurfaces.at(positiveFaceXY));
    // create an ApproachDescriptor with Boundary surfaces
    m_approachDescriptor = new Acts::
        GenericApproachDescriptor<const BoundarySurface<AbstractVolume>>(
            aSurfaces);
  } else {
    // create the new surfaces - positions first
    Acts::Vector3D     aspPosition(center() + 0.5 * thickness() * normal());
    Acts::Vector3D     asnPosition(center() - 0.5 * thickness() * normal());
    Acts::Transform3D* asnTransform
        = new Acts::Transform3D(Acts::Translation3D(asnPosition));
    Acts::Transform3D* aspTransform
        = new Acts::Transform3D(Acts::Translation3D(aspPosition));
    // create the vector
    std::vector<const Acts::Surface*> aSurfaces;
    aSurfaces.push_back(new Acts::DiscSurface(
        std::shared_ptr<Acts::Transform3D>(asnTransform), m_bounds));
    aSurfaces.push_back(new Acts::DiscSurface(
        std::shared_ptr<Acts::Transform3D>(aspTransform), m_bounds));
    // create an ApproachDescriptor with standard surfaces surfaces - these will
    // be deleted by the approach descriptor
    m_approachDescriptor
        = new Acts::GenericApproachDescriptor<const Acts::Surface>(aSurfaces);
  }
  for (auto& sIter : (m_approachDescriptor->containedSurfaces())) {
    sIter->associateLayer(*this);
  }
}
