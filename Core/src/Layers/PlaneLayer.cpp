// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// PlaneLayer.cpp, ACTS project
///////////////////////////////////////////////////////////////////

// Geometry module
#include "ACTS/Layers/PlaneLayer.hpp"
#include "ACTS/EventData/TrackParameters.hpp"
#include "ACTS/Layers/GenericApproachDescriptor.hpp"
#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Utilities/Helpers.hpp"

Acts::PlaneLayer::PlaneLayer(std::shared_ptr<Transform3D>         transform,
                             std::shared_ptr<const PlanarBounds>& pbounds,
                             std::unique_ptr<SurfaceArray>        surfaceArray,
                             double                               thickness,
                             std::unique_ptr<ApproachDescriptor>  ades,
                             LayerType                            laytyp)
  : PlaneSurface(transform, pbounds)
  , Layer(std::move(surfaceArray), thickness, std::move(ades), laytyp)
{
  // @todo create representing volume
  // register the layer to the surface
  Acts::PlaneSurface::associateLayer(*this);
  // deal with the approach descriptor
  if (!m_approachDescriptor && surfaceArray) buildApproachDescriptor();
  // register the layer if the approach descriptor was provided
  if (m_approachDescriptor) m_approachDescriptor->registerLayer(*this);
  // set the material if present
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

Acts::PlaneLayer::PlaneLayer(const PlaneLayer& play, const Transform3D& transf)
  : PlaneSurface(play, transf), Layer(play)
{
}

const Acts::PlaneSurface&
Acts::PlaneLayer::surfaceRepresentation() const
{
  return (*this);
}

Acts::PlaneSurface&
Acts::PlaneLayer::surfaceRepresentation()
{
  return (*this);
}

void
Acts::PlaneLayer::buildApproachDescriptor()
{
  // delete it
  m_approachDescriptor = nullptr;
  // delete the surfaces
  std::vector<const Acts::Surface*> aSurfaces;
  // get the appropriate transform, the center and the normal vector
  const Transform3D& lTransform = PlaneSurface::transform();
  RotationMatrix3D   lRotation  = lTransform.rotation();
  const Vector3D&    lCenter    = PlaneSurface::center();
  const Vector3D&    lVector    = PlaneSurface::normal();
  // create new surfaces
  Transform3D* apnTransform = new Transform3D(getTransformFromRotTransl(
      lRotation, (lCenter - 0.5 * Layer::m_layerThickness * lVector)));
  Transform3D* appTransform = new Transform3D(getTransformFromRotTransl(
      lRotation, (lCenter + 0.5 * Layer::m_layerThickness * lVector)));
  // create the new surfaces
  aSurfaces.push_back(new Acts::PlaneSurface(
      std::shared_ptr<Transform3D>(apnTransform), PlaneSurface::m_bounds));
  aSurfaces.push_back(new PlaneSurface(
      std::shared_ptr<Transform3D>(appTransform), PlaneSurface::m_bounds));
  // set the layer and make TrackingGeometry
  for (auto& sfPtr : aSurfaces) {
    auto& mutableSf = *(const_cast<Surface*>(sfPtr));
    mutableSf.associateLayer(*this);
  }
  m_approachDescriptor
      = std::make_unique<GenericApproachDescriptor<const Surface>>(aSurfaces);
}
