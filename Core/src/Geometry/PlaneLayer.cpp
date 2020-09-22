// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// Geometry module
#include "Acts/Geometry/PlaneLayer.hpp"

#include "Acts/Geometry/GenericApproachDescriptor.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Helpers.hpp"

#include <utility>

Acts::PlaneLayer::PlaneLayer(const Transform3D& transform,
                             std::shared_ptr<const PlanarBounds>& pbounds,
                             std::unique_ptr<SurfaceArray> surfaceArray,
                             double thickness,
                             std::unique_ptr<ApproachDescriptor> ades,
                             LayerType laytyp)
    : PlaneSurface(transform, pbounds),
      Layer(std::move(surfaceArray), thickness, std::move(ades), laytyp) {
  // @todo create representing volume
  // register the layer to the surface
  Acts::PlaneSurface::associateLayer(*this);
  // deal with the approach descriptor
  if (!m_approachDescriptor && m_surfaceArray) {
    buildApproachDescriptor();
  }
  // register the layer to the approach descriptor
  if (m_approachDescriptor) {
    approachDescriptor()->registerLayer(*this);
  }
}

const Acts::PlaneSurface& Acts::PlaneLayer::surfaceRepresentation() const {
  return (*this);
}

Acts::PlaneSurface& Acts::PlaneLayer::surfaceRepresentation() {
  return (*this);
}

void Acts::PlaneLayer::buildApproachDescriptor() {
  // delete it
  m_approachDescriptor.reset(nullptr);
  // delete the surfaces
  std::vector<std::shared_ptr<const Acts::Surface>> aSurfaces;
  // get the appropriate transform, the center and the normal vector

  //@todo fix with representing volume
  const Transform3D& lTransform = PlaneSurface::transform(GeometryContext());
  RotationMatrix3D lRotation = lTransform.rotation();
  const Vector3D& lCenter = PlaneSurface::center(GeometryContext());
  const Vector3D& lVector = Surface::normal(GeometryContext(), lCenter);
  // create new surfaces
  const Transform3D apnTransform = Transform3D(
      Translation3D(lCenter - 0.5 * Layer::m_layerThickness * lVector) *
      lRotation);
  const Transform3D appTransform = Transform3D(
      Translation3D(lCenter + 0.5 * Layer::m_layerThickness * lVector) *
      lRotation);
  // create the new surfaces
  aSurfaces.push_back(Surface::makeShared<Acts::PlaneSurface>(
      apnTransform, PlaneSurface::m_bounds));
  aSurfaces.push_back(Surface::makeShared<Acts::PlaneSurface>(
      appTransform, PlaneSurface::m_bounds));
  // set the layer and make TrackingGeometry
  for (auto& sfPtr : aSurfaces) {
    auto mutableSf = const_cast<Surface*>(sfPtr.get());
    mutableSf->associateLayer(*this);
  }
  // @todo check if we can provide the layer at surface creation
  m_approachDescriptor =
      std::make_unique<const GenericApproachDescriptor>(std::move(aSurfaces));
}