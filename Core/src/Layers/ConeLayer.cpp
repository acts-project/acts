// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// ConeLayer.cpp, ACTS project
///////////////////////////////////////////////////////////////////

#include "ACTS/Layers/ConeLayer.hpp"
#include "ACTS/EventData/TrackParameters.hpp"
#include "ACTS/Surfaces/ConeBounds.hpp"
#include "ACTS/Utilities/Definitions.hpp"

Acts::ConeLayer::ConeLayer(std::shared_ptr<Transform3D>      transform,
                           std::shared_ptr<const ConeBounds> cbounds,
                           std::unique_ptr<SurfaceArray>     surfaceArray,
                           double                            thickness,
                           OverlapDescriptor*                olap,
                           ApproachDescriptor*               ade,
                           LayerType                         laytyp)
  : ConeSurface(transform, cbounds)
  , Layer(std::move(surfaceArray), thickness, olap, ade, laytyp)
{
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

Acts::ConeLayer::ConeLayer(const ConeLayer& clay, const Transform3D& transf)
  : ConeSurface(clay, transf), Layer(clay)
{
}

const Acts::ConeSurface&
Acts::ConeLayer::surfaceRepresentation() const
{
  return (*this);
}
