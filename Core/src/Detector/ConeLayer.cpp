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

// Geometry module
#include "ACTS/Layers/ConeLayer.hpp"

#include "ACTS/EventData/TrackParameters.hpp"
#include "ACTS/Surfaces/ConeBounds.hpp"
#include "ACTS/Utilities/Definitions.hpp"
// Core module

Acts::ConeLayer::ConeLayer(std::shared_ptr<Acts::Transform3D>      transform,
                           std::shared_ptr<const Acts::ConeBounds> cbounds,
                           std::unique_ptr<SurfaceArray>           surfaceArray,
                           double                                  thickness,
                           Acts::OverlapDescriptor*                olap,
                           int                                     laytyp)
  : ConeSurface(transform, cbounds)
  , Layer(std::move(surfaceArray), thickness, olap, nullptr, laytyp)
{
}

Acts::ConeLayer::ConeLayer(const Acts::ConeLayer&   clay,
                           const Acts::Transform3D& transf)
  : ConeSurface(clay, transf), Layer(clay)
{
}

const Acts::ConeSurface&
Acts::ConeLayer::surfaceRepresentation() const
{
  return (*this);
}
