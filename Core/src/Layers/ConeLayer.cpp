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

Acts::ConeLayer::ConeLayer(std::shared_ptr<const Transform3D>  transform,
                           std::shared_ptr<const ConeBounds>   cbounds,
                           std::unique_ptr<SurfaceArray_old>       surfaceArray,
                           double                              thickness,
                           std::unique_ptr<ApproachDescriptor> ade,
                           LayerType                           laytyp)
  : ConeSurface(transform, cbounds)
  , Layer(std::move(surfaceArray), thickness, std::move(ade), laytyp)
{
}

const Acts::ConeSurface&
Acts::ConeLayer::surfaceRepresentation() const
{
  return (*this);
}

Acts::ConeSurface&
Acts::ConeLayer::surfaceRepresentation()
{
  return (*this);
}
