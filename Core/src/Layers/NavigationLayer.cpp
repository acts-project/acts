// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// NavigationLayer.cpp, ACTS project
///////////////////////////////////////////////////////////////////

#include "ACTS/Layers/NavigationLayer.hpp"
#include "ACTS/Surfaces/Surface.hpp"

Acts::NavigationLayer::NavigationLayer(
    std::unique_ptr<const Surface> surfaceRepresentation,
    double                         thickness)
  : Acts::Layer(), m_surfaceRepresentation(std::move(surfaceRepresentation))
{
  Layer::m_layerThickness = thickness;
  Layer::m_layerType      = navigation;
}

std::shared_ptr<const Acts::Layer>
Acts::NavigationLayer::cloneWithShift(const Acts::Transform3D& shift) const
{
  Surface* shiftedSurface = m_surfaceRepresentation->clone(&shift);
  return std::shared_ptr<const Acts::Layer>(
      new Acts::NavigationLayer(std::unique_ptr<Surface>(shiftedSurface), Layer::m_layerThickness));
}

Acts::NavigationLayer::~NavigationLayer()
{
}

bool
Acts::NavigationLayer::isOnLayer(const Vector3D& gp,
                                 const BoundaryCheck&  bcheck) const
{
  return m_surfaceRepresentation->isOnSurface(gp, bcheck);
}
