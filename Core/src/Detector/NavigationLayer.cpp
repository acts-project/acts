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

// Geometry module
#include "ACTS/Layers/NavigationLayer.hpp"
#include "ACTS/Surfaces/Surface.hpp"

// constructor with arguments
Acts::NavigationLayer::NavigationLayer(Acts::Surface* surfaceRepresentation, double thickness):
 Acts::Layer(),                                     
 m_surfaceRepresentation(surfaceRepresentation)
{
  Layer::m_layerThickness = thickness;
  // @TODO temporary - until GeoID service is in place
  assignGeoID(GeometryID(0));  
}


// constructor with shift
std::shared_ptr<const Acts::Layer> Acts::NavigationLayer::cloneWithShift(const Acts::Transform3D& shift) const
{
  Acts::Surface* shiftedSurface = m_surfaceRepresentation->clone(&shift);
  return std::shared_ptr<const Acts::Layer>( new Acts::NavigationLayer(shiftedSurface,Layer::m_layerThickness));
}

// destructor - only deletes surface representation
Acts::NavigationLayer::~NavigationLayer()
{ delete m_surfaceRepresentation; } 

bool Acts::NavigationLayer::isOnLayer(const Acts::Vector3D& gp, const BoundaryCheck& bcheck) const {
  return m_surfaceRepresentation->isOnSurface(gp, bcheck);
}
