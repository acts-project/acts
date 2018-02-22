// This file is part of the ACTS project.
//
// Copyright (C) 2016-2017 ACTS project team
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
  : Acts::Layer(nullptr)
  , m_surfaceRepresentation(std::move(surfaceRepresentation))
{
  Layer::m_layerThickness = thickness;
  Layer::m_layerType      = navigation;
}

Acts::NavigationLayer::~NavigationLayer()
{
}

bool
Acts::NavigationLayer::resolve(bool, bool, bool) const
{
  return false;
}
