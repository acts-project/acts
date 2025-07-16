// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/NavigationLayer.hpp"

namespace Acts {

NavigationLayer::NavigationLayer(
    std::shared_ptr<const Surface> surfaceRepresentation, double thickness)
    : Layer(nullptr),
      m_surfaceRepresentation(std::move(surfaceRepresentation)) {
  m_layerThickness = thickness;
  m_layerType = navigation;
}

NavigationLayer::~NavigationLayer() = default;

}  // namespace Acts
