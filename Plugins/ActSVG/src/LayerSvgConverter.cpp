// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/ActSVG/LayerSvgConverter.hpp"

#include "Acts/Geometry/Layer.hpp"
#include "Acts/Surfaces/Surface.hpp"

std::vector<actsvg::svg::object> Acts::Svg::layerSheets(
    const GeometryContext& gctx, const Layer& layer, const Style& surfaceStyle,
    const std::string& identification) {
  std::vector<actsvg::svg::object> sheets;

  if (layer.surfaceArray() != nullptr) {
    const auto& surfaces = layer.surfaceArray()->surfaces();
    // x-y view is needed 
    if (layer.surfaceRepresentation().type() ==
            Acts::Surface::SurfaceType::Disc or
        layer.surfaceRepresentation().type() ==
            Acts::Surface::SurfaceType::Plane) {




    }
  }

  return sheets;
}