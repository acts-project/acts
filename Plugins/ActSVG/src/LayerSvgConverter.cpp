// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/ActSVG/LayerSvgConverter.hpp"

#include "Acts/Geometry/Layer.hpp"
#include "Acts/Plugins/ActSVG/SurfaceArraySvgConverter.hpp"
#include "Acts/Plugins/ActSVG/SurfaceSvgConverter.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <set>
#include <sstream>

std::vector<actsvg::svg::object> Acts::Svg::layerSheets(
    const GeometryContext& gctx, const Layer& layer,
    const std::string& layerName, const Style& surfaceStyle) {
  // The local logger
  ACTS_LOCAL_LOGGER(
      getDefaultLogger("SurfaceArraySvgConverter", Logging::INFO));

  // The sheets
  std::vector<actsvg::svg::object> sheets;

  // The volume
  Acts::Svg::ProtoVolume volume;
  volume._name = layerName;
  ACTS_INFO("Processing layer: " << layerName);

  /// Convert the surface array into proto surfaces and a grid structure
  if (layer.surfaceArray() != nullptr) {
    auto [surfaces, grid] =
        convert(gctx, *(layer.surfaceArray()), surfaceStyle);
    volume._surfaces = surfaces;
    volume._surface_grid = grid;
  }

  // The sheet
  actsvg::svg::object module_sheet;
  actsvg::svg::object grid_sheet;

  const auto& layerSurface = layer.surfaceRepresentation();
  if (layerSurface.type() == Acts::Surface::Disc) {
    module_sheet = actsvg::display::endcap_sheet(
        layerName + "_modules", volume, {800, 800},
        actsvg::display::e_module_info);
    grid_sheet = actsvg::display::endcap_sheet(
        layerName + "_grid", volume, {800, 800}, actsvg::display::e_grid_info);
  } else if (layerSurface.type() == Acts::Surface::Cylinder){
    module_sheet = actsvg::display::barrel_sheet(
        layerName + "_modules", volume, {800, 800},
        actsvg::display::e_module_info);
    grid_sheet = actsvg::display::barrel_sheet(
        layerName + "_grid", volume, {800, 800}, actsvg::display::e_grid_info);
  }

  // Register
  sheets.push_back(module_sheet);
  sheets.push_back(grid_sheet);

  return sheets;
}