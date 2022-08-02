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
    const std::string& layerName, const Style& surfaceStyle,
    const std::array<std::array<ActsScalar, 2>, 2>& rangeRes) {
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
  actsvg::svg::object xy_layer;
  actsvg::svg::object zr_layer;

  const auto& layerSurface = layer.surfaceRepresentation();
  if (layerSurface.type() == Acts::Surface::Disc) {
    module_sheet = actsvg::display::endcap_sheet(
        layerName + "_modules", volume, {800, 800},
        actsvg::display::e_module_info);
    grid_sheet = actsvg::display::endcap_sheet(
        layerName + "_grid", volume, {800, 800}, actsvg::display::e_grid_info);
  } else if (layerSurface.type() == Acts::Surface::Cylinder) {
    // Draw the barrel type
    module_sheet = actsvg::display::barrel_sheet(
        layerName + "_modules", volume, {800, 800},
        actsvg::display::e_module_info);
    grid_sheet = actsvg::display::barrel_sheet(
        layerName + "_grid", volume, {800, 800}, actsvg::display::e_grid_info);
  }

  // The z_r view of things
  actsvg::views::z_r z_r_view;
  actsvg::views::x_y x_y_view;

  if (layer.surfaceArray() != nullptr) {
    // The x_y view of things
    xy_layer._tag = "g";
    xy_layer._id = layerName + "_xy_view";
    // The x_r view of things
    zr_layer._tag = "g";
    zr_layer._id = layerName + "_zr_view";
    unsigned int m = 0;
    // Potential labels
    Acts::ActsScalar avgRadius = 0.;
    Acts::ActsScalar avgZ = 0.;

    for (const auto& sf : layer.surfaceArray()->surfaces()) {
      // Surface center
      const Acts::Vector3 rCenter = sf->binningPosition(gctx, Acts::binR);
      const Acts::Vector3 sfCenter = sf->center(gctx);
      Acts::ActsScalar radius = Acts::VectorHelpers::perp(rCenter);
      Acts::ActsScalar phi = Acts::VectorHelpers::phi(rCenter);
      Acts::ActsScalar z = sfCenter.z();
      // Raw display surfaces for projects
      actsvg::proto::surface<std::vector<Acts::Vector3>> projSurface;
      projSurface._vertices = sf->polyhedronRepresentation(gctx, 1u).vertices;
      // Draw only if they fall into the range restriction - for phi
      if (phi >= rangeRes[1][0] and phi <= rangeRes[1][1]) {
        std::string m_zr_id = std::string("zr_") + std::to_string(m++);
        zr_layer.add_object(Acts::Svg::surfaceViewZR(projSurface, m_zr_id));
      } 
      // for z 
      if (z >= rangeRes[0][0] and z <= rangeRes[0][1]) {
        std::string m_xy_id = std::string("xy_") + std::to_string(m++);
        xy_layer.add_object(Acts::Svg::surfaceViewXY(projSurface, m_xy_id));
      }
    }
  }

  // Register
  sheets.push_back(module_sheet);
  sheets.push_back(grid_sheet);
  sheets.push_back(xy_layer);
  sheets.push_back(zr_layer);

  return sheets;
}
