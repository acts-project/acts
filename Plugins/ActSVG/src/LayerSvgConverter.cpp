// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/ActSVG/LayerSvgConverter.hpp"

#include "Acts/Geometry/Layer.hpp"
#include "Acts/Surfaces/SurfaceArray.hpp"
#include "ActsPlugins/ActSVG/SurfaceArraySvgConverter.hpp"
#include "ActsPlugins/ActSVG/SurfaceSvgConverter.hpp"

#include <stdexcept>

using namespace Acts;

std::vector<actsvg::svg::object> ActsPlugins::Svg::LayerConverter::convert(
    const GeometryContext& gctx, const Layer& layer,
    const LayerConverter::Options& cOptions) {
  // The sheets
  std::vector<actsvg::svg::object> sheets;

  // The volume
  ActsPlugins::Svg::ProtoVolume volume;
  volume._name = cOptions.name;

  /// Convert the surface array into proto surfaces and a grid structure
  if (layer.surfaceArray() != nullptr) {
    SurfaceArrayConverter::Options sacOptions;
    sacOptions.surfaceStyles = cOptions.surfaceStyles;
    auto [surfaces, grid, associations] = SurfaceArrayConverter::convert(
        gctx, *(layer.surfaceArray()), sacOptions);
    volume._surfaces = {surfaces};
    volume._surface_grid = grid;
    volume._grid_associations = {associations};
  }

  // The sheet. module_sheet and grid_sheet stay undefined: actsvg 0.4.57
  // removed the displays that filled them, and every consumer already skips
  // undefined sheets, which is what the moduleInfo/gridInfo options produced
  // when switched off.
  actsvg::svg::object module_sheet;
  actsvg::svg::object grid_sheet;
  actsvg::svg::object xy_layer;
  actsvg::svg::object zr_layer;

  if (cOptions.moduleInfo || cOptions.gridInfo) {
    throw std::invalid_argument(
        "LayerConverter: the module and grid sheets were removed in actsvg "
        "0.4.57 and cannot be produced");
  }

  // The z_r view of things
  actsvg::views::z_r z_r_view;
  actsvg::views::x_y x_y_view;

  if (layer.surfaceArray() != nullptr) {
    // The x_y view of things
    xy_layer._tag = "g";
    xy_layer._id = cOptions.name + "_xy_view";
    // The x_r view of things
    zr_layer._tag = "g";
    zr_layer._id = cOptions.name + "_zr_view";
    unsigned int m = 0;
    // Potential labels
    double avgRadius = 0.;

    for (const auto& sf : layer.surfaceArray()->surfaces()) {
      // Surface center
      const Vector3 rCenter = sf->referencePosition(gctx, AxisDirection::AxisR);
      const Vector3 sfCenter = sf->center(gctx);
      double radius = VectorHelpers::perp(rCenter);
      double phi = VectorHelpers::phi(rCenter);
      double z = sfCenter.z();
      // Get the average radius
      avgRadius += radius;
      // Raw display surfaces for projections
      actsvg::proto::surface<std::vector<Vector3>> projSurface;
      projSurface._vertices = sf->polyhedronRepresentation(gctx, 1u).vertices;
      // Draw only if they fall into the range restriction - for phi
      if (phi >= cOptions.phiRange[0] && phi <= cOptions.phiRange[1]) {
        std::string m_zr_id = std::string("zr_") + std::to_string(m++);
        zr_layer.add_object(ActsPlugins::Svg::View::zr(projSurface, m_zr_id));
      }
      // for z
      if (z >= cOptions.zRange[0] && z <= cOptions.zRange[1]) {
        std::string m_xy_id = std::string("xy_") + std::to_string(m++);
        xy_layer.add_object(ActsPlugins::Svg::View::xy(projSurface, m_xy_id));
      }
    }
    // Do the average
    avgRadius /= layer.surfaceArray()->surfaces().size();

    // Add a measure iuf requested
    if (cOptions.labelProjection) {
      double xEnd = avgRadius * std::cos(cOptions.labelGauge);
      double yEnd = avgRadius * std::sin(cOptions.labelGauge);
      xy_layer.add_object(measure(0., 0., xEnd, yEnd, "r", avgRadius, "mm"));
    }
  }

  // Register according to the enums
  sheets.push_back(module_sheet);
  sheets.push_back(grid_sheet);
  sheets.push_back(xy_layer);
  sheets.push_back(zr_layer);

  // Return the created sheets
  return sheets;
}
