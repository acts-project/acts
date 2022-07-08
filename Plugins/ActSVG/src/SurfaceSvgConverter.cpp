// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/ActSVG/SurfaceSvgConverter.hpp"

#include "Acts/Surfaces/PlanarBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"

Acts::Svg::ProtoSurface Acts::Svg::convert(const GeometryContext& gctx,
                                           const Surface& surface,
                                           unsigned int nSegs,
                                           const Style& surfaceStyle,
                                           bool templateSurface) {
  ProtoSurface pSurface;

  // Planar surfaces to be draws as polyhedron representations
  if (surface.type() == Acts::Surface::SurfaceType::Plane) {
    // set the type
    pSurface._type = ProtoSurface::type::e_rectangle;
    // This is not a template, hance properly placed
    if (not templateSurface) {
      Polyhedron surfaceHedron = surface.polyhedronRepresentation(gctx, nSegs);
      auto vertices3D = surfaceHedron.vertices;
      pSurface._vertices = vertices3D;
    } else {
      auto planarBounds =
          dynamic_cast<const Acts::PlanarBounds*>(&(surface.bounds()));
      if (planarBounds != nullptr) {
        auto vertices2D = planarBounds->vertices(nSegs);
        pSurface._vertices.reserve(vertices2D.size());
        for (const auto& v2 : vertices2D) {
          pSurface._vertices.push_back({v2[0], v2[1], 0.});
        }
      }
    }
    // Translate the bound values
    const auto& boundValues = surface.bounds().values();
    if (surface.bounds().type() ==
        Acts::SurfaceBounds::BoundsType::eRectangle) {
      pSurface._measures = {
          static_cast<actsvg::scalar>(0.5 * (boundValues[2] - boundValues[0])),
          static_cast<actsvg::scalar>(0.5 * (boundValues[3] - boundValues[1]))};
    }

  } else if (surface.type() == Acts::Surface::SurfaceType::Disc) {
  }

  // Attach the style
  pSurface._fill._fc = {surfaceStyle.fillColor,
                        static_cast<actsvg::scalar>(surfaceStyle.fillOpacity)};

  // Fill style
  pSurface._fill._fc._hl_rgb = surfaceStyle.highlightColor;
  pSurface._fill._fc._highlight = surfaceStyle.highlights;

  // Stroke sytle
  pSurface._stroke._sc = actsvg::style::color{surfaceStyle.strokeColor};
  pSurface._stroke._width = surfaceStyle.strokeWidth;

  return pSurface;
}