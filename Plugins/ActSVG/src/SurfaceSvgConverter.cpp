// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/ActSVG/SurfaceSvgConverter.hpp"

#include "Acts/Surfaces/AnnulusBounds.hpp"
#include "Acts/Surfaces/PlanarBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"

Acts::Svg::ProtoSurface Acts::Svg::convert(const GeometryContext& gctx,
                                           const Surface& surface,
                                           const Style& surfaceStyle,
                                           bool templateSurface) {
  ProtoSurface pSurface;

  // Planar surfaces to be draws as polyhedron representations
  if (surface.type() == Acts::Surface::SurfaceType::Plane or
      surface.bounds().type() == Acts::SurfaceBounds::BoundsType::eAnnulus) {
    // This is not a template, hance properly placed
    if (not templateSurface) {
      Polyhedron surfaceHedron = surface.polyhedronRepresentation(gctx, surfaceStyle.nSegments);
      auto vertices3D = surfaceHedron.vertices;
      pSurface._vertices = vertices3D;
    } else {
      auto planarBounds =
          dynamic_cast<const Acts::PlanarBounds*>(&(surface.bounds()));
      if (planarBounds != nullptr) {
        auto vertices2D = planarBounds->vertices(surfaceStyle.nSegments);
        pSurface._vertices.reserve(vertices2D.size());
        for (const auto& v2 : vertices2D) {
          pSurface._vertices.push_back({v2[0], v2[1], 0.});
        }
      } else {
        auto annulusBounds =
            dynamic_cast<const Acts::AnnulusBounds*>(&(surface.bounds()));
        if (annulusBounds != nullptr) {           
          auto vertices2D = annulusBounds->vertices(surfaceStyle.nSegments);
          pSurface._vertices.reserve(vertices2D.size());
          std::cout << "** Annulus with " << vertices2D.size() << " vertices" << std::endl;
          for (const auto& v2 : vertices2D) {
            pSurface._vertices.push_back({v2[0], v2[1], 0.});
          }
        }
      }
    }
    // Translate the bound values
    const auto& boundValues = surface.bounds().values();
    if (surface.bounds().type() ==
        Acts::SurfaceBounds::BoundsType::eRectangle) {
      // Set the type
      pSurface._type = ProtoSurface::type::e_rectangle;
      // Set the measure
      pSurface._measures = {
          static_cast<actsvg::scalar>(0.5 * (boundValues[2] - boundValues[0])),
          static_cast<actsvg::scalar>(0.5 * (boundValues[3] - boundValues[1]))};
    } else if (surface.bounds().type() ==
               Acts::SurfaceBounds::BoundsType::eTrapezoid) {
      // Set the type
      pSurface._type = ProtoSurface::type::e_trapez;
      // Set the measure
      pSurface._measures = {static_cast<actsvg::scalar>(boundValues[0]),
                            static_cast<actsvg::scalar>(boundValues[1]),
                            static_cast<actsvg::scalar>(boundValues[2])};
    } else if (surface.bounds().type() ==
               Acts::SurfaceBounds::BoundsType::eAnnulus) {
      // Set the type
      pSurface._type = ProtoSurface::type::e_annulus;
      // Set the measure 
      //for (const auto& bv : boundValues) {
      //  pSurface._measures.push_back(static_cast<actsvg::scalar>(bv));
      //}
    }
  } else if (surface.type() == Acts::Surface::SurfaceType::Disc) {
    // Set the type
    pSurface._type = ProtoSurface::type::e_disc;

    const auto& boundValues = surface.bounds().values();
    if (surface.bounds().type() == Acts::SurfaceBounds::BoundsType::eDisc) {
      pSurface._radii = {static_cast<actsvg::scalar>(boundValues[0]),
                         static_cast<actsvg::scalar>(boundValues[1])};
      pSurface._opening = {
          static_cast<actsvg::scalar>(boundValues[3] - boundValues[2]),
          static_cast<actsvg::scalar>(boundValues[3] + boundValues[2])};
    }
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