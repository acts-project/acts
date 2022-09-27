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

Acts::Svg::ProtoSurface Acts::Svg::SurfaceConverter::convert(
    const GeometryContext& gctx, const Surface& surface,
    const SurfaceConverter::Options& surfaceOptions) {
  // The local logger
  ACTS_LOCAL_LOGGER(
      getDefaultLogger("SurfaceSvgConverter", surfaceOptions.logLevel));

  ProtoSurface pSurface;

  // In case of non-template surfaces, the polyhedron does the trick
  if (not surfaceOptions.templateSurface) {
    ACTS_DEBUG("Not building template surface!");

    // Check if the surface is regular for disc and cylinder
    auto sfTransform = surface.transform(gctx);
    bool regular = sfTransform.rotation().col(2).isApprox(Vector3(0., 0., 1.));
    auto sfBoundsType = surface.bounds().type();
    if (regular and (sfBoundsType == Acts::SurfaceBounds::eCylinder or
                     sfBoundsType == Acts::SurfaceBounds::eDisc)) {
      actsvg::scalar zpos =
          static_cast<actsvg::scalar>(sfTransform.translation().z());
      // For drawing regular cylinders and discs, access the values
      const auto& boundValues = surface.bounds().values();
      if (sfBoundsType == Acts::SurfaceBounds::eDisc) {
        // The radii
        actsvg::scalar ri = static_cast<actsvg::scalar>(boundValues[0]);
        actsvg::scalar ro = static_cast<actsvg::scalar>(boundValues[1]);
        pSurface._radii = {ri, ro};
        pSurface._zparameters = {zpos, 0.};
      } else {
        actsvg::scalar ro = static_cast<actsvg::scalar>(boundValues[0]);
        actsvg::scalar hl = static_cast<actsvg::scalar>(boundValues[1]);
        pSurface._radii = {0., ro};
        pSurface._zparameters = {zpos, hl};
      }
      pSurface._opening = {
          static_cast<actsvg::scalar>(boundValues[3] - boundValues[2]),
          static_cast<actsvg::scalar>(boundValues[3] + boundValues[2])};

    } else {
      // Polyhedron surface for vertices needed anyways
      Polyhedron surfaceHedron =
          surface.polyhedronRepresentation(gctx, surfaceOptions.style.nSegments);
      auto vertices3D = surfaceHedron.vertices;
      pSurface._vertices = vertices3D;
    }
  } else {
    // In case it's a template surface, only the bounds matter
    // Check if planar bounds
    auto planarBounds =
        dynamic_cast<const Acts::PlanarBounds*>(&(surface.bounds()));
    if (planarBounds != nullptr) {
      auto vertices2D = planarBounds->vertices(surfaceOptions.style.nSegments);
      pSurface._vertices.reserve(vertices2D.size());
      for (const auto& v2 : vertices2D) {
        pSurface._vertices.push_back({v2[0], v2[1], 0.});
      }
    } else {
      // Or annulus bounds
      auto annulusBounds =
          dynamic_cast<const Acts::AnnulusBounds*>(&(surface.bounds()));
      if (annulusBounds != nullptr) {
        auto vertices2D = annulusBounds->vertices(surfaceOptions.style.nSegments);
        pSurface._vertices.reserve(vertices2D.size());
        for (const auto& v2 : vertices2D) {
          pSurface._vertices.push_back({v2[0], v2[1], 0.});
        }
      } else if (surface.type() == Acts::Surface::SurfaceType::Disc) {
        // Or disc bounds
        const auto& boundValues = surface.bounds().values();
        if (surface.bounds().type() == Acts::SurfaceBounds::BoundsType::eDisc) {
          // The radii
          actsvg::scalar ri = static_cast<actsvg::scalar>(boundValues[0]);
          actsvg::scalar ro = static_cast<actsvg::scalar>(boundValues[1]);
          pSurface._radii = {ri, ro};
          pSurface._opening = {
              static_cast<actsvg::scalar>(boundValues[3] - boundValues[2]),
              static_cast<actsvg::scalar>(boundValues[3] + boundValues[2])};

          actsvg::scalar pl = pSurface._opening[0];
          actsvg::scalar ph = pSurface._opening[1];

          pSurface._vertices = {
              {static_cast<actsvg::scalar>(ri * std::cos(pl)),
               static_cast<actsvg::scalar>(ri * std::sin(pl)), 0.},
              {static_cast<actsvg::scalar>(ro * std::cos(ph)),
               static_cast<actsvg::scalar>(ro * std::sin(ph)), 0.},
              {static_cast<actsvg::scalar>(ri * std::cos(pl)),
               static_cast<actsvg::scalar>(ri * std::sin(pl)), 0.},
              {static_cast<actsvg::scalar>(ro * std::cos(ph)),
               static_cast<actsvg::scalar>(ro * std::sin(ph)), 0.}};
        }
      }
    }
  }

  // Bound types and values
  const auto& boundValues = surface.bounds().values();
  auto bType = surface.bounds().type();
  auto bValues = surface.bounds().values();
  if (bType == Acts::SurfaceBounds::BoundsType::eRectangle) {
    // Set the surface bounds type for ACTSVG
    pSurface._type = ProtoSurface::type::e_rectangle;
    // Set the measure parameters, rectangle needs special treatment
    pSurface._measures = {
        static_cast<actsvg::scalar>(0.5 * (boundValues[2] - boundValues[0])),
        static_cast<actsvg::scalar>(0.5 * (boundValues[3] - boundValues[1]))};
  } else {
    // Set the measure parameters
    for (const auto& bv : boundValues) {
      pSurface._measures.push_back(static_cast<actsvg::scalar>(bv));
    }
    // Set the surface bounds type for ACTSVG
    if (bType == Acts::SurfaceBounds::BoundsType::eTrapezoid) {
      pSurface._type = ProtoSurface::type::e_trapez;
    } else if (bType == Acts::SurfaceBounds::BoundsType::eAnnulus) {
      pSurface._type = ProtoSurface::type::e_trapez;
    } else if (bType == Acts::SurfaceBounds::BoundsType::eDisc) {
      pSurface._type = ProtoSurface::type::e_disc;
    } else if (bType == Acts::SurfaceBounds::BoundsType::eCylinder) {
      pSurface._type = ProtoSurface::type::e_cylinder;
    }
  }


  auto [ fill, stroke ] = surfaceOptions.style.fillAndStroke();

  pSurface._fill = fill;
  pSurface._stroke = stroke;

  return pSurface;
}
