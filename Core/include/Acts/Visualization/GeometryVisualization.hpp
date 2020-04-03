// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/Polyhedron.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Visualization/IVisualization.hpp"

namespace Acts {

namespace Visualization {

/// Helper method to Surface objects
///
/// @param helper [in, out] The visualization helper
/// @param surface The surface to be drawn
/// @param gctx The geometry context for which it is drawn
/// @param transform An option additional transform
/// @param lseg The number of segments for a full arch (if needed)
/// @param triangulate The (optional) boolean diretive to triangulate the faces
/// @param color the (optional) color of the object to be written
static inline void drawSurface(
    IVisualization& helper, const Surface& surface, const GeometryContext& gctx,
    const Transform3D& transform = Transform3D::Identity(), size_t lseg = 72,
    bool triangulate = false,
    const IVisualization::ColorType& color = {120, 120, 120}) {
  // Drawing the polyhedron representation of surfaces
  Polyhedron phedron = surface.polyhedronRepresentation(gctx, lseg);
  phedron.move(transform);
  phedron.draw(helper, triangulate, color);
}

/// Helper method for volume objects
///
/// @tparam volume_t the templated volume class
///
/// @param helper [in, out] The visualization helper
/// @param volume The surface to be drawn
/// @param gctx The geometry context for which it is drawn
/// @param transform An option additional transform
/// @param lseg The number of segments for a full arch (if needed)
/// @param triangulate The (optional) boolean diretive to triangulate the faces
/// @param color the (optional) color of the object to be written
template <typename volume_t>
inline void drawVolume(IVisualization& helper, const volume_t& volume,
                       const GeometryContext& gctx,
                       const Transform3D& transform = Transform3D::Identity(),
                       size_t lseg = 72, bool triangulate = false,
                       const IVisualization::ColorType& color = {120, 120,
                                                                 120}) {
  // Drawing the polyhedron representation of surfaces
  auto bSurfaces = volume.boundarySurfaces();
  for (const auto& bs : bSurfaces) {
    drawSurface(helper, bs->surfaceRepresentation(), gctx, transform, lseg,
                triangulate, color);
  }
}

}  // namespace Visualization

}  // namespace Acts