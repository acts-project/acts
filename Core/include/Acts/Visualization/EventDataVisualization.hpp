// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Visualization/GeometryVisualization.hpp"
#include "Acts/Geometry/Polyhedron.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Visualization/IVisualization.hpp"

namespace Acts {

namespace Visualization {

/// Helper method to Surface objects
///
/// @param helper [in, out] The visualization helper
/// @param parameters The bound parameters to be drawn
/// @param gctx The geometry context for which it is drawn
/// @param momentumScale The scale of the momentum 
/// @param locErrorScale  The scale of the local error
/// @param angularErrorScale The sclae of the angular error
/// @param drawParameterSurface The indicator whether to draw the surface as well
/// @param lseg The number of segments for a full arch (if needed)
/// @param pcolor the (optional) color of the parameters to be written
/// @param scolor the (optional) color of the surface to be written
static inline void drawBoundParameters(
    IVisualization& helper, const SingleBoundTrackParameters& parameters, 
    const GeometryContext& gctx = GeometryContext(),
    double momentumScale = 1., double locErrorScale = 1., double angularErrorScale = 1.,
    bool drawParameterSurface = true,  size_t lseg = 72,
    const IVisualization::ColorType& pcolor = { 20, 120, 20}
    const IVisualization::ColorType& scolor = { 120, 20, 20}) {

  // First, if necessary, draw the surface
  if (drawParameterSurface){
      drawSurface()
  }

  // Drawing the polyhedron representation of surfaces
  Polyhedron phedron = surface.polyhedronRepresentation(gctx, lseg);
  phedron.draw(helper, false, color);
}