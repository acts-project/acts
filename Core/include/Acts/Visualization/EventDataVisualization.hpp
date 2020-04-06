// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/Polyhedron.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/UnitVectors.hpp"
#include "Acts/Visualization/GeometryVisualization.hpp"
#include "Acts/Visualization/IVisualization.hpp"

namespace Acts {
namespace Visualization {

/// Helper method to draw error ellipse
///
/// @param helper [in, out] The visualization helper
/// @param covariance The covariance matrix
/// @param referenceFrame The reference Frame
/// @param nsigma The sigmas to be drawn
/// @param locErrorScale The local Error scale
/// @param outofplane The out of plane drawning option
static inline void drawLocalCovariance(
    IVisualization& helper, const ActsSymMatrixD<2>& covariance,
    const RotationMatrix3D& referenceFrame, std::vector<int> nsigma = {3},
    double locErrorScale = 1, size_t lseg = 72,
    const IVisualization::ColorType& color = {20, 120, 20},
    double outofplane = 0.1) {
  double c00 = covariance(eLOC_0, eLOC_0);
  double c11 = covariance(eLOC_1, eLOC_1);
  double c01 = covariance(eLOC_0, eLOC_1);
  double r01 = c01 / (std::sqrt(c00) * std::sqrt(c11));
}

/// Helper method to Surface objects
///
/// @param helper [in, out] The visualization helper
/// @param parameters The bound parameters to be drawn
/// @param gctx The geometry context for which it is drawn
/// @param momentumScale The scale of the momentum
/// @param locErrorScale  The scale of the local error
/// @param angularErrorScale The sclae of the angular error
/// @param drawParameterSurface The indicator whether to draw the surface as
/// well
/// @param lseg The number of segments for a full arch (if needed)
/// @param pcolor the (optional) color of the parameters to be written
/// @param scolor the (optional) color of the surface to be written
template <typename parameters_t>
static inline void drawBoundParameters(
    IVisualization& helper, const parameters_t& parameters,
    const GeometryContext& gctx = GeometryContext(), double momentumScale = 1.,
    double locErrorScale = 1., double angularErrorScale = 1.,
    bool drawParameterSurface = true, size_t lseg = 72,
    const IVisualization::ColorType& pcolor = {20, 120, 20},
    const IVisualization::ColorType& scolor = {235, 198, 52}) {
  // First, if necessary, draw the surface
  if (drawParameterSurface) {
    drawSurface(helper, parameters.referenceSurface(), gctx,
                Transform3D::Identity(), lseg, false, scolor);
  }

  // Draw the parameter shaft and cone
  auto position = parameters.position();
  auto direction = parameters.momentum().normalized();
  double p = parameters.momentum().norm();
  drawArrowForward(helper, position, p * momentumScale * direction, 0.1, 0.1,
                   4., 72, pcolor);
}

}  // namespace Visualization
}  // namespace Acts