// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/StrawSurface.hpp"
#include "Acts/Surfaces/TrapezoidBounds.hpp"
#include "Acts/Utilities/UnitVectors.hpp"
#include "Acts/Visualization/EventDataVisualization.hpp"
#include "Acts/Visualization/IVisualization.hpp"

#include <fstream>
#include <sstream>
#include <string>

namespace Acts {
namespace EventDataVisualization {

/// Helper method to visualiza all types of surfaces
///
/// @param helper The visualziation helper
/// @param triangulate The directive whether to create triangular meshes
/// @param tag The test tag (mode) identification
/// @param suffix The file suffix for writing
/// @param msuffix the (optional) material file suffix
///
/// @return the total character count for regression testing, needs to be
/// updated if the output is about to change
static inline size_t test(IVisualization& helper) {
  std::stringstream cStream;
  cStream << std::setprecision(4);

  const IVisualization::ColorType pcolor = {20, 120, 20};
  const IVisualization::ColorType scolor = {235, 198, 52};

  auto gctx = GeometryContext();
  auto identity = std::make_shared<Transform3D>(Transform3D::Identity());

  // rectangle and plane
  auto rectangle = std::make_shared<RectangleBounds>(15., 15.);
  auto plane = Surface::makeShared<PlaneSurface>(identity, rectangle);

  // now create parameters on this surface
  // l_x, l_y, phi, theta, q/p (1/p), t
  std::array<double, 6> pars_array = {
      {-0.1234, 4.8765, 0.45, 0.888, 0.001, 21.}};
  BoundParameters::ParVector_t pars;
  pars << pars_array[0], pars_array[1], pars_array[2], pars_array[3],
      pars_array[4], pars_array[5];

  Acts::BoundSymMatrix cov;
  cov << 0.25, 0.0042, -0.00076, 6.156e-06, -2.11e-07, 0, 0.0042, 0.459,
      -0.000173, 0.000916, -4.017e-08, 0, -0.00076, -0.000173, 2.36 - 05,
      -2.76e-07, 1.12e-08, 0, 6.15e-06, 0.000916, -2.76e-07, 1.84e-06,
      -2.85e-11, 0, -2.11 - 07, -4.017e-08, 1.123e-08, -2.85 - 11, 1.26e-10, 0,
      0, 0, 0, 0, 0, 1;

  // constructor from parameter vector
  BoundParameters ataPlane(gctx, std::move(cov), pars, plane);

  Visualization::drawBoundParameters(helper, ataPlane, gctx, 0.02, 100., 10.,
                                     true, 72, pcolor, scolor);

  helper.write("EventData_BoundAtPlaneParameters");
  helper.write(cStream);

  helper.clear();

  return cStream.str().size();
}

}  // namespace EventDataVisualization
}  // namespace Acts