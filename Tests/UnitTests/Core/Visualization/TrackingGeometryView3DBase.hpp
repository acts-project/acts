// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Tests/CommonHelpers/CylindricalTrackingGeometry.hpp"
#include "Acts/Visualization/GeometryView3D.hpp"
#include "Acts/Visualization/IVisualization3D.hpp"

#include <fstream>
#include <sstream>
#include <string>

namespace Acts::TrackingGeometryView3DTest {

GeometryContext tgContext = GeometryContext();

Test::CylindricalTrackingGeometry cGeometry(tgContext);
auto tGeometry = cGeometry();

/// Helper method to visualize all types of surfaces
///
/// @param helper The visualization helper
/// @param triangulate The directive whether to create triangular meshes
/// @param tag The test tag (mode) identification
///
/// @return a string containing all written characters

static inline std::string run(IVisualization3D& helper, bool triangulate,
                              const std::string& tag) {
  std::stringstream cStream;

  ViewConfig viewSensitive = ViewConfig({0, 180, 240});
  viewSensitive.triangulate = triangulate;
  ViewConfig viewPassive = ViewConfig({240, 280, 0});
  viewPassive.triangulate = triangulate;
  ViewConfig viewVolume = ViewConfig({220, 220, 0});
  viewVolume.triangulate = triangulate;
  ViewConfig viewContainer = ViewConfig({220, 220, 0});
  viewContainer.triangulate = triangulate;
  ViewConfig viewGrid = ViewConfig({220, 0, 0});
  viewGrid.nSegments = 8;
  viewGrid.offset = 3.;
  viewGrid.triangulate = triangulate;

  const Acts::TrackingVolume& tgVolume = *(tGeometry->highestTrackingVolume());

  GeometryView3D::drawTrackingVolume(helper, tgVolume, tgContext, viewContainer,
                                     viewVolume, viewPassive, viewSensitive,
                                     viewGrid, true, tag);

  GeometryView3D::drawTrackingVolume(helper, tgVolume, tgContext, viewContainer,
                                     viewVolume, viewPassive, viewSensitive,
                                     viewGrid, false);
  helper.write(cStream);

  return cStream.str();
}

}  // namespace Acts::TrackingGeometryView3DTest
