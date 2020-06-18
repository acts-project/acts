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
#include "Acts/Visualization/GeometryView.hpp"
#include "Acts/Visualization/IVisualization.hpp"

#include <fstream>
#include <sstream>
#include <string>

namespace Acts {

namespace TrackingGeometryViewTest {

GeometryContext tgContext = GeometryContext();

Test::CylindricalTrackingGeometry cGeometry(tgContext);
auto tGeometry = cGeometry();

/// Helper method to visualiza all types of surfaces
///
/// @param helper The visualziation helper
/// @param triangulate The directive whether to create triangular meshes
/// @param tag The test tag (mode) identification
///
/// @return a string containing all written caracters

static inline std::string run(IVisualization& helper, bool triangulate,
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

  GeometryView::drawTrackingVolume(helper, tgVolume, tgContext, viewContainer,
                                   viewVolume, viewPassive, viewSensitive,
                                   viewGrid, true, tag);

  GeometryView::drawTrackingVolume(helper, tgVolume, tgContext, viewContainer,
                                   viewVolume, viewPassive, viewSensitive,
                                   viewGrid, false);
  helper.write(cStream);

  return cStream.str();
}

}  // namespace TrackingGeometryViewTest
}  // namespace Acts