// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Visualization/ViewConfig.hpp"

#include "Acts/Geometry/GeometryObject.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Surfaces/Surface.hpp"

#include <functional>

namespace Acts {

ViewConfig defaultGeometryColoring(const GeometryObject& geoObj) {
  if (geoObj.geometryId().boundary() != 0) {
    return ViewConfig{.color = {255, 165, 0}};
  }

  if (geoObj.geometryId().sensitive() != 0) {
    return ViewConfig{.color = {0, 255, 0}};
  }

  else {
    return ViewConfig{.visible = false};
  }
}

}  // namespace Acts
