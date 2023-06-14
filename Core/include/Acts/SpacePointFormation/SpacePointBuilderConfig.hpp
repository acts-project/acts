// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"

namespace Acts {
struct SpacePointBuilderConfig {
  /// Tracking geometry
  std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry;
  /// Perform the perpendicular projection for space point finding
  bool usePerpProj = false;

  SpacePointBuilderConfig() = default;
};

}  // namespace Acts
