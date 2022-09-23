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
  /// vertex position
  Vector3 vertex = {0., 0., 0.};
  /// Accepted squared difference in theta for two clusters
  double diffTheta2 = 1.;
  /// Accepted squared difference in phi for two clusters
  double diffPhi2 = 1.;
  /// Accepted distance between two clusters
  double diffDist = 100. * UnitConstants::mm;
  /// Allowed increase of strip length
  double stripLengthTolerance = 0.01;
  /// Allowed increase of strip length wrt gaps between strips
  double stripLengthGapTolerance = 0.01;
  /// Perform the perpendicular projection for space point finding
  bool usePerpProj = false;

  SpacePointBuilderConfig() = default;
};

}  // namespace Acts
