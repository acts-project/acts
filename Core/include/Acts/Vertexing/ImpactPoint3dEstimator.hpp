// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Vertexing/TrackAtVertex.hpp"

namespace Acts {

/// @class ImpactPoint3dEstimator
///
/// @brief Estimates point of closest approach in 3D
class ImpactPoint3dEstimator {
 public:
  ImpactPoint3dEstimator() = default;

  /// @brief Calculates 3D distance between a track and a 3D point
  ///
  /// @param params Track parameters
  /// @param refPos Position to calculate distance to
  ///
  /// @return Distance
  double calculateDistance(const BoundParameters& params,
                           const Vector3D& refPos) const;

  /// @brief Creates track parameters bound to plane
  /// at point of closest approach in 3d to given
  /// reference position. The parameters and errors
  /// are defined on the plane intersecting the track
  /// at point of closest approach, with track ortogonal
  /// to the plane and center of the plane defined as the
  /// given reference point (vertex).
  ///
  /// @param trk Track at vertex
  /// @param refPos Reference position (vertex)
  ///
  /// @return New track params
  // BoundParameters paramsAtIP3d(const TrackAtVertex& trk,
  //							const Vector3D& refPos)
  // const;
};

}  // namespace Acts
