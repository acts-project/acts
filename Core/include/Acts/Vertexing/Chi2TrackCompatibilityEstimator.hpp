// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <functional>
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Vertexing/TrackAtVertex.hpp"

namespace Acts {

/// @class Chi2TrackCompatibilityEstimator
///
/// @brief Estimates the compatibility of a track
///        to a vertex position based on the 3d
///        distance between the track and the vertex
template <typename input_track_t>
class Chi2TrackCompatibilityEstimator {
 public:
  /// Default constructor
  Chi2TrackCompatibilityEstimator() = default;

  /// @brief Estimates the compatibility value
  ///
  /// @param gctx The Geometry context
  /// @param track Track parameters at point of closest approach in
  ///   3d as retrieved by ImpactPoint3dEstimator::getParamsAtIP3d
  /// @param vertexPos The vertex position
  ///
  /// @return The compatibility value
  double getVtxCompatibility(const GeometryContext& gctx,
                             const BoundParameters* track,
                             const Vector3D& vertexPos) const;
};

}  // namespace Acts

#include "Acts/Vertexing/Chi2TrackCompatibilityEstimator.ipp"