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

template <typename input_track_t>
class Chi2TrackCompatibilityEstimator {
 public:
  /// Default constructor
  Chi2TrackCompatibilityEstimator() = default;

  double getCompatibility(const GeometryContext& gctx,
                          const BoundParameters& track,
                          const Vector3D& vertexPos) const;

  void setTrackCompatibility(const GeometryContext& gctx,
                             TrackAtVertex<input_track_t>& trackAtVertex,
                             const Vector3D& vertexPos,
                             const std::function<BoundParameters(input_track_t)>
                                 extractParameters) const;
};

}  // namespace Acts

#include "Acts/Vertexing/Chi2TrackCompatibilityEstimator.ipp"