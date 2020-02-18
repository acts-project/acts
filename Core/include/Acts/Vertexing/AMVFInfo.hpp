// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Definitions.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Vertexing/Vertex.hpp"

namespace Acts {

/// @brief Helper struct for storing vertex related information
  template<typename input_track_t>
  struct VertexInfo {
    // The linearization point
    Acts::SpacePointVector linPoint{Acts::SpacePointVector::Zero()};

    // The constraint vertex
    Acts::Vertex<input_track_t> constraintVertex;

    // Old position from last iteration
    Acts::SpacePointVector oldPosition;

    Acts::SpacePointVector seedPosition;

    // Needs relinearization bool
    bool relinearize;
  };

  /// @brief Helper struct for storing TrackAtVertex related
  template<typename input_track_t>
  struct TrackAtVertexInfo {
    // Links to vertices currently using the TrackAtVertex object
    std::vector<Vertex<input_track_t>*> linksToVertices;

    // Track parameters at point of closest approach in 3d as
    // retrieved by ImpactPoint3dEstimator::getParamsAtClosestApproach
    std::unique_ptr<const BoundParameters> ip3dParams;
  };

} // Acts namespace