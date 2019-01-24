// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Vertexing/TrackAtVertex.hpp"
#include "Acts/Vertexing/Vertex.hpp"

namespace Acts {

/// @brief Helper struct for storing MultiAdaptiveVertexing
/// related vertex infos
template <typename input_track_t>
struct MAVFVertexInfo {
  // The seed vertex position
  Acts::Vector3D seedPos;
  // The linearization point
  Acts::Vector3D linPoint;
  // The constraint vertex
  Acts::Vertex<input_track_t> constraintVertex;
  // Is initialized bool
  bool isInitialized;
  // Old position from last iteration
  Acts::Vector3D oldPosition;
  // Needs relinearization bool
  bool relinearize;
};

/// @brief Helper struct for storing MultiAdaptiveVertexing
/// related TrackAtVertex infos
template <typename input_track_t>
struct MAVFTrackAtVtxInfo {
  // Map to store a vector of links to all vertices that use the key
  // TrackAtVertex object
  std::map<TrackAtVertex<input_track_t>*, std::vector<Vertex<input_track_t>*>>
      linkMapToVertices;
};

}  // namespace Acts
