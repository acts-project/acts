// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Vertexing/Vertex.hpp"

#include <map>

namespace Acts {

/// @brief Helper struct for storing vertex related information
template <typename input_track_t>
struct VertexInfo {
  VertexInfo() = default;

  VertexInfo(const Acts::Vertex<input_track_t>& vtx, const Acts::Vector4& pos)
      : constraintVertex(vtx),
        linPoint(pos),
        oldPosition(pos),
        seedPosition(pos) {}

  // The constraint vertex
  Acts::Vertex<input_track_t> constraintVertex;

  // The linearization point
  Acts::Vector4 linPoint{Acts::Vector4::Zero()};

  // Old position from last iteration
  Acts::Vector4 oldPosition{Acts::Vector4::Zero()};

  // The seed position
  Acts::Vector4 seedPosition{Acts::Vector4::Zero()};

  // Needs relinearization bool
  bool relinearize = true;

  // Vector of all track currently held by vertex
  std::vector<const input_track_t*> trackLinks;

  std::map<const input_track_t*, const BoundTrackParameters> ip3dParams;
};

}  // namespace Acts
