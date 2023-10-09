// This file is part of the Acts project.
//
// Copyright (C) 2019-2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Vertexing/Vertex.hpp"

// TODO use unordered map?
#include <map>

namespace Acts {

/// @brief Helper struct for storing vertex related information
template <typename input_track_t>
struct VertexInfo {
  VertexInfo() = default;

  VertexInfo(const Acts::Vertex<input_track_t>& vtx, const Acts::Vector4& pos)
      : constraintVertex(vtx),
        linPoint(pos),
        position(pos),
        seedPosition(pos) {}

  // Vertex constraint
  Acts::Vertex<input_track_t> constraintVertex;

  // Point where all associated tracks were linearized (?)
  Acts::Vector4 linPoint{Acts::Vector4::Zero()};

  // Current vertex position
  Acts::Vector4 position{Acts::Vector4::Zero()};

  // The seed position (i.e., the first estimate for the vertex position as
  // obtained by the vertex seed finder)
  Acts::Vector4 seedPosition{Acts::Vector4::Zero()};

  // If set to true, the associated tracks need to be relinearized at a
  // different point in space
  bool relinearize = true;

  // Vector of all tracks that are currently assigned to vertex
  std::vector<const input_track_t*> trackLinks;

  std::map<const input_track_t*, const BoundTrackParameters> impactParams3D;
};

}  // namespace Acts
