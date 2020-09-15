// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Vertexing/Vertex.hpp"

#include <map>

namespace Acts {

/// @brief Helper struct for storing vertex related information
template <typename input_track_t>
struct VertexInfo {
  VertexInfo() = default;

  VertexInfo(const Acts::Vertex<input_track_t>& vtx, const Acts::Vector4D& pos)
      : constraintVertex(vtx),
        linPoint(pos),
        oldPosition(pos),
        seedPosition(pos) {}

  // The constraint vertex
  Acts::Vertex<input_track_t> constraintVertex;

  // The linearization point
  Acts::Vector4D linPoint{Acts::Vector4D::Zero()};

  // Old position from last iteration
  Acts::Vector4D oldPosition{Acts::Vector4D::Zero()};

  // The seed position
  Acts::Vector4D seedPosition{Acts::Vector4D::Zero()};

  // Needs relinearization bool
  bool relinearize = true;

  // Vector of all track currently held by vertex
  std::vector<const input_track_t*> trackLinks;

  std::map<const input_track_t*, const BoundTrackParameters> ip3dParams;
};

}  // namespace Acts