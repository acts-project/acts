// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Vertexing/KalmanVertexTrackUpdator.hpp"
#include "Acts/Vertexing/TrackAtVertex.hpp"
#include "Acts/Vertexing/Vertex.hpp"

namespace Acts {

/// @class SequentialVertexSmoother
///
/// @brief Takes tracks from vertex candidate
/// and updates them with knowledge of the reconstructed
/// vertex position
///
/// @tparam input_track_t Track object type

template <typename input_track_t>
class SequentialVertexSmoother {
 public:
  struct Config {
    KalmanVertexTrackUpdator<input_track_t> trackUpdator;
  };

  /// @brief Default constructor
  SequentialVertexSmoother(const Config& cfg = Config()) : m_cfg(cfg) {}

  /// @brief Updates all tracks at vertex
  /// with knowledge of the vertex position
  ///
  /// @param gctx The Geometry Context
  /// @param vtx The vertex
  void smooth(const GeometryContext& gctx, Vertex<input_track_t>* vtx) const;

 private:
  /// Configuration object
  const Config m_cfg;
};

}  // Namespace Acts

#include "SequentialVertexSmoother.ipp"
