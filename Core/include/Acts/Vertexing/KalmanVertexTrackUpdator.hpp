// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "Acts/Vertexing/KalmanVertexUpdator.hpp"
#include "Acts/Vertexing/TrackAtVertex.hpp"
#include "Acts/Vertexing/Vertex.hpp"

namespace Acts {

/// @class KalmanVertexTrackUpdator
///
/// @brief Refits a single track with the knowledge of
/// the vertex it has originated from
/// Based on R. Frühwirth et al.
/// Vertex reconstruction and track bundling at the lep collider using
/// robust Algorithms Computer Physics Comm.: 96 (1996) 189, chapter 2.1
///
/// @tparam input_track_t Track object type

template <typename input_track_t>
class KalmanVertexTrackUpdator {
 public:
  /// @struct Configuration struct
  struct Config {
    /// Kalman vertex updator
    KalmanVertexUpdator<input_track_t> vtx_updator;
  };

  /// Constructor
  KalmanVertexTrackUpdator(const Config& config = Config()) : m_cfg(config) {}

  /// @brief Refits a single track with the knowledge of
  /// the vertex it has originated from
  ///
  /// @param gctx The Geometry Context
  /// @param track Track to update
  /// @param vtx Vertex `track` belongs to
  void update(const GeometryContext& gctx, TrackAtVertex<input_track_t>& track,
              const Vertex<input_track_t>* vtx) const;

 private:
  /// Configuration object
  const Config m_cfg;

  /// @brief Function to correct 2-pi periodicity for phi and theta
  ///
  /// @param phiIn Phi
  /// @param thetaIn Theta
  ///
  /// @return Pair of (corrected phi, corrected theta)
  std::pair<double, double> correctPhiThetaPeriodicity(double phiIn,
                                                       double thetaIn) const;
};

}  // Namespace Acts

#include "Acts/Vertexing/KalmanVertexTrackUpdator.ipp"
