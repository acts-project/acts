// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Utilities/Result.hpp"
#include "Acts/Vertexing/GaussianTrackDensity.hpp"
#include "Acts/Vertexing/Vertex.hpp"
#include "Acts/Vertexing/VertexFitterConcept.hpp"
#include "Acts/Vertexing/VertexingOptions.hpp"

namespace Acts {

/// @class TrackDensityVertexFinder
///
/// @brief Finds a vertex seed based on the maximum
/// of a track density function.
/// Each track is modelled as a 2d density
/// function around its d0/z0 perigee parameter values.
/// The z seed position is then found as the position
/// of the maximum of all summed track density functions.
///
/// Ref. (1): https://cds.cern.ch/record/2670380
///
/// @tparam vfitter_t The vertex fitter type (needed to fulfill concept)
/// @tparam track_density_t The track density type
template <typename vfitter_t, typename track_density_t = GaussianTrackDensity>
class TrackDensityVertexFinder {
  // Provided vertex fitter type should comply with the VertexFitterConcept
  // to ensure providing an input track type InputTrack_t

  // static_assert(VertexFitterConcept<vfitter_t>,
  //              "Vertex fitter does not fulfill vertex fitter concept.");

 public:
  /// @brief The Config struct
  struct Config {
    // The track density estimator
    track_density_t trackDensityEstimator;
  };

  /// State struct for fulfilling interface
  struct State {};

  /// @brief Function that finds single vertex candidate
  ///
  /// @param trackVector Input track collection
  /// @param vertexingOptions Vertexing options
  /// @param state State for fulfilling interfaces
  ///
  /// @return Vector of vertices, filled with a single
  ///         vertex (for consistent interfaces)
  Result<std::vector<Vertex>> find(const std::vector<InputTrack>& trackVector,
                                   const VertexingOptions& vertexingOptions,
                                   State& state) const;

  /// @brief Constructor for user-defined InputTrack type
  ///
  /// @param cfg Configuration object
  /// @param func Function extracting BoundTrackParameters from InputTrack
  ///             object
  TrackDensityVertexFinder(
      const Config& cfg,
      const std::function<BoundTrackParameters(const InputTrack&)>& func)
      : m_cfg(cfg), m_extractParameters(func) {}

 private:
  Config m_cfg;

  /// @brief Function to extract track parameters,
  std::function<BoundTrackParameters(const InputTrack&)> m_extractParameters;
};

}  // namespace Acts

#include "TrackDensityVertexFinder.ipp"
