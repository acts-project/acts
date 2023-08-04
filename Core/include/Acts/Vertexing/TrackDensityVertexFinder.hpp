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
template <typename vfitter_t,
          typename track_density_t =
              GaussianTrackDensity<typename vfitter_t::InputTrack_t>>
class TrackDensityVertexFinder {
  // Provided vertex fitter type should comply with the VertexFitterConcept
  // to ensure providing an input track type InputTrack_t

  // static_assert(VertexFitterConcept<vfitter_t>,
  //              "Vertex fitter does not fulfill vertex fitter concept.");

  using InputTrack_t = typename vfitter_t::InputTrack_t;

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
  Result<std::vector<Vertex<InputTrack_t>>> find(
      const std::vector<const InputTrack_t*>& trackVector,
      const VertexingOptions<InputTrack_t>& vertexingOptions,
      State& state) const;

  /// @brief Constructor used if InputTrack_t type == BoundTrackParameters
  ///
  /// @param cfg Configuration object
  template <
      typename T = InputTrack_t,
      std::enable_if_t<std::is_same<T, BoundTrackParameters>::value, int> = 0>
  TrackDensityVertexFinder(const Config& cfg)
      : m_cfg(cfg), m_extractParameters([](T params) { return params; }) {}

  /// @brief Default constructor used if InputTrack_t type ==
  /// BoundTrackParameters
  template <
      typename T = InputTrack_t,
      std::enable_if_t<std::is_same<T, BoundTrackParameters>::value, int> = 0>
  TrackDensityVertexFinder()
      : m_extractParameters([](T params) { return params; }) {}

  /// @brief Constructor for user-defined InputTrack_t type =!
  /// BoundTrackParameters
  ///
  /// @param cfg Configuration object
  /// @param func Function extracting BoundTrackParameters from InputTrack_t
  /// object
  TrackDensityVertexFinder(
      const Config& cfg,
      const std::function<BoundTrackParameters(InputTrack_t)>& func)
      : m_cfg(cfg), m_extractParameters(func) {}

  /// @brief Constructor for user-defined InputTrack_t type =!
  /// BoundTrackParameters with default Config object
  ///
  /// @param func Function extracting BoundTrackParameters from InputTrack_t
  /// object
  TrackDensityVertexFinder(
      const std::function<BoundTrackParameters(InputTrack_t)>& func)
      : m_extractParameters(func) {}

 private:
  Config m_cfg;

  /// @brief Function to extract track parameters,
  /// InputTrack_t objects are BoundTrackParameters by default, function to be
  /// overwritten to return BoundTrackParameters for other InputTrack_t objects.
  std::function<BoundTrackParameters(InputTrack_t)> m_extractParameters;
};

}  // namespace Acts

#include "TrackDensityVertexFinder.ipp"
