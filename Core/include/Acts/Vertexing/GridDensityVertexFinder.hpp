// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Result.hpp"
#include "Acts/Vertexing/GaussianGridTrackDensity.hpp"
#include "Acts/Vertexing/Vertex.hpp"
#include "Acts/Vertexing/VertexingOptions.hpp"
#include "DummyVertexFitter.hpp"

namespace Acts {

/// @class TODO
template <int mainGridSize = 2000, int trkGridSize = 15,
          typename vfitter_t = DummyVertexFitter<>>
class GridDensityVertexFinder {
  using InputTrack_t = typename vfitter_t::InputTrack_t;
  using GridDensity = GaussianGridTrackDensity<mainGridSize, trkGridSize>;

 public:
  /// @brief The Config struct
  struct Config {
    Config(float zMinMax = 100)
        : zMinMax(zMinMax),
          gridDensity(typename GridDensity::Config(zMinMax)) {}
    // Min and max z value of big grid
    float zMinMax;  // mm

    GridDensity gridDensity;
  };

  /// @brief Function that finds single vertex candidate
  ///
  /// @param trackVector Input track collection
  /// @param vertexingOptions Vertexing options
  ///
  /// @return Vector of vertices, filled with a single
  ///         vertex (for consistent interfaces)
  Result<std::vector<Vertex<InputTrack_t>>> find(
      const std::vector<const InputTrack_t*>& trackVector,
      const VertexingOptions<InputTrack_t>& vertexingOptions) const;

  /// @brief Constructor used if InputTrack_t type == BoundParameters
  ///
  /// @param cfg Configuration object
  template <typename T = InputTrack_t,
            std::enable_if_t<std::is_same<T, BoundParameters>::value, int> = 0>
  GridDensityVertexFinder(const Config& cfg)
      : m_cfg(cfg), m_extractParameters([](T params) { return params; }) {}

  /// @brief Default constructor used if InputTrack_t type == BoundParameters
  template <typename T = InputTrack_t,
            std::enable_if_t<std::is_same<T, BoundParameters>::value, int> = 0>
  GridDensityVertexFinder()
      : m_extractParameters([](T params) { return params; }) {}

  /// @brief Constructor for user-defined InputTrack_t type =! BoundParameters
  ///
  /// @param cfg Configuration object
  /// @param func Function extracting BoundParameters from InputTrack_t object
  GridDensityVertexFinder(
      const Config& cfg,
      const std::function<BoundParameters(InputTrack_t)>& func)
      : m_cfg(cfg), m_extractParameters(func) {}

  /// @brief Constructor for user-defined InputTrack_t type =! BoundParameters
  /// with default Config object
  ///
  /// @param func Function extracting BoundParameters from InputTrack_t object
  GridDensityVertexFinder(
      const std::function<BoundParameters(InputTrack_t)>& func)
      : m_extractParameters(func) {}

 private:
  Config m_cfg;

  /// @brief Function to extract track parameters,
  /// InputTrack_t objects are BoundParameters by default, function to be
  /// overwritten to return BoundParameters for other InputTrack_t objects.
  ///
  /// @param InputTrack_t object to extract track parameters from
  std::function<BoundParameters(InputTrack_t)> m_extractParameters;
};

}  // namespace Acts

#include "GridDensityVertexFinder.ipp"
