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
#include "Acts/Vertexing/AdaptiveGridTrackDensity.hpp"
#include "Acts/Vertexing/DummyVertexFitter.hpp"
#include "Acts/Vertexing/Vertex.hpp"
#include "Acts/Vertexing/VertexingOptions.hpp"

#include <unordered_map>

namespace Acts {

/// @class AdaptiveGridDensityVertexFinder
/// @brief Vertex finder that makes use of a track density grid.
/// Each single track is modelled as a 2-dim Gaussian distribution grid
/// in the d0-z0 plane, but only the overlap with the z-axis (i.e. a 1-dim
/// density vector) needs to be calculated. All track contributions along the
/// beam axis (main density vector) are superimposed and the z-value of the bin
/// with the highest track density is returned as a vertex candidate.
/// Unlike the GridDensityVertexFinder, this seeder implements an adaptive
/// version where the density grid grows bigger with added tracks.
///
/// @tparam trkGridSize The 2-dim grid size of a single track, i.e.
/// a single track is modelled as a (trkGridSize x trkGridSize) grid
/// in the d0-z0 plane. Note: trkGridSize has to be an odd value.
template <int trkGridSize = 15, typename vfitter_t = DummyVertexFitter<>>
class AdaptiveGridDensityVertexFinder {
  // Assert odd trkGridSize
  static_assert(trkGridSize % 2);

  using InputTrack_t = typename vfitter_t::InputTrack_t;
  using GridDensity = AdaptiveGridTrackDensity<trkGridSize>;

 public:
  using DensityMap = typename GridDensity::DensityMap;

  /// @brief The Config struct
  struct Config {
    ///@param binSize Bin size of grid in mm
    Config(float binSize = 0.1)
        : gridDensity(typename GridDensity::Config(binSize)) {}
    ///@param gDensity The grid density
    Config(const GridDensity& gDensity) : gridDensity(gDensity) {}

    // The grid density object
    GridDensity gridDensity;

    // Cache the main grid and the density contributions (trackGrid and z-bin)
    // for every single track.
    // This option enables the possibility to calculate the entire main grid
    // only once in the first iteration. If tracks are removed from the track
    // collection, the individual track density contributions to the main grid
    // can just be removed without calculating the entire grid from scratch.
    bool cacheGridStateForTrackRemoval = true;

    // Maximum d0 impact parameter significance to use a track
    double maxD0TrackSignificance = 3.5;
    // Maximum z0 impact parameter significance to use a track
    double maxZ0TrackSignificance = 12.;
    // The actual corresponding cut values in the algorithm
    double d0SignificanceCut = maxD0TrackSignificance * maxD0TrackSignificance;
    double z0SignificanceCut = maxZ0TrackSignificance * maxZ0TrackSignificance;
    bool estimateSeedWidth = false;
  };

  /// @brief The State struct
  ///
  /// Only needed if cacheGridStateForTrackRemoval == true
  struct State {
    // Map from the z bin values to the corresponding track density
    DensityMap mainDensityMap;

    // Map from input track to corresponding track density map
    std::unordered_map<const InputTrack_t*, DensityMap> trackDensities;

    // Map to store bool if track has passed track selection or not
    std::unordered_map<const InputTrack_t*, bool> trackSelectionMap;

    // Store tracks that have been removed from track collection. These
    // tracks will be removed from the main grid
    std::vector<const InputTrack_t*> tracksToRemove;

    bool isInitialized = false;
  };

  /// @brief Function that finds single vertex candidate
  ///
  /// @param trackVector Input track collection
  /// @param vertexingOptions Vertexing options
  /// @param state The state object to cache the density grid
  /// and density contributions of each track, to be used
  /// if cacheGridStateForTrackRemoval == true
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
  AdaptiveGridDensityVertexFinder(const Config& cfg)
      : m_cfg(cfg), m_extractParameters([](T params) { return params; }) {}

  /// @brief Default constructor used if InputTrack_t type ==
  /// BoundTrackParameters
  template <
      typename T = InputTrack_t,
      std::enable_if_t<std::is_same<T, BoundTrackParameters>::value, int> = 0>
  AdaptiveGridDensityVertexFinder()
      : m_extractParameters([](T params) { return params; }) {}

  /// @brief Constructor for user-defined InputTrack_t type =!
  /// BoundTrackParameters
  ///
  /// @param cfg Configuration object
  /// @param func Function extracting BoundTrackParameters from InputTrack_t
  /// object
  AdaptiveGridDensityVertexFinder(
      const Config& cfg,
      const std::function<BoundTrackParameters(InputTrack_t)>& func)
      : m_cfg(cfg), m_extractParameters(func) {}

  /// @brief Constructor for user-defined InputTrack_t type =!
  /// BoundTrackParameters with default Config object
  ///
  /// @param func Function extracting BoundTrackParameters from InputTrack_t
  /// object
  AdaptiveGridDensityVertexFinder(
      const std::function<BoundTrackParameters(InputTrack_t)>& func)
      : m_extractParameters(func) {}

 private:
  /// @brief Checks if a track passes the selection criteria for seeding
  ///
  /// @param trk The track
  ///
  /// @return Bool track passes selection
  bool doesPassTrackSelection(const BoundTrackParameters& trk) const;

  // The configuration object
  const Config m_cfg;

  /// @brief Function to extract track parameters,
  /// InputTrack_t objects are BoundTrackParameters by default, function to be
  /// overwritten to return BoundTrackParameters for other InputTrack_t objects.
  std::function<BoundTrackParameters(InputTrack_t)> m_extractParameters;
};

}  // namespace Acts

#include "AdaptiveGridDensityVertexFinder.ipp"
