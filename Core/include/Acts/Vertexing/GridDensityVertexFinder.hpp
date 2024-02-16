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
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Utilities/Result.hpp"
#include "Acts/Vertexing/DummyVertexFitter.hpp"
#include "Acts/Vertexing/GaussianGridTrackDensity.hpp"
#include "Acts/Vertexing/IVertexFinder.hpp"
#include "Acts/Vertexing/Vertex.hpp"
#include "Acts/Vertexing/VertexingOptions.hpp"

#include <map>

namespace Acts {

// namespace detail {
// template <typename derived_t>
// }

/// @class GridDensityVertexFinder
/// @brief Vertex finder that makes use of a track density grid.
/// Each single track is modelled as a 2(!)-dim Gaussian distribution grid
/// in the d0-z0 plane, but only the overlap with the z-axis (i.e. a 1-dim
/// density vector) needs to be calculated. All track contributions along the
/// beam axis (main density grid) a superimposed and the z-value of the bin
/// with the highest track density is returned as a vertex candidate.
/// @tparam mainGridSize The size of the z-axis 1-dim main density grid
/// @tparam trkGridSize The 2(!)-dim grid size of a single track, i.e.
/// a single track is modelled as a (trkGridSize x trkGridSize) grid
/// in the d0-z0 plane. Note: trkGridSize has to be an odd value.
template <int mainGridSize = 2000, int trkGridSize = 15>
class GridDensityVertexFinder final : public IVertexFinder {
  // Assert odd trkGridSize
  static_assert(trkGridSize % 2);
  // Assert bigger main grid than track grid
  static_assert(mainGridSize > trkGridSize);

  using GridDensity = GaussianGridTrackDensity<mainGridSize, trkGridSize>;

 public:
  using MainGridVector = typename GridDensity::MainGridVector;
  using TrackGridVector = typename GridDensity::TrackGridVector;

  /// @brief The Config struct
  struct Config {
    ///@param zMinMax min and max z value of big z-axis grid
    Config(float zMinMax = 100)
        : gridDensity(typename GridDensity::Config(zMinMax)) {}
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

    // Function to extract parameters from InputTrack
    InputTrack::Extractor extractParameters;
  };

  /// @brief The State struct
  ///
  /// Only needed if cacheGridStateForTrackRemoval == true
  struct State {
    // The main density grid
    MainGridVector mainGrid = MainGridVector::Zero();
    // Map to store z-bin and track grid (i.e. the density contribution of
    // a single track to the main grid) for every single track
    std::map<InputTrack, std::pair<int, TrackGridVector>> binAndTrackGridMap;

    // Map to store bool if track has passed track selection or not
    std::map<InputTrack, bool> trackSelectionMap;

    // Store tracks that have been removed from track collection. These
    // track will be removed from the main grid
    std::vector<InputTrack> tracksToRemove;

    bool isInitialized = false;
  };

  /// @brief Function that finds single vertex candidate
  ///
  /// @param trackVector Input track collection
  /// @param vertexingOptions Vertexing options
  /// @param anyState The state object to cache the density grid
  /// and density contributions of each track, to be used
  /// if cacheGridStateForTrackRemoval == true
  ///
  /// @return Vector of vertices, filled with a single
  ///         vertex (for consistent interfaces)
  Result<std::vector<Vertex>> find(
      const std::vector<InputTrack>& trackVector,
      const VertexingOptions& vertexingOptions,
      IVertexFinder::State& anyState) const override;

  IVertexFinder::State makeState(
      const Acts::MagneticFieldContext& /*mctx*/) const override {
    return IVertexFinder::State{State{}};
  }

  void setTracksToRemove(
      IVertexFinder::State& anyState,
      const std::vector<InputTrack>& removedTracks) const override {
    auto& state = anyState.template as<State>();
    state.tracksToRemove = removedTracks;
  }

  /// @brief Constructor for user-defined InputTrack type
  ///
  /// @param cfg Configuration object
  /// @param func Function extracting BoundTrackParameters from InputTrack
  ///             object
  GridDensityVertexFinder(const Config& cfg) : m_cfg(cfg) {
    if (!m_cfg.extractParameters.connected()) {
      throw std::invalid_argument(
          "GridDensityVertexFinder: "
          "No track parameter extractor provided.");
    }
  }

 private:
  /// @brief Checks if a track passes the selection criteria for seeding
  ///
  /// @param trk The track
  ///
  /// @return Bool track passes selection
  bool doesPassTrackSelection(const BoundTrackParameters& trk) const;

  // The configuration object
  const Config m_cfg;
};

}  // namespace Acts

#include "GridDensityVertexFinder.ipp"
