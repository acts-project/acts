// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Vertexing/TrackDensity.hpp"

namespace Acts {

/// @class GaussianTrackDensity
///
/// @brief Class to model tracks as 2D Gaussian-shaped density functions
/// based on their d0 and z0 perigee parameters (mean value)
/// and covariance matrices (determining the width of the function)
class GaussianTrackDensity {
 public:
  /// @brief The State struct
  struct State {
    State(unsigned int nTracks) : trackDensityState(nTracks) {}
    // The track density
    // Defaulted to Gaussian shaped density function
    TrackDensity trackDensity;

    // The track density state
    TrackDensity::State trackDensityState;
  };

  /// @brief The Config struct
  struct Config {
    // Maximum d0 impact parameter significance to use a track
    double d0MaxSignificance = 3.5;

    // Maximum z0 impact parameter significance to use a track
    double z0MaxSignificance = 12.;
  };

  /// Default constructor
  GaussianTrackDensity() = default;

  /// Constructor with config
  GaussianTrackDensity(const Config& cfg) : m_cfg(cfg) {}

  /// @brief Calculates the global maximum
  ///
  /// @param trackList The list of tracks
  /// @param state The GaussianTrackDensity state
  ///
  /// @return The z position of the maximum
  double globalMaximum(const std::vector<Acts::BoundParameters>& trackList,
                       State& state) const;

  /// @brief Calculates the global maximum with width
  ///
  /// @param trackList The list of tracks
  /// @param state The GaussianTrackDensity state
  ///
  /// @return The z position of the maximum and its width
  std::pair<double, double> globalMaximumWithWidth(
      const std::vector<Acts::BoundParameters>& trackList, State& state) const;

 private:
  /// The configuration
  Config m_cfg;

  /// @brief Adds tracks based on a significance cut
  ///
  /// @param trackList The list of tracks
  /// @param state The GaussianTrackDensity state
  void addTracks(const std::vector<Acts::BoundParameters>& trackList,
                 State& state) const;
};
}  // namespace Acts
