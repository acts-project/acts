// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <map>
#include <set>
#include "Acts/EventData/TrackParameters.hpp"

namespace Acts {

/// @class TrackDensity
///
/// @brief Class to model tracks as 2D density functions based on
/// their d0 and z0 perigee parameters (mean value) and covariance
/// matrices (determining the width of the function)
class TrackDensity {
 public:
  /// @brief Struct to store information for a single track
  struct TrackEntry {
    TrackEntry() = default;
    TrackEntry(double c0in, double c1in, double c2in, double zMin, double zMax)
        : c0(c0in), c1(c1in), c2(c2in), lowerBound(zMin), upperBound(zMax) {}

    // Dummy constructor for binary search
    TrackEntry(double z) : lowerBound(z), upperBound(z) {}

    // Cached information for a single track
    // z-independent term in exponent
    double c0 = 0;
    // linear coefficient in exponent
    double c1 = 0;
    // quadratic coefficient in exponent
    double c2 = 0;
    double lowerBound = 0;
    double upperBound = 0;
  };

  /// @brief Functor to compare two z0 values
  struct predPerigee {
    bool operator()(const BoundParameters& left,
                    const BoundParameters& right) const {
      return left.parameters()[ParID_t::eLOC_Z0] <
             right.parameters()[ParID_t::eLOC_Z0];
    }
  };

  /// @brief Functor to compare two TrackEntry values based
  /// on their lower limits (low to high)
  struct predEntryByMin {
    bool operator()(const TrackEntry& left, const TrackEntry& right) const {
      return left.lowerBound < right.lowerBound;
    }
  };

  /// @brief Functor to compare two TrackEntry values based on
  /// their upper limits (low to high)
  struct predEntryByMax {
    bool operator()(const TrackEntry& left, const TrackEntry& right) const {
      return left.upperBound < right.upperBound;
    }
  };

  using TrackMap = std::map<BoundParameters, TrackEntry, predPerigee>;
  using LowerMap = std::map<TrackEntry, BoundParameters, predEntryByMax>;
  using UpperMap = std::map<TrackEntry, BoundParameters, predEntryByMin>;

  /// @brief The Config struct
  struct Config {
    // Assumed shape of density function:
    // Gaussian (true) or parabolic (false)
    bool isGaussianShaped = true;
  };

  /// @brief The State struct
  struct State {
    double maxZRange = 0;

    // Maps to cache track information
    std::set<BoundParameters, predPerigee> trackSet;
    LowerMap lowerMap;
    UpperMap upperMap;
  };

  /// Default constructor
  TrackDensity() = default;

  /// Constructor with config
  TrackDensity(const Config& cfg) : m_cfg(cfg) {}

  /// @brief Add a track to the set being considered
  ///
  /// @param state The track density state
  /// @param trk Track parameters.
  /// @param d0SignificanceCut Significance cut on d0.
  /// @param z0SignificanceCut Significance cut on z0.
  void addTrack(State& state, const BoundParameters& trk,
                double d0SignificanceCut, double z0SignificanceCut) const;

  /// @brief Calculates z position of global maximum with Gaussian width
  /// for density function.
  /// Strategy:
  /// The global maximum must be somewhere near a track.
  /// Since we can calculate the first and second derivatives, at each point we
  /// can determine a) whether the function is curved up (minimum) or down
  /// (maximum) b) the distance to nearest maximum, assuming either Newton
  /// (parabolic) or Gaussian local behavior.
  /// For each track where the second derivative is negative, find step to
  /// nearest maximum, take that step and then do one final refinement. The
  /// largest density encountered in this procedure (after checking all tracks)
  /// is considered the maximum.
  ///
  /// @param state The track density state
  ///
  /// @return Pair of position of global maximum and Gaussian width
  std::pair<double, double> globalMaximumWithWidth(State& state) const;

  /// @brief Calculates the z position of the global maximum
  ///
  /// @param state The track density state
  ///
  /// @return z position of the global maximum
  double globalMaximum(State& state) const;

  /// @brief Evaluate the density function at the specified
  /// coordinate along the beamline
  ///
  /// @param state The track density state
  /// @param z z-position along the beamline
  ///
  /// @return The track density value
  double trackDensity(State& state, double z) const;

  /// @brief Evaluate the density function and its two first
  /// derivatives at the specified coordinate along the beamline
  ///
  /// @param state The track density state
  /// @param z z-position along the beamline
  /// @param[out] firstDerivative The first derivative
  /// @param[out] secondDerivative The second derivative
  ///
  /// @return The track density value
  double trackDensity(State& state, double z, double& firstDerivative,
                      double& secondDerivative) const;

 private:
  /// The configuration
  Config m_cfg;

  /// @brief Update the current maximum values
  ///
  /// @param newZ The new z value
  /// @param newValue The new value at z position
  /// @param newSecondDerivative The new second derivative
  /// @param[out] maxZ The max z value
  /// @param[out] maxValue The max value at z position
  /// @param[out] maxSecondDerivative The max second derivative
  void updateMaximum(double newZ, double newValue, double newSecondDerivative,
                     double& maxZ, double& maxValue,
                     double& maxSecondDerivative) const;

  /// @brief Calculates the step size
  ///
  /// @param y Position value
  /// @param dy First derivative
  /// @param ddy Second derivative
  ///
  /// @return The step size
  double stepSize(double y, double dy, double ddy) const;
};

}  // namespace Acts