// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/TransformationHelpers.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Utilities/Helpers.hpp"

#include <unordered_map>

namespace Acts::MultiTrajectoryHelpers {

/// @brief Struct for brief trajectory summary info
///
struct TrajectoryState {
  /// Total number of track states in trajectory
  std::size_t nStates = 0;
  /// Number of track states with measurements
  std::size_t nMeasurements = 0;
  /// Number of track states marked as outliers
  std::size_t nOutliers = 0;
  /// Number of hole states in trajectory
  std::size_t nHoles = 0;
  /// Sum of chi-squared values from measurements
  double chi2Sum = 0;
  /// Chi-squared contribution from each measurement
  std::vector<double> measurementChi2 = {};
  /// Chi-squared contribution from each outlier
  std::vector<double> outlierChi2 = {};
  /// Number of degrees of freedom in track fit
  std::size_t NDF = 0;
  /// Volume identifiers for measurement surfaces
  std::vector<unsigned int> measurementVolume = {};
  /// Layer identifiers for measurement surfaces
  std::vector<unsigned int> measurementLayer = {};
  /// Volume identifiers for outlier surfaces
  std::vector<unsigned int> outlierVolume = {};
  /// Layer identifiers for outlier surfaces
  std::vector<unsigned int> outlierLayer = {};
  /// Number of hits shared with other tracks
  std::size_t nSharedHits = 0;
};

// Container for trajectory summary info at a specific volume
using VolumeTrajectoryStateContainer =
    std::unordered_map<GeometryIdentifier::Value, TrajectoryState>;

/// @brief Getter for global trajectory info
///
/// @param multiTraj The MultiTrajectory object
/// @param tipIndex The entry index of trajectory to investigate
///
/// @return The trajectory summary info
template <typename traj_t>
TrajectoryState trajectoryState(const traj_t& multiTraj, std::size_t tipIndex) {
  TrajectoryState trajState;
  multiTraj.visitBackwards(tipIndex, [&](const auto& state) {
    // Get the volume Id of this surface
    const auto& geoID = state.referenceSurface().geometryId();
    const auto& volume = geoID.volume();
    const auto& layer = geoID.layer();
    trajState.nStates++;
    auto typeFlags = state.typeFlags();
    if (typeFlags.isHole()) {
      trajState.nHoles++;
    } else if (typeFlags.isOutlier()) {
      trajState.nOutliers++;
      trajState.outlierChi2.push_back(state.chi2());
      trajState.outlierVolume.push_back(volume);
      trajState.outlierLayer.push_back(layer);
    } else if (typeFlags.isMeasurement()) {
      if (typeFlags.isSharedHit()) {
        trajState.nSharedHits++;
      }
      trajState.nMeasurements++;
      trajState.measurementChi2.push_back(state.chi2());
      trajState.measurementVolume.push_back(volume);
      trajState.measurementLayer.push_back(layer);
      trajState.chi2Sum += state.chi2();
      trajState.NDF += state.calibratedSize();
    }
  });
  return trajState;
}

/// @brief Getter for trajectory info for different sub-detectors
///
/// @tparam source_link_t Type of source link
///
/// @param multiTraj The MultiTrajectory object
/// @param tipIndex The entry index of trajectory to investigate
/// track states at different sub-detectors.
/// @param volumeIds The container for sub-detector Ids
///
/// @return The trajectory summary info at different sub-detectors (i.e.
/// different volumes)
template <typename traj_t>
VolumeTrajectoryStateContainer trajectoryState(
    const traj_t& multiTraj, std::size_t tipIndex,
    const std::vector<GeometryIdentifier::Value>& volumeIds) {
  VolumeTrajectoryStateContainer trajStateContainer;
  multiTraj.visitBackwards(tipIndex, [&](const auto& state) {
    // Get the volume Id of this surface
    const auto& geoID = state.referenceSurface().geometryId();
    const auto& volume = geoID.volume();
    const auto& layer = geoID.layer();
    // Check if the track info for this sub-detector is requested
    if (!rangeContainsValue(volumeIds, volume)) {
      return true;
    }
    // The trajectory state for this volume
    auto& trajState = trajStateContainer[volume];
    trajState.nStates++;
    trajState.NDF += state.calibratedSize();
    auto typeFlags = state.typeFlags();
    if (typeFlags.isHole()) {
      trajState.nHoles++;
    } else if (typeFlags.isOutlier()) {
      trajState.nOutliers++;
      trajState.outlierChi2.push_back(state.chi2());
      trajState.outlierVolume.push_back(volume);
      trajState.outlierLayer.push_back(layer);
    } else if (typeFlags.isMeasurement()) {
      if (typeFlags.isSharedHit()) {
        trajState.nSharedHits++;
      }
      trajState.nMeasurements++;
      trajState.measurementChi2.push_back(state.chi2());
      trajState.measurementVolume.push_back(volume);
      trajState.measurementLayer.push_back(layer);
      trajState.chi2Sum += state.chi2();
    }
    return true;
  });
  return trajStateContainer;
}

/// @brief Transforms the filtered parameters from a @c TrackStateProxy to free
/// parameters
///
/// @tparam track_state_proxy_t Type of the @c TrackStateProxy
/// @param [in] gctx Geometry context
/// @param [in] trackStateProxy TrackStateProxy
///
/// @return Free parameters representation of the filtered parameters
template <typename track_state_proxy_t>
FreeVector freeFiltered(const GeometryContext& gctx,
                        const track_state_proxy_t& trackStateProxy) {
  return transformBoundToFreeParameters(trackStateProxy.referenceSurface(),
                                        gctx, trackStateProxy.filtered());
}

/// @brief Transforms the smoothed parameters from a @c TrackStateProxy to free
/// parameters
///
/// @tparam track_state_proxy_t Type of the @c TrackStateProxy
/// @param [in] gctx Geometry context
/// @param [in] trackStateProxy TrackStateProxy
///
/// @return Free parameters representation of the smoothed parameters
template <typename track_state_proxy_t>
FreeVector freeSmoothed(const GeometryContext& gctx,
                        const track_state_proxy_t& trackStateProxy) {
  return transformBoundToFreeParameters(trackStateProxy.referenceSurface(),
                                        gctx, trackStateProxy.smoothed());
}
}  // namespace Acts::MultiTrajectoryHelpers
