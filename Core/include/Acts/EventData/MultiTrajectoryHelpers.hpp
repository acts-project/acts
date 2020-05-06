// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
#include <functional>
#include <unordered_map>

#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/Geometry/Layer.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Surfaces/Surface.hpp"

namespace Acts {

namespace MultiTrajectoryHelpers {

/// @brief Struct for brief trajectory summary info
/// @TODO: add nSharedHits
///
struct TrajectoryState {
  size_t nStates = 0;
  size_t nMeasurements = 0;
  size_t nOutliers = 0;
  size_t nHoles = 0;
  double chi2Sum = 0;
};

using TrajectoryStateContainer =
    std::unordered_map<std::string, TrajectoryState>;

/// @brief Getter for global trajectory info
///
/// @tparam source_link_t Type of source link
///
/// @param multiTraj The MultiTrajectory object
/// @param entryIndex The entry index of trajectory to investigate
///
/// @return The trajectory summary info
template <typename source_link_t>
TrajectoryState trajectoryState(
    const Acts::MultiTrajectory<source_link_t>& multiTraj,
    const size_t& entryIndex) {
  TrajectoryState trajState;
  multiTraj.visitBackwards(entryIndex, [&](const auto& state) {
    trajState.nStates++;
    trajState.chi2Sum += state.chi2();
    auto typeFlags = state.typeFlags();
    if (typeFlags.test(Acts::TrackStateFlag::MeasurementFlag)) {
      trajState.nMeasurements++;
    } else if (typeFlags.test(Acts::TrackStateFlag::OutlierFlag)) {
      trajState.nOutliers++;
    } else if (typeFlags.test(Acts::TrackStateFlag::HoleFlag)) {
      trajState.nHoles++;
    }
  });
  return trajState;
}

/// @brief Getter for trajectory info for different sub-detectors
///
/// @tparam source_link_t Type of source link
///
/// @param multiTraj The MultiTrajectory object
/// @param entryIndex The entry index of trajectory to investigate
/// @param subDetName The container for sub-detector names
/// track states at different sub-detectors.
///
/// @return The trajectory summary info at different sub-detectors
template <typename source_link_t>
TrajectoryStateContainer trajectoryState(
    const Acts::MultiTrajectory<source_link_t>& multiTraj,
    const size_t& entryIndex, const std::vector<std::string>& subDetName) {
  TrajectoryStateContainer trajStateContainer;
  multiTraj.visitBackwards(entryIndex, [&](const auto& state) {
    // Get the tracking volume name this surface is associated with
    const Surface* surface = &state.referenceSurface();
    std::string volumeName;
    if (surface->associatedLayer() != nullptr) {
      const Layer* layer = surface->associatedLayer();
      if (layer->trackingVolume() != nullptr) {
        volumeName = layer->trackingVolume()->volumeName();
      }
    }
    if (volumeName.empty()) {
      // Skip the state if the volume name is not found
      return true;
    }
    auto it = std::find(subDetName.begin(), subDetName.end(), volumeName);
    if (it == subDetName.end()) {
      // Skip the state if track info for this sub-detector is not requested
      return true;
    }
    // The trajectory state corresponding to this classifier
    auto& trajState = trajStateContainer[volumeName];
    trajState.nStates++;
    trajState.chi2Sum += state.chi2();
    auto typeFlags = state.typeFlags();
    if (typeFlags.test(Acts::TrackStateFlag::MeasurementFlag)) {
      trajState.nMeasurements++;
    } else if (typeFlags.test(Acts::TrackStateFlag::OutlierFlag)) {
      trajState.nOutliers++;
    } else if (typeFlags.test(Acts::TrackStateFlag::HoleFlag)) {
      trajState.nHoles++;
    }
    return true;
  });
  return trajStateContainer;
}

}  // namespace MultiTrajectoryHelpers

}  // namespace Acts
