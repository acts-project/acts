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

using SubDetectorClassifier = std::function<bool(const Acts::GeometryID&)>;
using SubDetectorClassifierContainer =
    std::unordered_map<std::string, SubDetectorClassifier>;
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
/// @param subDetClassifierContainer The container for classifiers to classfify
/// track states at different sub-detectors.
///
/// @return The trajectory summary info at different sub-detectors
template <typename source_link_t>
TrajectoryStateContainer trajectoryState(
    const Acts::MultiTrajectory<source_link_t>& multiTraj,
    const size_t& entryIndex,
    const SubDetectorClassifierContainer& subDetClassifierContainer) {
  TrajectoryStateContainer trajStateContainer;
  multiTraj.visitBackwards(entryIndex, [&](const auto& state) {
    // Get the geometry identifier of the reference surface
    auto geoID = state.referenceSurface().geoID();
    // Loop over all the classifiers
    for (const auto& [detName, detClassifier] : subDetClassifierContainer) {
      if (not detClassifier(geoID)) {
        continue;
      }
      // The trajectory state corresponding to this classifier
      auto& trajState = trajStateContainer[detName];
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
    }  // end of loop for the classifiers
  });
  return trajStateContainer;
}

}  // namespace MultiTrajectoryHelpers

}  // namespace Acts
