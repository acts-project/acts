// This file is part of the Acts project.
//
// Copyright (C) 2018-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"
#include "Acts/Utilities/TypeTraits.hpp"

namespace Acts {

template <typename traj_t>
void voidKalmanCalibrator(
    const GeometryContext& /*gctx*/,
    typename MultiTrajectory<traj_t>::TrackStateProxy /*trackState*/) {
  throw std::runtime_error{"VoidKalmanCalibrator should not ever execute"};
}

template <typename traj_t>
Result<void> voidKalmanUpdater(
    const GeometryContext& /*gctx*/,
    typename MultiTrajectory<traj_t>::TrackStateProxy trackState,
    NavigationDirection /*direction*/, LoggerWrapper /*logger*/) {
  trackState.filtered() = trackState.predicted();
  trackState.filteredCovariance() = trackState.predictedCovariance();
  return Result<void>::success();
}

template <typename traj_t>
Result<void> voidKalmanSmoother(const GeometryContext& /*gctx*/,
                                MultiTrajectory<traj_t>& trackStates,
                                size_t entry, LoggerWrapper /*logger*/) {
  trackStates.applyBackwards(entry, [](auto trackState) {
    trackState.smoothed() = trackState.filtered();
    trackState.smoothedCovariance() = trackState.filteredCovariance();
  });

  return Result<void>::success();
}

template <typename traj_t>
bool voidOutlierFinder(
    typename MultiTrajectory<traj_t>::ConstTrackStateProxy /*trackState*/) {
  return false;
}

template <typename traj_t>
bool voidReverseFilteringLogic(
    typename MultiTrajectory<traj_t>::ConstTrackStateProxy /*trackState*/) {
  return false;
}

}  // namespace Acts
