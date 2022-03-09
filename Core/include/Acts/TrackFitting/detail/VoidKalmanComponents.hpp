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

inline void voidKalmanCalibrator(
    const GeometryContext& /*gctx*/,
    MultiTrajectory::TrackStateProxy /*trackState*/) {
  throw std::runtime_error{"VoidKalmanCalibrator should not ever execute"};
}

inline Result<void> voidKalmanUpdater(
    const GeometryContext& /*gctx*/,
    MultiTrajectory::TrackStateProxy trackState,
    NavigationDirection /*direction*/, LoggerWrapper /*logger*/) {
  trackState.filtered() = trackState.predicted();
  trackState.filteredCovariance() = trackState.predictedCovariance();
  return Result<void>::success();
}

inline Result<void> voidKalmanSmoother(const GeometryContext& /*gctx*/,
                                       MultiTrajectory& trackStates,
                                       size_t entry, LoggerWrapper /*logger*/) {
  trackStates.applyBackwards(entry, [](const auto trackState) {
    trackState.smoothed() = trackState.filtered();
    trackState.smoothedCovariance() = trackState.filteredCovariance();
  });

  return Result<void>::success();
}

inline bool voidOutlierFinder(
    MultiTrajectory::ConstTrackStateProxy /*trackState*/) {
  return false;
}

inline bool voidReverseFilteringLogic(
    MultiTrajectory::ConstTrackStateProxy /*trackState*/) {
  return false;
}

}  // namespace Acts
