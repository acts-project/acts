// This file is part of the Acts project.
//
// Copyright (C) 2018-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Direction.hpp"
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Utilities/CalibrationContext.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"
#include "Acts/Utilities/TypeTraits.hpp"

namespace Acts::detail {

template <typename traj_t>
void voidFitterCalibrator(const GeometryContext& /*gctx*/,
                          const CalibrationContext& /*cctx*/,
                          const SourceLink& /*sourceLink*/,
                          typename traj_t::TrackStateProxy /*trackState*/) {
  throw std::runtime_error{"voidFitterCalibrator should not ever execute"};
}

template <typename traj_t>
Result<void> voidFitterUpdater(const GeometryContext& /*gctx*/,
                               typename traj_t::TrackStateProxy trackState,
                               Direction /*direction*/,
                               const Logger& /*logger*/) {
  trackState.filtered() = trackState.predicted();
  trackState.filteredCovariance() = trackState.predictedCovariance();
  return Result<void>::success();
}

template <typename traj_t>
Result<void> voidFitterSmoother(const GeometryContext& /*gctx*/,
                                traj_t& trackStates, size_t entry,
                                const Logger& /*logger*/) {
  trackStates.applyBackwards(entry, [](auto trackState) {
    trackState.smoothed() = trackState.filtered();
    trackState.smoothedCovariance() = trackState.filteredCovariance();
  });

  return Result<void>::success();
}

template <typename traj_t>
bool voidOutlierFinder(typename traj_t::ConstTrackStateProxy /*trackState*/) {
  return false;
}

template <typename traj_t>
bool voidReverseFilteringLogic(
    typename traj_t::ConstTrackStateProxy /*trackState*/) {
  return false;
}

inline const Surface* voidSurfaceAccessor(const SourceLink& /*sourceLink*/) {
  throw std::runtime_error{"voidSurfaceAccessor should not ever execute"};
  return nullptr;
}
}  // namespace Acts::detail
