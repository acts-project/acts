// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/AnyTrackStateProxy.hpp"
#include "Acts/EventData/TrackProxyConcept.hpp"
#include "Acts/EventData/TrackStatePropMask.hpp"
#include "Acts/TrackFitting/KalmanFitterError.hpp"
#include "Acts/TrackFitting/detail/SmoothedConditioning.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"

#include <algorithm>
#include <iterator>
#include <stdexcept>
#include <vector>

namespace Acts {

/// @addtogroup track_fitting
/// @{

/// Reject a measurement of a smoothed track without a refit. The state
/// becomes an outlier and the smoothed states become those of a fit with this
/// outlier, exact in the linear model around the stored trajectory.
///
/// Works in place; copy a track whose states are shared with other tracks.
/// The track must come from one filter pass (Kalman fitter, CKF), not from
/// joined tracklets or the GSF. Only the smoothed states and the summary
/// change; the later filtered states keep the measurement, so a second
/// rejection on the same track is not exact. The chi2 drops by the unbiased
/// chi2 of the measurement.
///
/// @tparam track_proxy_t The track proxy type
/// @param track The track, smoothed
/// @param state The measurement state to reject
/// @param logger The logger
/// @return Success or the first error; on error the smoothed states may
///         have changed
template <TrackProxyConcept track_proxy_t>
Result<void> rejectMeasurement(track_proxy_t& track,
                               typename track_proxy_t::TrackStateProxy state,
                               const Logger& logger = getDummyLogger()) {
  if (&state.container() != &track.container().trackStateContainer()) {
    throw std::invalid_argument("rejectMeasurement: state is not on the track");
  }
  if (!state.typeFlags().isMeasurement()) {
    throw std::invalid_argument("rejectMeasurement: state is no measurement");
  }
  if (!state.hasSmoothed() || !state.hasPredicted()) {
    ACTS_DEBUG("State " << state.index()
                        << " lacks smoothed or predicted parameters");
    return KalmanFitterError::SmoothFailed;
  }

  std::vector<AnyMutableTrackStateProxy> states = detail::statesInFilterOrder(
      track.container().trackStateContainer(), track.tipIndex());
  const auto it = std::ranges::find_if(
      states, [&](const auto& ts) { return ts.index() == state.index(); });
  if (it == states.end()) {
    throw std::invalid_argument("rejectMeasurement: state is not on the track");
  }

  auto chi2 = detail::rejectSmoothedMeasurement(
      states, static_cast<std::size_t>(std::distance(states.begin(), it)),
      logger);
  if (!chi2.ok()) {
    return chi2.error();
  }

  state.typeFlags().setIsOutlier();
  state.shareFrom(TrackStatePropMask::Predicted, TrackStatePropMask::Filtered);
  track.chi2() -= static_cast<float>(*chi2);
  track.nDoF() -= state.calibratedSize();
  track.nMeasurements()--;
  track.nOutliers()++;
  return Result<void>::success();
}

/// @}

}  // namespace Acts
