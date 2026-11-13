// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/AnyTrackStateProxy.hpp"
#include "Acts/EventData/SubspaceHelpers.hpp"
#include "Acts/EventData/TrackProxyConcept.hpp"
#include "Acts/EventData/TrackStatePropMask.hpp"
#include "Acts/TrackFitting/KalmanFitterError.hpp"
#include "Acts/TrackFitting/detail/SmoothedConditioning.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"

#include <algorithm>
#include <array>
#include <iterator>
#include <span>
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

/// Remove the prior of the filter from a smoothed track without a refit, so
/// that the smoothed states depend on the measurements alone. The prior, the
/// predicted state of the first state, is removed like a measurement of
/// @p parameters. Parameters without measurement information, usually the
/// time, must stay out; their prior is kept, conditional on the others.
///
/// Works in place, with the requirements of @ref rejectMeasurement. The chi2
/// drops by the unbiased chi2 of the prior, the degrees of freedom by the
/// number of parameters.
///
/// @tparam track_proxy_t The track proxy type
/// @param track The track, smoothed
/// @param parameters The parameters whose prior is removed
/// @param logger The logger
/// @return Success or the first error; on error the smoothed states may
///         have changed
template <TrackProxyConcept track_proxy_t>
Result<void> removePrior(track_proxy_t& track,
                         std::span<const BoundIndices> parameters,
                         const Logger& logger = getDummyLogger()) {
  if (parameters.empty() ||
      !checkSubspaceIndices(parameters, eBoundSize, parameters.size())) {
    throw std::invalid_argument("removePrior: invalid parameters");
  }
  if (track.nDoF() < parameters.size()) {
    ACTS_DEBUG("Too few degrees of freedom to remove the prior");
    return KalmanFitterError::SmoothFailed;
  }
  std::vector<AnyMutableTrackStateProxy> states = detail::statesInFilterOrder(
      track.container().trackStateContainer(), track.tipIndex());
  if (states.empty() || !states.front().hasSmoothed() ||
      !states.front().hasPredicted()) {
    ACTS_DEBUG("The first state lacks smoothed or predicted parameters");
    return KalmanFitterError::SmoothFailed;
  }

  auto chi2 = detail::removeSmoothedPrior(states, parameters, logger);
  if (!chi2.ok()) {
    return chi2.error();
  }
  track.chi2() -= static_cast<float>(*chi2);
  track.nDoF() -= static_cast<unsigned int>(parameters.size());
  return Result<void>::success();
}

/// Remove the prior of the filter on all parameters but the time.
/// @see removePrior
/// @tparam track_proxy_t The track proxy type
/// @param track The track, smoothed
/// @param logger The logger
/// @return Success or the first error
template <TrackProxyConcept track_proxy_t>
Result<void> removePrior(track_proxy_t& track,
                         const Logger& logger = getDummyLogger()) {
  static constexpr std::array parameters = {eBoundLoc0, eBoundLoc1, eBoundPhi,
                                            eBoundTheta, eBoundQOverP};
  return removePrior(track, parameters, logger);
}

/// @}

}  // namespace Acts
