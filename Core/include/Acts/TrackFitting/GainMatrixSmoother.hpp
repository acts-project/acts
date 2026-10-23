// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/AnyTrackStateProxy.hpp"
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/Types.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/TrackFitting/KalmanFitterError.hpp"
#include "Acts/TrackFitting/SmootherAdjoint.hpp"
#include "Acts/Utilities/Delegate.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"

#include <cassert>
#include <cstddef>
#include <system_error>

namespace Acts {

/// Kalman trajectory smoother based on gain matrix formalism.
///
/// This implements not a single smoothing step, but the full backwards
/// smoothing procedure for a filtered, forward trajectory using the stored
/// linearization.
/// @ingroup track_fitting
class GainMatrixSmoother {
 public:
  /// Whether to check the covariance matrices if they are semi-positive and if
  /// not attempt to correct them.
  bool doCovCheckAndAttemptFix = false;

  /// Run the Kalman smoothing for one trajectory.
  ///
  /// @param[in] gctx The geometry context to be used
  /// @param[in,out] trajectory The trajectory to be smoothed
  /// @param[in] entryIndex The index of state to start the smoothing
  /// @param[in] logger Where to write logging information to
  /// @return Success or failure of the smoothing procedure
  template <typename traj_t>
  Result<void> operator()(const GeometryContext& gctx, traj_t& trajectory,
                          std::size_t entryIndex,
                          const Logger& logger = getDummyLogger()) const {
    return smooth(gctx, trajectory, entryIndex, kTrackIndexInvalid, nullptr,
                  logger);
  }

  /// Run the Kalman smoothing for a segment of one trajectory and carry the
  /// adjoint out.
  ///
  /// The walk starts at @p entryIndex and follows the previous links. With a
  /// valid @p stopIndex it smooths that state last; otherwise it smooths to
  /// the end of the chain. With @p carryOut the gain-matrix form of the
  /// adjoint after the last smoothed state is transported through that
  /// state's jacobian and returned, so it lives on the node the last state
  /// was predicted from. Apply it with @ref applySmootherAdjoint. The carry
  /// is the same quantity the @ref MbfSmoother returns, so the two smoothers
  /// stay interchangeable.
  ///
  /// @tparam traj_t The MultiTrajectory type
  /// @param gctx The geometry context to be used
  /// @param trajectory The trajectory to be smoothed
  /// @param entryIndex The index of state to start the smoothing
  /// @param stopIndex The index of the last state to smooth, or
  ///        @ref kTrackIndexInvalid for the end of the chain
  /// @param carryOut The adjoint after the last smoothed state, or nullptr
  /// @param logger Where to write logging information to
  /// @return The status of the smoothing
  template <typename traj_t>
  Result<void> smooth(const GeometryContext& gctx, traj_t& trajectory,
                      std::size_t entryIndex, std::size_t stopIndex,
                      SmootherAdjoint* carryOut,
                      const Logger& logger = getDummyLogger()) const {
    static_cast<void>(gctx);

    ACTS_VERBOSE("Invoked GainMatrixSmoother on entry index: " << entryIndex);

    // For the last state: smoothed is filtered - also: switch to next
    ACTS_VERBOSE("Getting previous track state");
    auto prev_ts = trajectory.getTrackState(entryIndex);

    prev_ts.shareFrom(TrackStatePropMask::Filtered,
                      TrackStatePropMask::Smoothed);

    bool reachedStop = (prev_ts.index() == stopIndex);

    // make sure there is more than one track state
    if (!prev_ts.hasPrevious() || reachedStop) {
      ACTS_VERBOSE("Only one track state given, smoothing terminates early");
    } else {
      ACTS_VERBOSE("Start smoothing from previous track state at index: "
                   << prev_ts.previous());

      // default-constructed error represents success, i.e. an invalid error
      // code
      std::error_code error;
      trajectory.applyBackwards(prev_ts.previous(), [&, this](auto ts) {
        // should have filtered and predicted, this should also include the
        // covariances.
        assert(ts.hasFiltered());
        assert(ts.hasPredicted());

        // previous trackstate should have smoothed and predicted
        assert(prev_ts.hasSmoothed());
        assert(prev_ts.hasPredicted());
        assert(prev_ts.hasJacobian());

        ACTS_VERBOSE("Calculate smoothing matrix:");
        ACTS_VERBOSE("Filtered covariance:\n" << ts.filteredCovariance());
        ACTS_VERBOSE("Jacobian:\n" << prev_ts.jacobian());

        // ensure the track state has a smoothed component
        ts.addComponents(TrackStatePropMask::Smoothed);

        if (auto res = calculate(AnyMutableTrackStateProxy{ts},
                                 AnyConstTrackStateProxy{prev_ts}, logger);
            !res.ok()) {
          error = res.error();
          return false;
        }

        prev_ts = ts;
        reachedStop = (ts.index() == stopIndex);
        return !reachedStop;  // continue execution
      });

      if (error) {
        return Result<void>::failure(error);
      }
    }

    if (stopIndex != kTrackIndexInvalid && !reachedStop) {
      return KalmanFitterError::SmoothFailed;
    }

    if (carryOut != nullptr) {
      auto adjoint = calculateAdjoint(AnyConstTrackStateProxy{prev_ts}, logger);
      if (!adjoint.ok()) {
        return adjoint.error();
      }
      *carryOut = *adjoint;
    }

    return Result<void>::success();
  }

  using GetParameters =
      Acts::Delegate<TrackStateTraits<kMeasurementSizeMax, false>::Parameters(
          void*)>;
  /// Type alias for delegate to get track state covariance matrix
  using GetCovariance =
      Acts::Delegate<TrackStateTraits<kMeasurementSizeMax, false>::Covariance(
          void*)>;

  /// Calculate smoothed parameters for a single track state using gain matrix
  /// formalism.
  ///
  /// @param ts Current track state to be smoothed
  /// @param prev_ts Previous track state (in forward direction)
  /// @param logger Logger for verbose output
  /// @return Success or failure of the smoothing calculation
  Result<void> calculate(AnyMutableTrackStateProxy ts,
                         AnyConstTrackStateProxy prev_ts,
                         const Logger& logger) const;

  /// Calculate the adjoint on the node a smoothed state was predicted from.
  ///
  /// With the predicted estimate `(xp, Pp)`, the smoothed estimate `(xs, Ps)`
  /// and the jacobian `F` of the state, the adjoint is
  /// `bigLambda = F^T Pp^-1 (Pp - Ps) Pp^-1 F` and
  /// `smallLambda = -F^T Pp^-1 (xs - xp)`. Applied to the estimate of the
  /// node with @ref applySmootherAdjoint this is the Rauch-Tung-Striebel
  /// step onto that node.
  ///
  /// @param ts The last smoothed state, with predicted, smoothed and jacobian
  /// @param logger The logger
  /// @return The adjoint on the node the state was predicted from, or
  ///         @ref KalmanFitterError::SmoothFailed if the predicted covariance
  ///         is not positive definite or the result is not finite
  Result<SmootherAdjoint> calculateAdjoint(AnyConstTrackStateProxy ts,
                                           const Logger& logger) const;
};

}  // namespace Acts
