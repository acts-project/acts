// This file is part of the Acts project.
//
// Copyright (C) 2018-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/detail/covariance_helper.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/TrackFitting/KalmanFitterError.hpp"
#include "Acts/Utilities/Delegate.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"

namespace Acts {

/// Kalman trajectory smoother based on gain matrix formalism.
///
/// This implements not a single smoothing step, but the full backwards
/// smoothing procedure for a filtered, forward trajectory using the stored
/// linearization.
class GainMatrixSmoother {
  using TrackStateTraits =
      TrackStateTraits<MultiTrajectoryTraits::MeasurementSizeMax, false>;

  // struct InternalTrackState {
  // TrackStateTraits::Parameters filtered;
  // TrackStateTraits::Covariance filteredCovariance;
  // TrackStateTraits::Parameters smoothed;
  // TrackStateTraits::Parameters prevPredicted;
  // TrackStateTraits::Covariance prevPredictedCovariance;
  // TrackStateTraits::Parameters prevSmoothed;
  // TrackStateTraits::Covariance prevSmoothedCovariance;
  // TrackStateTraits::Covariance prevJacobian;
  // };

 public:
  /// Run the Kalman smoothing for one trajectory.
  ///
  /// @param[in] gctx The geometry context for the smoothing
  /// @param[in,out] trajectory The trajectory to be smoothed
  /// @param[in] entryIndex The index of state to start the smoothing
  /// @param[in] logger Where to write logging information to
  template <typename D>
  Result<void> operator()(const GeometryContext& gctx,
                          MultiTrajectory<D>& trajectory, size_t entryIndex,
                          LoggerWrapper logger = getDummyLogger()) const {
    (void)gctx;

    using TrackStateProxy = typename MultiTrajectory<D>::TrackStateProxy;

    GetParameters filtered;
    GetCovariance filteredCovariance;
    GetParameters smoothed;
    GetParameters predicted;
    GetCovariance predictedCovariance;
    GetCovariance smoothedCovariance;
    GetCovariance jacobian;

    filtered.connect([](const void*, void* ts) {
      return static_cast<TrackStateProxy*>(ts)->filtered();
    });
    filteredCovariance.connect([](const void*, void* ts) {
      return static_cast<TrackStateProxy*>(ts)->filteredCovariance();
    });

    smoothed.connect([](const void*, void* ts) {
      return static_cast<TrackStateProxy*>(ts)->smoothed();
    });
    smoothedCovariance.connect([](const void*, void* ts) {
      return static_cast<TrackStateProxy*>(ts)->smoothedCovariance();
    });

    predicted.connect([](const void*, void* ts) {
      return static_cast<TrackStateProxy*>(ts)->predicted();
    });
    predictedCovariance.connect([](const void*, void* ts) {
      return static_cast<TrackStateProxy*>(ts)->predictedCovariance();
    });

    jacobian.connect([](const void*, void* ts) {
      return static_cast<TrackStateProxy*>(ts)->jacobian();
    });

    // return calculate(entryIndex, filtered, filteredCovariance, smoothed,
    // predicted, predictedCovariance, smoothedCovariance,
    // jacobian, logger);

    ACTS_VERBOSE("Invoked GainMatrixSmoother on entry index: " << entryIndex);

    // For the last state: smoothed is filtered - also: switch to next
    ACTS_VERBOSE("Getting previous track state");
    auto prev_ts = trajectory.getTrackState(entryIndex);

    prev_ts.smoothed() = prev_ts.filtered();
    prev_ts.smoothedCovariance() = prev_ts.filteredCovariance();

    // make sure there is more than one track state
    if (!prev_ts.hasPrevious()) {
      ACTS_VERBOSE("Only one track state given, smoothing terminates early");
      return Result<void>::success();
    }

    ACTS_VERBOSE("Start smoothing from previous track state at index: "
                 << prev_ts.previous());

    // @TODO: Put back into own compilation unit!

    // default-constructed error represents success, i.e. an invalid error code
    std::error_code error;
    trajectory.applyBackwards(prev_ts.previous(), [&, this](auto ts) {
      // should have filtered and predicted, this should also include the
      // covariances.
      assert(ts.hasFiltered());
      assert(ts.hasPredicted());
      assert(ts.hasJacobian());

      // previous trackstate should have smoothed and predicted
      assert(prev_ts.hasSmoothed());
      assert(prev_ts.hasPredicted());

      ACTS_VERBOSE("Calculate smoothing matrix:");
      ACTS_VERBOSE("Filtered covariance:\n" << ts.filteredCovariance());
      ACTS_VERBOSE("Jacobian:\n" << ts.jacobian());
      ACTS_VERBOSE("Prev. predicted covariance\n"
                   << prev_ts.predictedCovariance() << "\n, inverse: \n"
                   << prev_ts.predictedCovariance().inverse());

      // if (auto res = calculate(
      // InternalTrackState{
      // ts.filtered(), ts.filteredCovariance(), ts.smoothed(),
      // prev_ts.predicted(), prev_ts.predictedCovariance(),
      // prev_ts.smoothed(), prev_ts.smoothedCovariance(),
      // prev_ts.jacobian()},
      // logger);
      // !res.ok()) {
      // error = res.error();
      // return false;
      // }

      if (auto res = calculate(&ts, &prev_ts, filtered, filteredCovariance,
                               smoothed, predicted, predictedCovariance,
                               smoothedCovariance, jacobian, logger);
          !res.ok()) {
        error = res.error();
        return false;
      }

      prev_ts = ts;
      return true;  // continue execution
    });

    return error ? Result<void>::failure(error) : Result<void>::success();
  }

  using GetParameters = Acts::Delegate<TrackStateTraits::Parameters(void*)>;
  using GetCovariance = Acts::Delegate<TrackStateTraits::Covariance(void*)>;

  Result<void> calculate(void* ts, void* prev_ts, const GetParameters& filtered,
                         const GetCovariance& filteredCovariance,
                         const GetParameters& smoothed,
                         const GetParameters& predicted,
                         const GetCovariance& predictedCovariance,
                         const GetCovariance& smoothedCovariance,
                         const GetCovariance& jacobian,
                         LoggerWrapper logger) const;

  // Result<void> calculate(InternalTrackState trackState,
  // LoggerWrapper logger) const {
  // // Gain smoothing matrix
  // // NB: The jacobian stored in a state is the jacobian from previous
  // // state to this state in forward propagation
  // BoundMatrix G = trackState.filteredCovariance *
  // trackState.prevJacobian.transpose() *
  // trackState.prevPredictedCovariance.inverse();

  // if (G.hasNaN()) {
  // return KalmanFitterError::SmoothFailed;
  // }

  // ACTS_VERBOSE("Gain smoothing matrix G:\n" << G);

  // ACTS_VERBOSE("Calculate smoothed parameters:");
  // ACTS_VERBOSE("Filtered parameters: " << trackState.filtered.transpose());
  // ACTS_VERBOSE(
  // "Prev. smoothed parameters: " << trackState.prevSmoothed.transpose());
  // ACTS_VERBOSE(
  // "Prev. predicted parameters: " << trackState.prevPredicted.transpose());

  // // Calculate the smoothed parameters
  // trackState.smoothed = trackState.filtered + G * (trackState.prevSmoothed -
  // trackState.prevPredicted);

  // ACTS_VERBOSE(
  // "Smoothed parameters are: " << trackState.smoothed.transpose());
  // ACTS_VERBOSE("Calculate smoothed covariance:");
  // ACTS_VERBOSE("Prev. smoothed covariance:\n"
  // << trackState.prevSmoothedCovariance);

  // // And the smoothed covariance
  // trackState.prevSmoothedCovariance =
  // trackState.filteredCovariance -
  // G *
  // (trackState.prevPredictedCovariance -
  // trackState.prevSmoothedCovariance) *
  // G.transpose();

  // return Result<void>::success();
  // }
};

}  // namespace Acts
