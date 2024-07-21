// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/TrackFitting/KalmanFitterError.hpp"
#include "Acts/Utilities/Delegate.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"

#include <cassert>
#include <cstddef>
#include <system_error>

namespace Acts {

/// Kalman trajectory smoother based on the Modified Bryson–Frazier smoother.
///
/// This implements not a single smoothing step, but the full backwards
/// smoothing procedure for a filtered, forward trajectory using the stored
/// linearization.
class MbfSmoother {
 public:
  struct InternalTrackState {
    // This is used to build a covariance matrix view in the .cpp file
    unsigned int calibratedSize;
    const double* calibrated;
    const double* calibratedCovariance;
    TrackStateTraits<MultiTrajectoryTraits::MeasurementSizeMax,
                     false>::Projector projector;
    TrackStateTraits<MultiTrajectoryTraits::MeasurementSizeMax,
                     false>::Covariance jacobian;

    TrackStateTraits<MultiTrajectoryTraits::MeasurementSizeMax,
                     false>::Parameters predicted;
    TrackStateTraits<MultiTrajectoryTraits::MeasurementSizeMax,
                     false>::Covariance predictedCovariance;
    TrackStateTraits<MultiTrajectoryTraits::MeasurementSizeMax,
                     false>::Parameters filtered;
    TrackStateTraits<MultiTrajectoryTraits::MeasurementSizeMax,
                     false>::Covariance filteredCovariance;
    TrackStateTraits<MultiTrajectoryTraits::MeasurementSizeMax,
                     false>::Parameters smoothed;
    TrackStateTraits<MultiTrajectoryTraits::MeasurementSizeMax,
                     false>::Covariance smoothedCovariance;
  };

  /// Run the Kalman smoothing for one trajectory.
  ///
  /// @param[in,out] trajectory The trajectory to be smoothed
  /// @param[in] entryIndex The index of state to start the smoothing
  /// @param[in] logger Where to write logging information to
  template <typename traj_t>
  Result<void> operator()(const GeometryContext& /*gctx*/, traj_t& trajectory,
                          std::size_t entryIndex,
                          const Logger& /*logger*/ = getDummyLogger()) const {
    using TrackStateProxy = typename traj_t::TrackStateProxy;

    TrackStateProxy start_ts = trajectory.getTrackState(entryIndex);

    BoundMatrix big_lambda_hat = BoundMatrix::Zero();
    BoundVector small_lambda_hat = BoundVector::Zero();

    trajectory.applyBackwards(start_ts.index(), [&](TrackStateProxy ts) {
      // ensure the track state has a smoothed component
      ts.addComponents(TrackStatePropMask::Smoothed);

      InternalTrackState internalTrackState{
          ts.calibratedSize(),
          // Note that we pass raw pointers here which are used in the correct
          // shape later
          ts.effectiveCalibrated().data(),
          ts.effectiveCalibratedCovariance().data(),
          ts.projector(),
          ts.jacobian(),
          ts.predicted(),
          ts.predictedCovariance(),
          ts.filtered(),
          ts.filteredCovariance(),
          ts.smoothed(),
          ts.smoothedCovariance(),
      };

      calculateSmoothed(internalTrackState, big_lambda_hat, small_lambda_hat);

      if (!ts.hasPrevious()) {
        return;
      }

      if (!ts.typeFlags().test(TrackStateFlag::MeasurementFlag)) {
        visitNonMeasurement(internalTrackState, big_lambda_hat,
                            small_lambda_hat);

        return;
      }

      visitMeasurement(internalTrackState, big_lambda_hat, small_lambda_hat);
    });

    return Result<void>::success();
  }

 private:
  void calculateSmoothed(InternalTrackState& ts,
                         const BoundMatrix& big_lambda_hat,
                         const BoundVector& small_lambda_hat) const;
  void visitNonMeasurement(const InternalTrackState& ts,
                           BoundMatrix& big_lambda_hat,
                           BoundVector& small_lambda_hat) const;
  void visitMeasurement(const InternalTrackState& ts,
                        BoundMatrix& big_lambda_hat,
                        BoundVector& small_lambda_hat) const;
};

}  // namespace Acts
