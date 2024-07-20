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
class MBFSmoother {
 public:
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

      ts.smoothedCovariance() =
          ts.filteredCovariance() -
          ts.filteredCovariance() * big_lambda_hat * ts.filteredCovariance();
      ts.smoothed() =
          ts.filtered() - ts.filteredCovariance() * small_lambda_hat;

      if (!ts.hasPrevious()) {
        return;
      }

      if (!ts.typeFlags().test(TrackStateFlag::MeasurementFlag)) {
        const auto F = ts.jacobian();

        big_lambda_hat = F.transpose() * big_lambda_hat * F;
        small_lambda_hat = F.transpose() * small_lambda_hat;

        return;
      }

      visit_measurement(ts.calibratedSize(), [&](auto N) -> void {
        constexpr std::size_t kMeasurementSize = decltype(N)::value;

        const auto H =
            ts.projector()
                .template topLeftCorner<kMeasurementSize, eBoundSize>()
                .eval();

        const auto S = (H * ts.predictedCovariance() * H.transpose() +
                        ts.template calibratedCovariance<kMeasurementSize>())
                           .eval();
        // TODO Sinv could be cached by the filter step
        const auto Sinv = S.inverse().eval();

        // TODO K could be cached by the filter step
        const auto K = (ts.predictedCovariance() * H.transpose() * Sinv).eval();

        const auto C = (BoundMatrix::Identity() - K * H).eval();
        const auto y =
            (ts.template calibrated<kMeasurementSize>() - H * ts.predicted())
                .eval();

        const auto big_lambda_tilde =
            (H.transpose() * Sinv * H + C.transpose() * big_lambda_hat * C)
                .eval();
        const auto small_lambda_tilde =
            (-H.transpose() * Sinv * y + C.transpose() * small_lambda_hat)
                .eval();

        const auto F = ts.jacobian();

        big_lambda_hat = F.transpose() * big_lambda_tilde * F;
        small_lambda_hat = F.transpose() * small_lambda_tilde;
      });
    });

    return Result<void>::success();
  }
};

}  // namespace Acts
