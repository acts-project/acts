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
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"

#include <cstddef>
#include <system_error>
#include <utility>

namespace Acts {

/// Kalman trajectory smoother based on the Modified Bryson–Frazier (mBF)
/// smoother.
///
/// @ingroup track_fitting
///
/// The benefit of the mBF smoother is that it does not require the inverse of
/// the full covariance matrix, but only the inverse of the residual covariance
/// matrix which can be cached by the filter step. The same holds for the
/// Kalman gain matrix.
///
/// This implements not a single smoothing step, but the full backwards
/// smoothing procedure for a filtered, forward trajectory using the stored
/// linearization. Coordinate residuals and covariance canonicalization use
/// each state's surface parametrization and the supplied geometry context.
/// If no reference surface is attached (e.g. in standalone test states), the
/// interpretation is canonical global spherical. This is independent of the
/// state's position or proximity to the beam. No near-axis cutoff is applied;
/// covariance canonicalization is undefined exactly at theta = 0 or pi.
/// Traversal runs backwards from entryIndex. On failure, states later along
/// forward propagation may already have smoothed results; earlier states are
/// not visited.
///
/// See
/// [Wikipedia](https://en.wikipedia.org/wiki/Kalman_filter#Modified_Bryson%E2%80%93Frazier_smoother)
/// for more information.
class MbfSmoother {
 public:
  /// Run the Kalman smoothing for one trajectory.
  ///
  /// @param[in] gctx The geometry context to be used
  /// @param[in,out] trajectory The trajectory to be smoothed
  /// @param[in] entryIndex The index of state to start the smoothing
  /// @param[in] logger Where to write logging information to
  /// @return Success or failure of the MBF smoothing procedure
  template <typename traj_t>
  Result<void> operator()(const GeometryContext& gctx, traj_t& trajectory,
                          std::size_t entryIndex,
                          const Logger& logger = getDummyLogger()) const {
    static_cast<void>(logger);

    using TrackStateProxy = typename traj_t::TrackStateProxy;

    TrackStateProxy startTs = trajectory.getTrackState(entryIndex);

    // Notation consistent with the Wikipedia article
    // https://en.wikipedia.org/wiki/Kalman_filter
    BoundMatrix bigLambdaHat = BoundMatrix::Zero();
    BoundVector smallLambdaHat = BoundVector::Zero();

    std::error_code error;
    trajectory.applyBackwards(startTs.index(), [&](TrackStateProxy ts) {
      // ensure the track state has a smoothed component
      ts.addComponents(TrackStatePropMask::Smoothed);

      AnyMutableTrackStateProxy internalTrackState(ts);

      // Smooth the current state before transporting backward information.
      const auto smoothed = calculateSmoothed(gctx, internalTrackState,
                                              bigLambdaHat, smallLambdaHat);
      if (!smoothed.ok()) {
        error = smoothed.error();
        return false;
      }

      // We smoothed the last state - no need to update the lambdas
      if (!ts.hasPrevious()) {
        return false;
      }

      // Update the lambdas depending on the type of track state
      auto transported = Result<void>::success();
      if (ts.typeFlags().isMeasurement()) {
        transported = visitMeasurement(gctx, AnyConstTrackStateProxy{ts},
                                       bigLambdaHat, smallLambdaHat);
      } else {
        transported = visitNonMeasurement(std::as_const(ts).jacobian(),
                                          bigLambdaHat, smallLambdaHat);
      }
      if (!transported.ok()) {
        error = transported.error();
        return false;
      }
      return true;
    });

    return error ? Result<void>::failure(error) : Result<void>::success();
  }

 private:
  /// Calculate the smoothed parameters and covariance.
  Result<void> calculateSmoothed(const GeometryContext& gctx,
                                 AnyMutableTrackStateProxy& ts,
                                 const BoundMatrix& bigLambdaHat,
                                 const BoundVector& smallLambdaHat) const;

  /// Visit a non-measurement track state and update the lambdas.
  Result<void> visitNonMeasurement(
      const AnyConstTrackStateProxy::ConstCovarianceMap& jacobian,
      BoundMatrix& bigLambdaHat, BoundVector& smallLambdaHat) const;

  /// Visit a measurement track state and update the lambdas.
  Result<void> visitMeasurement(const GeometryContext& gctx,
                                const AnyConstTrackStateProxy& ts,
                                BoundMatrix& bigLambdaHat,
                                BoundVector& smallLambdaHat) const;

  template <std::size_t N>
  Result<void> visitMeasurementImpl(const GeometryContext& gctx,
                                    const AnyConstTrackStateProxy& ts,
                                    BoundMatrix& bigLambdaHat,
                                    BoundVector& smallLambdaHat) const;
};

}  // namespace Acts
