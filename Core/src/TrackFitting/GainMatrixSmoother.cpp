// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/TrackFitting/GainMatrixSmoother.hpp"

namespace Acts {

Result<void> GainMatrixSmoother::calculate(
    void* ts, void* prev_ts, const GetParameters& filtered,
    const GetCovariance& filteredCovariance, const GetParameters& smoothed,
    const GetParameters& predicted, const GetCovariance& predictedCovariance,
    const GetCovariance& smoothedCovariance, const GetCovariance& jacobian,
    LoggerWrapper logger) const {
  static constexpr double epsilon = 1e-13;
  auto regularization = BoundMatrix::Identity() * epsilon;

  ACTS_VERBOSE("Prev. predicted covariance\n"
               << predictedCovariance(prev_ts) << "\n, inverse: \n"
               << predictedCovariance(prev_ts).inverse()
               << "\n, regularized inverse: \n"
               << (predictedCovariance(prev_ts) + regularization).inverse());

  // Gain smoothing matrix
  // NB: The jacobian stored in a state is the jacobian from previous
  // state to this state in forward propagation
  BoundMatrix G = filteredCovariance(ts) * jacobian(prev_ts).transpose() *
                  (predictedCovariance(prev_ts) + regularization).inverse();

  if (G.hasNaN()) {
    // error = KalmanFitterError::SmoothFailed;  // set to error
    // return false;                             // abort execution
    return KalmanFitterError::SmoothFailed;
  }

  ACTS_VERBOSE("Gain smoothing matrix G:\n" << G);

  ACTS_VERBOSE("Calculate smoothed parameters:");
  ACTS_VERBOSE("Filtered parameters: " << filtered(ts).transpose());
  ACTS_VERBOSE("Prev. smoothed parameters: " << smoothed(prev_ts).transpose());
  ACTS_VERBOSE(
      "Prev. predicted parameters: " << predicted(prev_ts).transpose());

  // Calculate the smoothed parameters
  smoothed(ts) = filtered(ts) + G * (smoothed(prev_ts) - predicted(prev_ts));

  ACTS_VERBOSE("Smoothed parameters are: " << smoothed(ts).transpose());
  ACTS_VERBOSE("Calculate smoothed covariance:");
  ACTS_VERBOSE("Prev. smoothed covariance:\n" << smoothedCovariance(prev_ts));

  // And the smoothed covariance
  smoothedCovariance(ts) =
      filteredCovariance(ts) -
      G * (predictedCovariance(prev_ts) - smoothedCovariance(prev_ts)) *
          G.transpose();

  // Check if the covariance matrix is semi-positive definite.
  // If not, make one (could do more) attempt to replace it with the
  // nearest semi-positive def matrix,
  // but it could still be non semi-positive
  BoundSymMatrix smoothedCov = smoothedCovariance(ts);
  if (not detail::covariance_helper<BoundSymMatrix>::validate(smoothedCov)) {
    ACTS_DEBUG(
        "Smoothed covariance is not positive definite. Could result in "
        "negative covariance!");
  }
  // Reset smoothed covariance
  smoothedCovariance(ts) = smoothedCov;
  ACTS_VERBOSE("Smoothed covariance is: \n" << smoothedCovariance(ts));

  return Result<void>::success();
}
}  // namespace Acts
