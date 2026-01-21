// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/TrackFitting/GainMatrixSmoother.hpp"

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/TrackParameterHelpers.hpp"
#include "Acts/EventData/detail/CovarianceHelper.hpp"
#include "Acts/TrackFitting/KalmanFitterError.hpp"

namespace Acts {

Result<void> GainMatrixSmoother::calculate(AnyMutableTrackStateProxy ts,
                                           AnyConstTrackStateProxy prev_ts,
                                           const Logger& logger) const {
  auto prevPredictedCovariance = prev_ts.predictedCovariance();
  auto filteredCovariance = ts.filteredCovariance();
  auto prevJacobian = prev_ts.jacobian();

  ACTS_VERBOSE("Prev. predicted covariance\n"
               << prevPredictedCovariance << "\n, inverse: \n"
               << prevPredictedCovariance.inverse());

  // Gain smoothing matrix
  // NB: The jacobian stored in a state is the jacobian from previous
  // state to this state in forward propagation
  BoundMatrix G = filteredCovariance * prevJacobian.transpose() *
                  prevPredictedCovariance.inverse();

  if (G.hasNaN()) {
    ACTS_VERBOSE("Gain smoothing matrix G has NaNs");

    ACTS_VERBOSE("Filtered covariance:\n" << filteredCovariance);
    ACTS_VERBOSE("Jacobian:\n" << prevJacobian);
    ACTS_VERBOSE("Predicted covariance:\n" << prevPredictedCovariance);
    ACTS_VERBOSE("Inverse of predicted covariance:\n"
                 << prevPredictedCovariance.inverse());

    ACTS_VERBOSE("Gain smoothing matrix G:\n" << G);

    return KalmanFitterError::SmoothFailed;
  }

  ACTS_VERBOSE("Gain smoothing matrix G:\n" << G);

  auto prevPredicted = prev_ts.predicted();
  auto prevSmoothed = prev_ts.smoothed();

  auto filtered = ts.filtered();
  auto smoothed = ts.smoothed();
  auto smoothedCovariance = ts.smoothedCovariance();
  auto prevSmoothedCovariance = prev_ts.smoothedCovariance();

  ACTS_VERBOSE("Calculate smoothed parameters:");
  ACTS_VERBOSE("Filtered parameters: " << filtered.transpose());
  ACTS_VERBOSE("Prev. smoothed parameters: " << prevSmoothed.transpose());
  ACTS_VERBOSE("Prev. predicted parameters: " << prevPredicted.transpose());

  // Calculate the smoothed parameters
  smoothed =
      filtered + G * subtractBoundParameters(prevSmoothed, prevPredicted);
  // Normalize phi and theta
  smoothed = normalizeBoundParameters(smoothed);

  ACTS_VERBOSE("Smoothed parameters are: " << smoothed.transpose());
  ACTS_VERBOSE("Calculate smoothed covariance:");
  ACTS_VERBOSE("Prev. smoothed covariance:\n" << prevSmoothedCovariance);

  // And the smoothed covariance
  smoothedCovariance =
      filteredCovariance +
      G * (prevSmoothedCovariance - prevPredictedCovariance) * G.transpose();

  if (doCovCheckAndAttemptFix) {
    // Check if the covariance matrix is semi-positive definite.
    // If not, make one (could do more) attempt to replace it with the
    // nearest semi-positive def matrix,
    // but it could still be non semi-positive
    BoundSquareMatrix smoothedCov = smoothedCovariance;
    if (!detail::CovarianceHelper<BoundSquareMatrix>::validate(smoothedCov)) {
      ACTS_DEBUG(
          "Smoothed covariance is not positive definite. Could result in "
          "negative covariance!");
    }
    // Reset smoothed covariance
    smoothedCovariance = smoothedCov;
  }

  ACTS_VERBOSE("Smoothed covariance is: \n" << smoothedCovariance);

  return Result<void>::success();
}

}  // namespace Acts
