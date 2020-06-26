// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <boost/range/adaptors.hpp>
#include <memory>
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/detail/covariance_helper.hpp"
#include "Acts/Fitter/KalmanFitterError.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"

namespace Acts {

/// @brief Kalman smoother implementation based on Gain matrix formalism
///
/// @tparam parameters_t Type of the track parameters
/// @tparam jacobian_t Type of the Jacobian
class GainMatrixSmoother {
 public:
  /// @brief Gain Matrix smoother implementation
  ///

  /// Constructor with (non-owning) logger
  /// @param logger a logger instance
  GainMatrixSmoother(
      std::shared_ptr<const Logger> logger = std::shared_ptr<const Logger>(
          getDefaultLogger("GainMatrixSmoother", Logging::INFO).release()))
      : m_logger(std::move(logger)) {}

  /// Operater for Kalman smoothing
  ///
  /// @tparam source_link_t The type of source link
  ///
  /// @param gctx The geometry context for the smoothing
  /// @param trajectory The trajectory to be smoothed
  /// @param entryIndex The index of state to start the smoothing
  /// @param globalTrackParamsCovPtr The pointer to global track parameters
  /// covariance matrix
  ///
  /// @return The smoothed track parameters at the first measurement state
  template <typename source_link_t>
  Result<void> operator()(const GeometryContext& gctx,
                          MultiTrajectory<source_link_t>& trajectory,
                          size_t entryIndex) const {
    ACTS_VERBOSE("Invoked GainMatrixSmoother on entry index: " << entryIndex);
    using namespace boost::adaptors;

    using Matrix =
        ActsSymMatrixD<MultiTrajectory<source_link_t>::ParametersSize>;

    // For the last state: smoothed is filtered - also: switch to next
    ACTS_VERBOSE("Getting previous track state");
    auto prev_ts = trajectory.getTrackState(entryIndex);

    prev_ts.smoothed() = prev_ts.filtered();
    prev_ts.smoothedCovariance() = prev_ts.filteredCovariance();

    // Smoothing gain matrix
    Matrix G;

    // make sure there is more than one track state
    std::optional<std::error_code> error{std::nullopt};  // assume ok
    if (prev_ts.previous() == Acts::detail_lt::IndexData::kInvalid) {
      ACTS_VERBOSE("Only one track state given, smoothing terminates early");
    } else {
      ACTS_VERBOSE("Start smoothing from previous track state at index: "
                   << prev_ts.previous());

      trajectory.applyBackwards(prev_ts.previous(), [&prev_ts, &G, &error,
                                                     this](auto ts) {
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

        // Gain smoothing matrix
        // NB: The jacobian stored in a state is the jacobian from previous
        // state to this state in forward propagation
        G = ts.filteredCovariance() * prev_ts.jacobian().transpose() *
            prev_ts.predictedCovariance().inverse();

        if (G.hasNaN()) {
          error = KalmanFitterError::SmoothFailed;  // set to error
          return false;                             // abort execution
        }

        ACTS_VERBOSE("Gain smoothing matrix G:\n" << G);

        ACTS_VERBOSE("Calculate smoothed parameters:");
        ACTS_VERBOSE("Filtered parameters: " << ts.filtered().transpose());
        ACTS_VERBOSE(
            "Prev. smoothed parameters: " << prev_ts.smoothed().transpose());
        ACTS_VERBOSE(
            "Prev. predicted parameters: " << prev_ts.predicted().transpose());

        // Calculate the smoothed parameters
        ts.smoothed() =
            ts.filtered() + G * (prev_ts.smoothed() - prev_ts.predicted());

        ACTS_VERBOSE("Smoothed parameters are: " << ts.smoothed().transpose());

        ACTS_VERBOSE("Calculate smoothed covariance:");
        ACTS_VERBOSE("Prev. smoothed covariance:\n"
                     << prev_ts.smoothedCovariance());

        // And the smoothed covariance
        ts.smoothedCovariance() =
            ts.filteredCovariance() -
            G * (prev_ts.predictedCovariance() - prev_ts.smoothedCovariance()) *
                G.transpose();

        // Check if the covariance matrix is semi-positive definite.
        // If not, make one (could do more) attempt to replace it with the
        // nearest semi-positive def matrix,
        // but it could still be non semi-positive
        Matrix smoothedCov = ts.smoothedCovariance();
        if (not detail::covariance_helper<Matrix>::validate(smoothedCov)) {
          ACTS_DEBUG(
              "Smoothed covariance is not positive definite. Could result in "
              "negative covariance!");
        }
        // Reset smoothed covariance
        ts.smoothedCovariance() = smoothedCov;
        ACTS_VERBOSE("Smoothed covariance is: \n" << ts.smoothedCovariance());

        prev_ts = ts;
        return true;  // continue execution
      });
    }
    if (error) {
      // error is set, return result
      return *error;
    }

    // construct parameters from last track state
    return Result<void>::success();
  }

  /// Pointer to a logger that is owned by the parent, KalmanFilter
  std::shared_ptr<const Logger> m_logger{nullptr};

  /// Getter for the logger, to support logging macros
  const Logger& logger() const {
    assert(m_logger);
    return *m_logger;
  }
};
}  // namespace Acts