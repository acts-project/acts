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
#include "Acts/Utilities/ParameterDefinitions.hpp"
#include "Acts/Utilities/Result.hpp"

namespace Acts {

/// @brief Kalman smoother implementation based on Gain matrix formalism
///
/// @tparam parameters_t Type of the track parameters
/// @tparam jacobian_t Type of the Jacobian
template <typename parameters_t>
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

  template <typename source_link_t>
  Result<parameters_t> operator()(const GeometryContext& gctx,
                                  MultiTrajectory<source_link_t>& trajectory,
                                  size_t entryIndex,
                                  GlobalBoundSymMatrix& optGlobalCov) const {
    ACTS_VERBOSE("Invoked GainMatrixSmoother on entry index: " << entryIndex);
    using namespace boost::adaptors;

    // using ParVector_t = typename parameters_t::ParVector_t;
    using CovMatrix_t = typename parameters_t::CovMatrix_t;
    using gain_matrix_t = CovMatrix_t;

    // The bound parameters size
    constexpr size_t parametersSize = eBoundParametersSize;

    // For the last state: smoothed is filtered - also: switch to next
    ACTS_VERBOSE("Getting previous track state");
    auto prev_ts = trajectory.getTrackState(entryIndex);

    prev_ts.smoothed() = prev_ts.filtered();
    prev_ts.smoothedCovariance() = prev_ts.filteredCovariance();

    // Fill the covariance of last state
    if (optGlobalCov.size() != 0) {
      size_t globalCovSize = optGlobalCov.rows();
      assert(globalCovSize % parametersSize == 0);
      ACTS_VERBOSE("Size of global covariance matrix is: " << globalCovSize);
      optGlobalCov.block<parametersSize, parametersSize>(
          globalCovSize - parametersSize, globalCovSize - parametersSize) =
          prev_ts.smoothedCovariance();
    }

    // Smoothing gain matrix
    gain_matrix_t G;

    // make sure there is more than one track state
    std::optional<std::error_code> error{std::nullopt};  // assume ok
    if (prev_ts.previous() == Acts::detail_lt::IndexData::kInvalid) {
      ACTS_VERBOSE("Only one track state given, smoothing terminates early");
    } else {
      ACTS_VERBOSE("Start smoothing from previous track state at index: "
                   << prev_ts.previous());

      size_t nSmoothed = 1;
      trajectory.applyBackwards(prev_ts.previous(), [&prev_ts, &G, &error,
                                                     &nSmoothed,
                                                     &parametersSize,
                                                     &optGlobalCov,
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
        CovMatrix_t smoothedCov = ts.smoothedCovariance();
        if (not detail::covariance_helper<CovMatrix_t>::validate(smoothedCov)) {
          ACTS_DEBUG(
              "Smoothed covariance is not positive definite. Could result in "
              "negative covariance!");
        }
        // Reset smoothed covariance
        ts.smoothedCovariance() = smoothedCov;
        ACTS_VERBOSE("Smoothed covariance is: \n" << ts.smoothedCovariance());

        if (optGlobalCov.size() != 0) {
          // Fill global track parameters covariance matrix
          size_t globalCovSize = optGlobalCov.rows();
          size_t nStates = globalCovSize / parametersSize;
          // Fill the diagonal element
          size_t iRow = globalCovSize - parametersSize * (nSmoothed + 1);
          optGlobalCov.block<parametersSize, parametersSize>(iRow, iRow) =
              ts.smoothedCovariance();
          // Fill the correlation between this state and already smoothed states
          for (size_t iSmoothed = 1; iSmoothed <= nSmoothed; iSmoothed++) {
            size_t iCol = iRow + parametersSize * iSmoothed;
            CovMatrix_t prev_correlation =
                optGlobalCov.block<parametersSize, parametersSize>(
                    iRow + parametersSize, iCol);
            CovMatrix_t correlation = G * prev_correlation;
            ACTS_VERBOSE("Fill block of size ("
                         << parametersSize << ", " << parametersSize
                         << "), starting at (" << iRow << ", " << iCol
                         << ") for track parameters correlation:\n"
                         << correlation);
            optGlobalCov.block<parametersSize, parametersSize>(iRow, iCol) =
                correlation;
            ACTS_VERBOSE("Fill block of size ("
                         << parametersSize << ", " << parametersSize
                         << "), starting at (" << iCol << ", " << iRow
                         << ") for track parameters correlation:\n"
                         << correlation.transpose());
            optGlobalCov.block<parametersSize, parametersSize>(iCol, iRow) =
                correlation.transpose();
          }
        }

        prev_ts = ts;
        nSmoothed++;
        return true;  // continue execution
      });
    }
    if (error) {
      // error is set, return result
      return *error;
    }

    // construct parameters from last track state
    return prev_ts.smoothedParameters(gctx);
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
