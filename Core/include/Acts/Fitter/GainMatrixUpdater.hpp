// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <memory>
#include <variant>
#include "Acts/EventData/Measurement.hpp"
#include "Acts/EventData/MeasurementHelpers.hpp"
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Fitter/KalmanFitterError.hpp"
#include "Acts/Fitter/detail/VoidKalmanComponents.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"

namespace Acts {

/// @brief Update step of Kalman Filter using gain matrix formalism
///
/// @tparam parameters_t Type of the parameters to be updated
/// @tparam jacobian_t Type of the Transport jacobian
///
template <typename parameters_t>
class GainMatrixUpdater {
  using jacobian_t = typename parameters_t::CovMatrix_t;

 public:
  /// Explicit constructor
  ///
  /// @param calibrator is the calibration struct/class that converts
  /// uncalibrated measurements into calibrated ones
  /// @param logger a logger instance
  GainMatrixUpdater(
      calibrator_t calibrator = calibrator_t(),
      std::shared_ptr<const Logger> logger = std::shared_ptr<const Logger>(
          getDefaultLogger("GainMatrixUpdater", Logging::INFO).release()))
      : m_logger(std::move(logger)), m_mCalibrator(std::move(calibrator)) {}

  /// @brief Public call operator for the boost visitor pattern
  ///
  /// @tparam track_state_t Type of the track state for the update
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param trackState the measured track state
  ///
  /// @return Bool indicating whether this update was 'successful'
  /// @note Non-'successful' updates could be holes or outliers,
  ///       which need to be treated differently in calling code.
  template <typename track_state_t>
  Result<void> operator()(const GeometryContext& /*gctx*/,
                          track_state_t trackState) const {
    ACTS_VERBOSE("Invoked GainMatrixUpdater");
    // let's make sure the types are consistent
    using SourceLink = typename track_state_t::SourceLink;
    using TrackStateProxy =
        typename MultiTrajectory<SourceLink>::TrackStateProxy;
    static_assert(std::is_same_v<track_state_t, TrackStateProxy>,
                  "Given track state type is not a track state proxy");

    // we should definitely have an uncalibrated measurement here
    assert(trackState.hasUncalibrated());
    // there should be a calibrated measurement
    assert(trackState.hasCalibrated());
    // we should have predicted state set
    assert(trackState.hasPredicted());
    // filtering should not have happened yet, but is allocated, therefore set
    assert(trackState.hasFiltered());

    // read-only handles. Types are eigen maps to backing storage
    const auto predicted = trackState.predicted();
    // const auto predicted = *trackState.parameter.predicted;
    const auto predicted_covariance = trackState.predictedCovariance();
    // const CovMatrix_t& predicted_covariance = *predicted.covariance();

    // ParVector_t filtered_parameters;
    // CovMatrix_t filtered_covariance;
    // read-write handles. Types are eigen maps into backing storage.
    // This writes directly into the trajectory storage
    auto filtered = trackState.filtered();
    auto filtered_covariance = trackState.filteredCovariance();

    visit_measurement(
        trackState.calibrated(), trackState.calibratedCovariance(),
        trackState.calibratedSize(),
        [&](const auto calibrated, const auto calibrated_covariance) {
          constexpr size_t measdim = decltype(calibrated)::RowsAtCompileTime;
          using cov_t = ActsSymMatrixD<measdim>;
          using par_t = ActsVectorD<measdim>;

          const ActsMatrixD<measdim, BoundParsDim> H =
              trackState.projector()
                  .template topLeftCorner<measdim, BoundParsDim>();

          const ActsMatrixD<BoundParsDim, measdim> K =
              predicted_covariance * H.transpose() *
              (H * predicted_covariance * H.transpose() + calibrated_covariance)
                  .inverse();

          filtered = predicted + K * (calibrated - H * predicted);
          filtered_covariance =
              (ActsSymMatrixD<
                   MultiTrajectory<SourceLink>::ParametersSize>::Identity() -
               K * H) *
              predicted_covariance;

          // calculate filtered residual
          par_t residual(trackState.calibratedSize());
          residual = (calibrated - H * filtered);

          trackState.chi2() =
              (residual.transpose() *
               ((cov_t::Identity() - H * K) * calibrated_covariance).inverse() *
               residual)
                  .value();
        });

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
