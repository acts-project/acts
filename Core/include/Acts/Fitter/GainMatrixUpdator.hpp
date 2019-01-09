// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <boost/variant.hpp>
#include <memory>
#include "Acts/EventData/Measurement.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Fitter/detail/VoidKalmanComponents.hpp"
#include "Acts/Utilities/Definitions.hpp"

namespace Acts {

/// @brief Update step of Kalman Filter using gain matrix formalism
///
/// @tparam parameters_t Type of the parameters to be updated
/// @tparam jacobian_t Type of the Transport jacobian
/// @tparam calibrator_t A measurement calibrator (can be void)
///
/// This is implemented as a boost vistor pattern for use of the
/// boost variant container
template <typename parameters_t, typename calibrator_t = VoidKalmanComponents>
class GainMatrixUpdator
{

  using jacobian_t = typename parameters_t::CovMatrix_t;

public:
  /// Explicit constructor
  ///
  /// @param calibrator is the calibration struct/class that converts
  /// uncalibrated measurements into calibrated ones
  GainMatrixUpdator(calibrator_t calibrator = calibrator_t())
    : m_mCalibrator(std::move(calibrator))
  {
  }

  /// @brief Public call operator for the boost visitor pattern
  ///
  /// @tparam track_state_t Type of the track state for the update
  ///
  /// @param trackState the measured track state
  ///
  /// @return Bool indicating whether this update was 'successful'
  /// @note Non-'successful' updates could be holes or outliers,
  ///       which need to be treated differently in calling code.
  template <typename track_state_t>
  bool
  operator()(track_state_t& trackState) const
  {
    using CovMatrix_t = typename parameters_t::CovMatrix_t;
    using ParVector_t = typename parameters_t::ParVector_t;

    // we should definitely have an uncalibrated measurement here
    assert(trackState.measurement.uncalibrated);
    // there should be no calibrated measurement
    assert(!trackState.measurement.calibrated);
    // we should have predicted state set
    assert(trackState.parameter.predicted);
    // filtring should not have happened yet
    assert(!trackState.parameter.filtered);

    // read-only prediction handle
    const parameters_t& predicted = *trackState.parameter.predicted;

    const CovMatrix_t& predicted_covariance = *predicted.covariance();

    ParVector_t filtered_parameters;
    CovMatrix_t filtered_covariance;

    // we need to remove type-erasure on the measurement type
    // to access its methods
    boost::apply_visitor(
        [&](const auto& uncalibrated) {
          // type of measurement
          using meas_t = typename std::remove_const<
              typename std::remove_reference<decltype(uncalibrated)>::type>::
              type;
          // type of projection
          using projection_t = typename meas_t::Projection_t;
          // type of gain matrix (transposed projection)
          using gain_matrix_t = ActsMatrixD<projection_t::ColsAtCompileTime,
                                            projection_t::RowsAtCompileTime>;

          // Calibrate the measurement
          meas_t calibrated = m_mCalibrator(uncalibrated, predicted);

          // Take the projector (measurement mapping function)
          const projection_t& H = calibrated.projector();

          // The Kalman gain matrix
          gain_matrix_t K = predicted_covariance * H.transpose()
              * (H * predicted_covariance * H.transpose()
                 + calibrated.covariance())
                    .inverse();

          // filtered new parameters after update
          filtered_parameters
              = predicted.parameters() + K * calibrated.residual(predicted);

          // updated covariance after filtering
          filtered_covariance
              = (CovMatrix_t::Identity() - K * H) * predicted_covariance;

          // plug calibrated measurement back into track state
          trackState.measurement.calibrated = std::move(calibrated);

        },
        *trackState.measurement.uncalibrated);

    // Create new filtered parameters and covariance
    parameters_t filtered(
        std::make_unique<const CovMatrix_t>(std::move(filtered_covariance)),
        filtered_parameters,
        predicted.referenceSurface().getSharedPtr());

    trackState.parameter.filtered = std::move(filtered);

    // always succeed, no outlier logic yet
    return true;
  }

private:
  /// The measurement calibrator
  calibrator_t m_mCalibrator;
};

}  // namespace Acts