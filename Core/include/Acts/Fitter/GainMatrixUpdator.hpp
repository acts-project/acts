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
template <typename parameters_t,
          typename jacobian_t,
          typename calibrator_t = VoidKalmanComponents>
class GainMatrixUpdator
{

public:
  // Shorthand
  using predicted_state_t = std::tuple<parameters_t, jacobian_t, double>;

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
  /// @tparam predicted_state_t Type of the track state prediction
  ///
  /// @param m the measured track state
  /// @param predicted the predicted track state
  ///
  /// @return The optional parameters - indicating if the update happened
  template <typename track_state_t>
  boost::optional<parameters_t>
  operator()(track_state_t& m, predicted_state_t predicted) const
  {
    GainMatrixUpdatorImpl impl(m_mCalibrator, std::move(predicted));
    return boost::apply_visitor(impl, m);
  }

private:
  /// The measurement calibrator
  calibrator_t m_mCalibrator;

  /// @brief GainMatrix updator implementation
  struct GainMatrixUpdatorImpl
      : public boost::static_visitor<boost::optional<parameters_t>>
  {
    /// @brief Explicit constructor of the GainMatrix updator
    ///
    /// @param calibrator The calibration struct/class that converts
    /// uncalibrated measurements into calibrated ones
    /// @param predictedState The tuple of predicted parameters, jacobian, path
    explicit GainMatrixUpdatorImpl(const calibrator_t& calibrator,
                                   predicted_state_t   predictedState)
      : m_mCalibrator(&calibrator), m_pState(std::move(predictedState))
    {
    }

    /// @brief Call operator for the visitor pattern
    ///
    /// @tparam measurement_t Type of the measurement to be used
    /// @param rawMeasurement The measurement
    ///
    /// @todo Include the incremental chi-2 here or do some outlier steering ?
    ///
    /// @return The filtered parameters
    template <typename track_state_t>
    boost::optional<parameters_t>
    operator()(track_state_t& trackState) const
    {
      // Covariance matrix initialization
      static const ActsSymMatrixD<Acts::NGlobalPars> unit
          = ActsSymMatrixD<Acts::NGlobalPars>::Identity();

      // Predicted Parameters
      const auto& predicted = std::get<parameters_t>(m_pState);

      // Calibrate the measurement
      auto cMeasurement = (*m_mCalibrator)(
          trackState.measurement.uncalibrated.get(), predicted);
      const auto* pCov_trk = predicted.covariance();

      // Take the projector (measurement mapping function)
      const auto& H = cMeasurement.projector();

      // The Kalman gain matrix
      ActsMatrixD<Acts::NGlobalPars, track_state_t::size()> K = (*pCov_trk)
          * H.transpose()
          * (H * (*pCov_trk) * H.transpose() + cMeasurement.covariance())
                .inverse();

      // New parameters after update
      typename parameters_t::ParVector_t newParValues
          = predicted.parameters() + K * cMeasurement.residual(predicted);

      // New covaraincd after update
      typename parameters_t::CovMatrix_t newCov = (unit - K * H) * (*pCov_trk);

      // Create a new filtered parameters and covariance
      parameters_t filtered(
          std::make_unique<const typename parameters_t::CovMatrix_t>(
              std::move(newCov)),
          newParValues,
          predicted.referenceSurface().getSharedPtr());

      // Set (and move) everything
      trackState.measurement.calibrated = std::move(cMeasurement);
      trackState.parametric.predicted   = std::move(predicted);
      trackState.parametric.filtered    = std::move(filtered);
      trackState.parametric.jacobian
          = std::move(std::get<jacobian_t>(m_pState));
      trackState.parametric.pathLength = std::get<double>(m_pState);

      // Return the optional filtered state
      return trackState.parametric.filtered;
    }

  private:
    /// The Calibator
    const calibrator_t* m_mCalibrator;

    /// Predicted state (parameters, jacobian, path)
    predicted_state_t m_pState;
  };
};

}  // namespace Acts