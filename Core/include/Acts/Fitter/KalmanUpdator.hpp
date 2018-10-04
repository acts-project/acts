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
#include "Acts/Utilities/Definitions.hpp"

namespace Acts {

/// @brief Update step of Kalman Filter using gain matrix formalism
///
///
/// This is implemented as a boost vistor pattern for use of the
/// boost variant container
template <typename parameters_t,
          typename predicted_state_t,
          typename calibrator_t>
class GainMatrixUpdator
{

public:
  /// Explicit constructor
  ///
  GainMatrixUpdator(calibrator_t calibrator)
    : m_mCalibrator(std::move(calibrator))
  {
  }

  /// @brief Public call operator for the boost visitor pattern
  ///
  ///
  /// @return The optional parameters - indicating if the update happened
  template <typename track_state_t>
  const parameters_t*
  operator()(track_state_t& m, predicted_state_t predicted) const
  {
    GainMatrixUpdatorImpl impl(m_mCalibrator, predicted);
    return boost::apply_visitor(impl, m);
  }

private:
  /// The calibrator
  calibrator_t m_mCalibrator;

  struct GainMatrixUpdatorImpl
      : public boost::static_visitor<const parameters_t*>
  {
  public:
    /// @brief Explicit constructor of the GainMatrix updator
    ///
    /// @param pars The predicted parameters
    explicit GainMatrixUpdatorImpl(const calibrator_t& calibrator,
                                   predicted_state_t   predictedState)
      : m_mCalibrator(&calibrator), m_pState(std::move(predictedState))
    {
    }

    /// @brief Private call operator for the visitor pattern
    ///
    /// @tparam measurement_t Type of the measurement to be used
    /// @param rawMeasurement The measurement
    ///
    /// @todo Include the incremental chi-2 here or do some outlier steering ?
    ///
    /// @return The updated parameters
    template <typename track_state_t>
    const parameters_t*
    operator()(track_state_t& trackState) const
    {
      // Covariance matrix initialization
      static const ActsSymMatrixD<Acts::NGlobalPars> unit
          = ActsSymMatrixD<Acts::NGlobalPars>::Identity();

      // Predicated Parameters
      const auto& predicted = m_pState.parameters;
      const auto& jacobian  = m_pState.jacobian;

      // Calibrate the measurement
      auto cMeasurement = m_calibrator(trackState.measurement.get(), predicted);
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

      // Create a new updated parameters and covariance
      parameters_t updated(
          std::make_unique<const typename parameters_t::CovMatrix_t>(
              std::move(newCov)),
          newParValues,
          predicted.referenceSurface());

      // Prepare the return parameters pointers
      const parameters_t* updatedPtr = &updated;

      // Set & move everything
      trackState.calibratedMeasurement = std::move(cMeasurement);
      trackState.predicted             = std::move(m_pState.parameters);
      trackState.updated               = std::move(updated);
      trackState.jacobian              = std::move(m_pState.jacobian);
      trackState.pathLength            = m_pState.pathLength;

      // Return the pointer for stepping update
      return updatedPtr;
    }

  private:
    /// Predicted state (parameters, jacobian, path)
    predicted_state_t m_pState;

    /// The Calibator
    const calibrator_t* m_mCalibrator;
  };
};

}  // namespace Acts