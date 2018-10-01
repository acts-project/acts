// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

namespace Acts {

/// @brief void measurement calibrator
struct VoidKalmanComponents
{
  /// @brief Public call mimicking a calibrator
  ///
  /// @tparam measurement_t Type of the measurement
  /// @tparam parameter_t Type of the parameters for calibration
  ///
  /// @param m Measurement to be moved through
  /// @param pars Parameters to be used for calibration
  ///
  /// @return void-calibrated measurement
  template <typename measurement_t, typename parameters_t>
  measurement_t
  operator()(measurement_t m, const parameters_t& /*pars*/) const
  {
    return std::move(m);
  }

  /// @brief Public call mimicking an updator
  ///
  /// @tparam measurement_t Type of the measurement to be used
  /// @param m The measurement
  /// @param pars The predicted parameters
  ///
  /// @return The copied predicted parameters
  template <typename measurement_t, typename parameters_t>
  parameters_t
  operator()(const measurement_t& /*m*/, const parameters_t& pars) const
  {
    return pars;
  }

  /// @brief void measurement converter only moves the
  /// the measurement through for further processing
  ///
  /// @tparam measurement_container_t Type of the measurement
  ///
  /// @param ms Measurements to be moved through
  ///
  /// @return moved measurements
  template <typename measurements_t>
  measurements_t
  operator()(measurements_t ms) const
  {
    return std::move(ms);
  }
};

}  // namespace Acts