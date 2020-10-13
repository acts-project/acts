// This file is part of the Acts project.
//
// Copyright (C) 2018-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/SourceLinkConcept.hpp"
#include "Acts/Utilities/Result.hpp"
#include "Acts/Utilities/TypeTraits.hpp"

namespace Acts {

/// @brief void Measurement calibrator and converter
struct VoidKalmanComponents {
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
  Result<measurement_t> operator()(measurement_t measurement,
                                   const parameters_t& /* parameters */) const {
    return measurement;
  }
};

/// @brief void Kalman updater
struct VoidKalmanUpdater {
  /// @brief Public call mimicking an updater
  ///
  /// @tparam measurement_t Type of the measurement to be used
  /// @tpredicted_state_t Type of the (bound) predicted state
  ///
  /// @param m The measurement
  /// @param predicted The predicted parameters
  ///
  /// @return The copied predicted parameters
  template <typename track_state_t, typename predicted_state_t>
  auto operator()(track_state_t& /* trackState */,
                  const predicted_state_t& predicted) const {
    return &(predicted.parameters);
  }
};

/// @brief void Kalman smoother
struct VoidKalmanSmoother {
  /// @brief Public call mimicking an updater
  ///
  /// @tparam track_states_t Type of the track states
  ///
  /// @param states The track states to be smoothed
  ///
  /// @return The resulting
  template <typename parameters_t, typename track_states_t>
  const parameters_t* operator()(track_states_t& /* trackStates */) const {
    return nullptr;
  }
};

/// @brief void outlier finder
struct VoidOutlierFinder {
  /// @brief Public call mimicking an outlier finder
  ///
  /// @tparam track_state_t Type of the track state
  ///
  /// @param trackState The trackState to investigate
  ///
  /// @return Whether it's outlier or not
  template <typename track_state_t>
  constexpr bool operator()(const track_state_t& /* trackState */) const {
    return false;
  }
};

}  // namespace Acts
