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
#include "Acts/EventData/TrackParameters.hpp"

namespace Acts {

/// @brief Kalman smoother implementation based on Gain matrix formalism
///
/// @tparam parameters_t Type of the track parameters
/// @tparam jacobian_t Type of the Jacobian
template <typename parameters_t>
class GainMatrixSmoother {
  using jacobian_t = typename parameters_t::CovMatrix_t;

 public:
  /// @brief Gain Matrix smoother implementation
  ///

  template <typename track_states_t>
  boost::optional<parameters_t> operator()(
      const GeometryContext& gctx, track_states_t& filteredStates) const {
    using namespace boost::adaptors;

    using track_state_t = typename track_states_t::value_type;
    using ParVector_t = typename parameters_t::ParVector_t;
    using CovMatrix_t = typename parameters_t::CovMatrix_t;
    using gain_matrix_t = CovMatrix_t;

    // smoothed parameter vector and covariance matrix
    ParVector_t smoothedPars;
    CovMatrix_t smoothedCov;

    // For the last state: smoothed is filtered - also: switch to next
    auto* prev_ts = &filteredStates.back();
    assert(prev_ts->parameter.filtered);
    prev_ts->parameter.smoothed = *prev_ts->parameter.filtered;

    // Smoothing gain matrix
    gain_matrix_t G;

    // Loop and smooth the remaining states
    for (track_state_t& ts :
         filteredStates | sliced(0, filteredStates.size() - 1) | reversed) {
      // The current state
      assert(ts.parameter.filtered);
      assert(ts.parameter.predicted);
      assert(ts.parameter.jacobian);
      assert(ts.parameter.predicted->covariance() != nullptr);
      assert(ts.parameter.filtered->covariance() != nullptr);

      assert(prev_ts->parameter.smoothed);
      assert(prev_ts->parameter.predicted);

      // clang-format off

      // Gain smoothing matrix
      G = (*ts.parameter.filtered->covariance())
          * ts.parameter.jacobian->transpose()
          * (*prev_ts->parameter.predicted->covariance()).inverse();
      // Calculate the smoothed parameters
      smoothedPars = ts.parameter.filtered->parameters()
                     + G * (prev_ts->parameter.smoothed->parameters()
                            - prev_ts->parameter.predicted->parameters());

      // And the smoothed covariance
      smoothedCov = (*ts.parameter.filtered->covariance())
                    - G * (*(prev_ts->parameter.predicted->covariance())
                           - (*prev_ts->parameter.smoothed->covariance()))
                           * G.transpose();

      // clang-format on

      // Create smoothed track parameters
      ts.parameter.smoothed = parameters_t(
          gctx, std::make_unique<CovMatrix_t>(std::move(smoothedCov)),
          smoothedPars, ts.referenceSurface().getSharedPtr());

      // Point prev state to current state
      prev_ts = &ts;
    }
    // The result is the pointer to the last smoothed state - for the cache
    return prev_ts->parameter.smoothed;
  }
};
}  // namespace Acts
