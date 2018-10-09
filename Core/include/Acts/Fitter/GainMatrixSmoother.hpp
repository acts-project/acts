// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <boost/variant.hpp>
#include <memory>
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/detail/trackstate_manipulation.hpp"

namespace Acts {

/// @brief Kalman smoother implementation based on Gain matrix formalism
///
///
/// @tparam parameters_t Type of the track parameters
/// @tparam jacobian_t Type of the Jacobian
template <typename parameters_t, typename jacobian_t>
class GainMatrixSmoother
{

public:
  /// @todo write documentation

  template <typename track_states_t>
  boost::optional<parameters_t>
  operator()(track_states_t& filteredStates) const
  {
    // The reverse iteration
    auto rit = filteredStates.rbegin();
    // For the last measurement the filtered state and the smoothed state are

    using ps_t = ParametricState<parameters_t, jacobian_t>;
    // smoothed parameter vector and covariance matrix
    typename parameters_t::ParVector_t smoothedPars;
    typename parameters_t::CovMatrix_t smoothedCov;

    // For the last state: smoothed is filtered - also: swicth to next
    auto& lState    = detail::parametricState<ps_t>(*rit++);
    lState.smoothed = lState.filtered.get();

    // Smoothing gain matrix
    using GMatrix = ActsMatrixD<Acts::NGlobalPars, Acts::NGlobalPars>;
    GMatrix G;

    // Loop and smooth ofer the remaining
    for (; rit != filteredStates.rend(); ++rit) {

      // The current state
      auto& cState = detail::parametricState<ps_t>(*rit);

      /// Gain smoothing matrix
      G = (*cState.filtered.get().covariance())
          * cState.jacobian.get().transpose()
          * (*lState.predicted.get().covariance()).inverse();
      // Calculate the smoothed parameters
      smoothedPars = cState.filtered.get().parameters()
          + G * (lState.smoothed.get().parameters()
                 - lState.predicted.get().parameters());

      // And the smoothed covariance
      smoothedCov = (*cState.filtered.get().covariance())
          - G * (*(lState.predicted.get().covariance())
                 - (*lState.smoothed.get().covariance()))
              * G.transpose();

      // Create smoothed track parameters
      cState.smoothed = parameters_t(
          std::make_unique<const decltype(smoothedCov)>(std::move(smoothedCov)),
          smoothedPars,
          cState.filtered.get().referenceSurface());

      // Point last state to current state
      lState = cState;
    }
    // The result is the pointer to the last smoothed state - for the cache
    return lState.smoothed;
  }
};
}
