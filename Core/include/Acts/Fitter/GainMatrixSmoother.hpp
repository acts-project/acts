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

    using ps_t        = ParametricState<parameters_t, jacobian_t>;
    using ParVector_t = typename parameters_t::ParVector_t;
    using CovMatrix_t = typename parameters_t::CovMatrix_t;
    // smoothed parameter vector and covariance matrix
    ParVector_t smoothedPars;
    CovMatrix_t smoothedCov;

    // For the last state: smoothed is filtered - also: switch to next
    auto& lState    = detail::getParametricState<ps_t>(*rit);
    lState.smoothed = lState.filtered.get();

    // Smoothing gain matrix
    using GMatrix = ActsMatrixD<Acts::NGlobalPars, Acts::NGlobalPars>;
    GMatrix G;

    // Loop and smooth the remaining states
    for (++rit; rit != filteredStates.rend(); ++rit) {

      // The current state
      auto& cState = detail::getParametricState<ps_t>(*rit);

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
      parameters_t smoothed(
          std::make_unique<CovMatrix_t>(std::move(smoothedCov)),
          smoothedPars,
          cState.filtered.get().referenceSurface().getSharedPtr());

      cState.smoothed = std::move(smoothed);

      // Point last state to current state
      lState = cState;
    }
    // The result is the pointer to the last smoothed state - for the cache
    return lState.smoothed;
  }
};
}
