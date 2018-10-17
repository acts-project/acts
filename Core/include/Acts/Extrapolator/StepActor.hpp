// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <cmath>
#include "Acts/Extrapolator/detail/InteractionFormulas.hpp"
#include "Acts/Utilities/Helpers.hpp"

namespace Acts {

using cstep = detail::ConstrainedStep;

// TODO: Logging
struct StepActor
{
  template <typename propagator_state_t>
  void
  operator()(propagator_state_t& state) const
  {
    // if we are on target, everything should have been done
    if (state.navigation.targetReached) {
      return;
    }
    // if switched off, then return - alows run-time configuration
    if (!multipleScattering && !energyLoss) {
      return;
    }
    // No action at first step
    if (state.stepping.pathAccumulated == 0.) {
      if (state.navigation.currentVolume
          && state.navigation.currentVolume->material()
          && state.stepping.stepSize > maxStepSize) {
        state.stepping.stepSize.update(maxStepSize, cstep::user);
      }
      return;
    }

    if (state.navigation.currentVolume
        && state.navigation.currentVolume->material()) {
      std::shared_ptr<const Material> matVol
          = state.navigation.currentVolume->material();

      // to integrate process noise, we need to transport
      // the covariance to the current position in space
      if (state.stepping.covTransport) {
        state.stepping.covarianceTransport(false);
      }


      }

      // apply the energy loss
      if (energyLoss) {
        if (state.stepping.stepSize.value(cstep::user) > maxStepSize) {
          state.stepping.stepSize.update(maxStepSize, cstep::user);
        }

      }
		else 
		{
      if (energyLoss
          && state.options.maxStepSize
              > state.stepping.stepSize.value(cstep::user)) {
        state.stepping.stepSize.release(cstep::user);
      }
    }
  }
};
}
