// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

namespace Acts {

// TODO: Logging
struct StepActor
{
	// Configurations for 
    /// Boolean flag for energy loss while stepping
    bool energyLossFlag = true;
    /// Tolerance for the error of the integration
	double tolerance = 5e-5;
    /// Boolean flag for inclusion of d(dEds)d(q/p) into energy loss
	bool includeGgradient = false;
	/// Cut-off value for the momentum
    double momentumCutOff = 0.;
    
  template <typename propagator_state_t>
  void
  operator()(propagator_state_t& state) const
  {
    // if we are on target, everything should have been done
    if (state.navigation.targetReached) {
      return;
    }

    if (state.navigation.currentVolume
        && state.navigation.currentVolume != state.stepping.volume) 
	{
		state.stepping.mass = state.options.mass;
		state.stepping.volume = state.navigation.currentVolume;
		state.stepping.energyLossFlag = energyLossFlag;
		state.stepping.tolerance = tolerance;
		state.stepping.includeGgradient = includeGgradient;
		state.stepping.momentumCutOff = momentumCutOff;
	}
  }
};
}
