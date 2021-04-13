// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Definitions/Algebra.hpp"

namespace Acts {

struct DenseScatteringInteractor {
  struct Result {
	  ActsScalar s = 0.;
	  ActsScalar qop0 = 0.;
	  ActsScalar t0 = 0.;
	  TrackingVolume const* currentVolume = nullptr;
	  ActsScalar sAccumulated = 0.;
	  ActsScalar sInX0 = 0.;
	  Vector3 dir0 = Vector3::Zero();
	  Vector3 previousPos = Vector3::Zero();
  };
  using result_type = Result;

  template <typename propagator_state_t, typename stepper_t>
  void operator()(propagator_state_t& state, const stepper_t& stepper,
                  result_type& result) const {
	
	// The volume changed
	if(result.currentVolume != state.navigation.currentVolume)
	{
		if(result.s != 0.)
		{
			BoundSymMatrix msCovariance = evaluateMultipleScattering(state.options.absPdgCode, state.options.mass, stepper.charge(state.stepping), 
		stepper.momentum(state.stepping), stepper.time(state.stepping), result);
			state.stepping.cov += msCovariance;
		}
		
		// Set the current volume
		result.currentVolume = state.navigation.currentVolume;
		
		// Reset the tracker variables
		/// Should occur whenever the volume changes or a transport occurs
		reset(stepper.charge(state.stepping), stepper.momentum(state.stepping), stepper.time(state.stepping), stepper.direction(state.stepping), result);
	} else {
	
	// Same volume with material: accumulate path
	if(result.currentVolume != nullptr && result.currentVolume->volumeMaterial() != nullptr)
	{
		ActsScalar deltaS = state.stepping.pathAccumulated - result.sAccumulated;
		result.s += deltaS;
		result.sInX0 += deltaS / result.currentVolume->volumeMaterial()->material(result.previousPos).X0();
	}
	
	// If the covariance was transported then the MS part is added
	if(state.stepping.jacTransport == FreeMatrix::Identity() || state.navigation.targetReached)
	{
		BoundSymMatrix msCovariance = evaluateMultipleScattering(state.options.absPdgCode, state.options.mass, stepper.charge(state.stepping), 
		stepper.momentum(state.stepping), stepper.time(state.stepping), result);
			state.stepping.cov += msCovariance;
		
		reset(stepper.charge(state.stepping), stepper.momentum(state.stepping), stepper.time(state.stepping), stepper.direction(state.stepping), result);
	}
}
	result.sAccumulated = state.stepping.pathAccumulated;
	result.previousPos = stepper.position(state.stepping);
  }

 private:

  void 
  reset(ActsScalar q, ActsScalar momentum, ActsScalar time, Vector3 direction, result_type& result) const;
  
  BoundSymMatrix evaluateMultipleScattering(int pdg, float mass, float q, ActsScalar momentum, ActsScalar time,
                    result_type& result) const;
};  // namespace Acts
}  // namespace Acts
