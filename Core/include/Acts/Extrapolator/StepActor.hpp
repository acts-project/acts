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

// TODO: This requires a test if there actually exists sth like a covariance matrix in the stepper
// TODO: Pure post step update
struct StepActor
{
  /// multiple scattering switch on/off
  bool multipleScattering = true;
  /// the scattering struct
  detail::HighlandScattering scattering;

  /// energy loss switch on/off
  bool energyLoss = true;
  /// the energy loss struct
  detail::IonisationLoss ionisationloss;
	

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
		if(state.stepping.pathAccumulated == 0.)
		 return;

	  if(state.navigation.currentVolume && state.navigation.currentVolume->material())
	  {
		  std::shared_ptr<const Material> matVol = state.navigation.currentVolume->material();
		  
			// to integrate process noise, we need to transport
			// the covariance to the current position in space
			if (state.stepping.covTransport) {
			  state.stepping.covarianceTransport(true);
			}
			
			const double thickness = state.stepping.stepSize; // TODO: Assuming there are no further corrections necessary
			state.stepping.stepSize = 1. * units::_mm; // TODO: hardcoded, needs to be released if volume is left
		  
		    // the momentum at current position
			const double p     = state.stepping.p;
			const double m     = state.options.mass;
			const double E     = std::sqrt(p * p + m * m);
			const double lbeta = p / E;
		  
		    // apply the multiple scattering
			// - only when you do covariance transport
			if (multipleScattering && state.stepping.covTransport) {
			  // thickness in X0 from without path correction
			  double tInX0 = thickness / matVol->X0();
			  // retrieve the scattering contribution
			  double sigmaScat = scattering(p, lbeta, tInX0);
			  double sinTheta
				  = std::sin(VectorHelpers::theta(state.stepping.direction()));
			  double sigmaDeltaPhiSq
				  = sigmaScat * sigmaScat / (sinTheta * sinTheta);
			  double sigmaDeltaThetaSq = sigmaScat * sigmaScat;
			  // good in any case for positive direction
			  if (state.stepping.navDir == forward) {
				// just add the multiple scattering component
				state.stepping.cov(ePHI, ePHI)
					+= state.stepping.navDir * sigmaDeltaPhiSq;
				state.stepping.cov(eTHETA, eTHETA)
					+= state.stepping.navDir * sigmaDeltaThetaSq;
			  } else {
				// we check if the covariance stays positive
				double sEphi   = state.stepping.cov(ePHI, ePHI);
				double sEtheta = state.stepping.cov(eTHETA, eTHETA);
				if (sEphi > sigmaDeltaPhiSq && sEtheta > sigmaDeltaThetaSq) {
				  // noise removal is not applied if covariance would fall below 0
				  state.stepping.cov(ePHI, ePHI) -= sigmaDeltaPhiSq;
				  state.stepping.cov(eTHETA, eTHETA) -= sigmaDeltaThetaSq;
				}
			  }
			}
			// apply the energy loss
			if (energyLoss) {
			  // calculate gamma
			  const double lgamma = E / m;
			  // energy loss and straggling - per unit length
			  std::pair<double, double> eLoss
				  = ionisationloss(m, lbeta, lgamma, *matVol, thickness);
			  // apply the energy loss
			  const double dEdl = state.stepping.navDir * eLoss.first;
			  const double dE   = thickness * dEdl;
			  // check for energy conservation, and only apply momentum change
			  // when kinematically allowed
			  if (E + dE > m) {
				// calcuate the new momentum
				const double newP = std::sqrt((E + dE) * (E + dE) - m * m);
				// update the state/momentum
				state.stepping.p = std::copysign(newP, state.stepping.p);
			  }
			  // transfer this into energy loss straggling and appply to covariance:
			  // do that even if you had not applied energy loss do to
			  // the kineamtic limit to catch the cases of deltE < MOP/MPV
			  if (state.stepping.covTransport) {
				// calculate the straggling
				const double sigmaQoverP
					= thickness * eLoss.second / (lbeta * p * p);
				// good in any case for positive direction
				if (state.stepping.navDir == forward) {
				  state.stepping.cov(eQOP, eQOP)
					  += state.stepping.navDir * sigmaQoverP * sigmaQoverP;
				} else {
				  // check that covariance entry doesn't become neagive
				  double sEqop = state.stepping.cov(eQOP, eQOP);
				  if (sEqop > sigmaQoverP * sigmaQoverP) {
					state.stepping.cov(eQOP, eQOP)
						+= state.stepping.navDir * sigmaQoverP * sigmaQoverP;
				  }
				}
			  }
			}
	  }
  }
};
}



