// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <cmath>
#include <sstream>
#include <utility>

#include "Acts/Material/MaterialProperties.hpp"
#include "Acts/Material/Interactions.hpp"
#include "Acts/Surfaces/Surface.hpp"

namespace Acts {
namespace detail {
	struct Data
	{
		template<typename propagator_state_t, typename stepper_t>
		Data(const Surface* sSurface, const propagator_state_t& state, const stepper_t& stepper) 
			: surface(sSurface), pos(stepper.position(state.stepping)), time(stepper.time(state.stepping)), dir(stepper.direction(state.stepping)), 
			momentum(stepper.momentum(state.stepping)), q(stepper.charge(state.stepping)), qOverP(q / momentum), mass(state.options.mass), pdg(state.options.absPdgCode), 
			performCovarianceTransport(state.stepping.covTransport), nav(state.stepping.navDir) {}
    		
		const Surface* surface;
		const Vector3D pos;
		const double time;
		const Vector3D dir;
		const double momentum;
		const double q;
		const double qOverP;
		const double mass;
		const int pdg;
		const bool performCovarianceTransport;
		const NavigationDirection nav;
		
		MaterialProperties slab;
		double pathCorrection;
		Vector3D variances;
		double Eloss;
		
		double nextP;
	};
	
  void covarianceContributions(Data& d, bool multipleScattering, bool energyLoss)
  {
	// Compute contributions from interactions
	if(multipleScattering)
	{
	  // TODO use momentum before or after energy loss in backward mode?
	  const auto theta0 =
		  computeMultipleScatteringTheta0(d.slab, d.pdg, d.mass, d.qOverP, d.q);
	  // sigmaPhi = theta0 / sin(theta)
	  const auto sigmaPhi = theta0 * (d.dir.norm() / d.dir.z());
	  d.variances.x() = sigmaPhi * sigmaPhi;
	  // sigmaTheta = theta0
	  d.variances.y() = theta0 * theta0;
	}
	// TODO just ionisation loss or full energy loss?
	if(energyLoss) {
	  const auto sigmaQOverP =
		  computeEnergyLossLandauSigmaQOverP(d.slab, d.pdg, d.mass, d.qOverP, d.q);
	  d.variances.z() = sigmaQOverP * sigmaQOverP;
	}
  }
 
	void evaluatePointwiseMaterialInteraction(Data& d, bool multipleScattering, bool energyLoss)
	{
		if (energyLoss) {
		  d.Eloss = computeEnergyLossBethe(d.slab, d.pdg, d.mass, d.qOverP, d.q);
		}
		// Compute contributions from interactions
		if (d.performCovarianceTransport) {
			covarianceContributions(d, multipleScattering, energyLoss);
		}
	}
	
  template<typename propagator_state_t>
  bool evaluateMaterialProperties(const propagator_state_t& state, Data& d,
      MaterialUpdateStage updateStage = fullUpdate) {  
    // We are at the start surface
    if (d.surface == state.navigation.startSurface) {
      updateStage = postUpdate;
      // Or is it the target surface ?
    } else if (d.surface == state.navigation.targetSurface) {
      updateStage = preUpdate;
    }

	// Retrieve the material properties
	d.slab = state.navigation.currentSurface->surfaceMaterial()->materialProperties(d.pos, d.nav, updateStage);

    // Correct the material properties for non-zero incidence
    d.pathCorrection =
        d.surface->pathCorrection(state.geoContext, d.pos, d.dir);
    d.slab.scaleThickness(d.pathCorrection);

    // Get the surface material & properties from them
    return d.slab;
  }

  void changeVariance(double& variance, const double change, const NavigationDirection nav)
  {
	 variance = std::max(0., variance + std::copysign(change, nav));
  }
  
  template<typename propagator_state_t, typename stepper_t>
  void updateState(propagator_state_t& state, const stepper_t& stepper, Data& d)
  {
	// in forward(backward) propagation, energy decreases(increases) and variances increase(decrease)
    const auto nextE = std::sqrt(d.mass * d.mass + d.momentum * d.momentum) - std::copysign(d.Eloss, d.nav);
    // put particle at rest if energy loss is too large
	d.nextP = (d.mass < nextE) ? std::sqrt(nextE * nextE - d.mass * d.mass) : 0;
	// update track parameters and covariance
	stepper.update(state.stepping, d.pos, d.dir, d.nextP, d.time);
	changeVariance(state.stepping.cov(ePHI, ePHI), d.variances.x(), d.nav);
	changeVariance(state.stepping.cov(eTHETA, eTHETA), d.variances.y(), d.nav);
	changeVariance(state.stepping.cov(eQOP, eQOP), d.variances.z(), d.nav);
  }
}  // namespace detail
}  // end of namespace Acts