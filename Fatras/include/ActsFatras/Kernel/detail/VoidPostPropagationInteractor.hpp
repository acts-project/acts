// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsFatras/Kernel/SimulationResult.hpp"
#include "ActsFatras/EventData/Particle.hpp"

namespace ActsFatras {

/// Void class that allows simulating effects that e.g. destroy the particle. Since these types of effects will stop the simulation in any case there is no need to implement these in the physics list used by the @c Interactor.
/// This class allows to abort the propagation based on conditions that are set based on the initial condition (i.e. the particle properties before the propagation). After the propagation is finished, this class can manipulate the resulting particles once more (i.e. by calculating a particle decay).
/// @note Any conditions used in this class have no effect upon the @c Interactor.
struct VoidPostPropagationInteractor {

	/// @brief Aborter condition used in the propagation. This method allows to trigger its effect(s) if its conditions are fulfilled.
	constexpr bool operator()(const Particle::Scalar /*x0*/, const Particle::Scalar /*l0*/, const Particle::Scalar /*time*/) const {
	  return false;
	}
		
	/// @brief Configuration method. Allows to configure its aborting conditions based upon the initial state
	template <typename generator_t>
	void setAbortConditions(generator_t& /*generator*/, const Particle& /*particle*/) {
	}
	 
	/// @brief The actual effect that is triggered after the propagation is finished.
	template <typename generator_t>
	void operator()(generator_t& /*generator*/, SimulationResult& /*result*/) const {
	 }
};
}