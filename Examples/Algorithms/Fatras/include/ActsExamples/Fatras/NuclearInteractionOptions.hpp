// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "Acts/Definitions/Units.hpp"
#include "ActsExamples/Fatras/FatrasAlgorithm.hpp"
#include "ActsExamples/Utilities/OptionsFwd.hpp"
#include "ActsFatras/Kernel/Process.hpp"
#include "ActsFatras/Physics/StandardPhysicsLists.hpp"
#include "ActsFatras/Physics/NuclearInteraction/NuclearInteraction.hpp"
#include "ActsFatras/Physics/NuclearInteraction/Parameters.hpp"

#include <utility>

#include <boost/program_options.hpp>

namespace ActsExamples {
namespace Options {

/// Add Fatras options.
///
/// @param desc The options description to add options to
void addNuclearInteractionOptions(Description& desc);

/// Reads the parametrisation and provides the parametrisation
ActsFatras::detail::MultiParticleParametrisation readParametrisations(
    const std::string& fileName);
      
/// Read Fatras options to create the algorithm config.
///
/// @tparam simulator_t type of the simulation kernel
/// @param vars         the variables to read from
/// @param simulator    the simulation kernel
template <typename simulator_t>
void readNuclearInteractionConfig(
    const boost::program_options::variables_map& variables, simulator_t& simulator) {

  ACTS_LOCAL_LOGGER(
      Acts::getDefaultLogger("NuclearInteractionOptions", Acts::Logging::INFO))
  
  const auto nuclearInteractionParametrisation = variables["fatras-nuclear-interaction-parametrisation"].as<std::string>();

	if(nuclearInteractionParametrisation.empty()) {
		ACTS_WARNING("No parametrisation for the nuclear interaction provided.");
		return;
	}
	ACTS_VERBOSE("Nuclear interaction parametrisation file" << nuclearInteractionParametrisation << " provided.");
	
	auto& chargedNuclearInteraction = simulator.charged.pointlike.template get<ActsFatras::NuclearInteraction>();
	auto& neutralNuclearInteraction = simulator.neutral.pointlike.template get<ActsFatras::NuclearInteraction>();

	const auto mpp = readParametrisations(nuclearInteractionParametrisation);
	ACTS_VERBOSE("Parametrisations for nuclear interaction from " << mpp.size() << " particles provided");
	
	chargedNuclearInteraction.multiParticleParameterisation = mpp;
	neutralNuclearInteraction.multiParticleParameterisation = mpp;
}

}  // namespace Options
}  // namespace ActsExamples