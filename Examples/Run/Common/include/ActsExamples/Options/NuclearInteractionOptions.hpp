// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Units.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Utilities/OptionsFwd.hpp"
#include "ActsFatras/Physics/NuclearInteraction/NuclearInteraction.hpp"
#include "ActsFatras/Physics/NuclearInteraction/NuclearInteractionParameters.hpp"

#include <utility>

#include "boost/program_options.hpp"

namespace ActsExamples {
namespace Options {

/// Add Fatras options.
///
/// @param desc The options description to add options to
void addNuclearInteractionOptions(Description& desc);

/// Reads the parametrisation and provides the parametrisation
ActsFatras::detail::MultiParticleNuclearInteractionParametrisation
readParametrisations(const std::string& fileName);

/// Read Fatras options to create the algorithm config.
///
/// @param variables the variables to read from
std::string readNuclearInteractionConfig(
    const boost::program_options::variables_map& variables);

/// Read the parametrisations.
///
/// @tparam simulator_t type of the simulation kernel
/// @param [in] nuclearInteractionParametrisation File name of the
/// parametrisations
/// @param [in, out] simulator The simulation kernel
template <typename simulator_t>
void setNuclearInteractionParametrisations(
    const std::string& nuclearInteractionParametrisation,
    simulator_t& simulator) {
  if (nuclearInteractionParametrisation.empty()) {
    return;
  }

  auto& chargedNuclearInteraction =
      simulator.charged.interactions
          .template get<ActsFatras::NuclearInteraction>();
  auto& neutralNuclearInteraction =
      simulator.neutral.interactions
          .template get<ActsFatras::NuclearInteraction>();

  const auto mpp = readParametrisations(nuclearInteractionParametrisation);

  chargedNuclearInteraction.multiParticleParameterisation = mpp;
  neutralNuclearInteraction.multiParticleParameterisation = mpp;
}

}  // namespace Options
}  // namespace ActsExamples
