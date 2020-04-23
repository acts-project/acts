// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <boost/program_options.hpp>

#include "ACTFW/Generators/EventGenerator.hpp"

namespace FW {
namespace Options {

/// Add options for a particle-gun-like event generator.
void addParticleGunOptions(boost::program_options::options_description& opt);

/// Create the event generator config from particle-gun options.
EventGenerator::Config readParticleGunOptions(
    const boost::program_options::variables_map& vm);

}  // namespace Options
}  // namespace FW
